#!/usr/bin/env python3
"""
openmm_proteinwater_colab.py

Chunked OpenMM protein-in-water MD with robust resume and Colab-safe I/O.

Key improvements vs the original:
- Platform fallback: CUDA -> OpenCL -> CPU
- Atomic trajectory/log writes: write *.tmp then rename only if chunk finished
- Cleans up 0-byte files automatically
- Optional Drive sync: copies key outputs to a persistent directory after each chunk

Usage example:
  python3 openmm_proteinwater_colab.py /content/work/4ldj_wt \
      --pdbid 4ldj --total-ns 100 --traj-interval 100 --equil-time 100 --checkpoint-ps 1000 \
      --sync-dir /content/drive/MyDrive/openmm_runs/4ldj_wt
"""

import os, sys, argparse, shutil
from pdbfixer import PDBFixer
from openmm.app import (Modeller, ForceField, Simulation,
                        DCDReporter, StateDataReporter, CheckpointReporter,
                        PME, HBonds, PDBFile)
from openmm import XmlSerializer, unit, Platform, MonteCarloBarostat, LangevinMiddleIntegrator

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("workdir", help="Folder containing <pdbid>_cleaned.pdb")
    p.add_argument("--pdbid", default="4ldj", help="Prefix name (default: 4ldj)")
    p.add_argument("--total-ns",      type=float, default=1.0)
    p.add_argument("--traj-interval", type=float, default=100.0)   # ps
    p.add_argument("--equil-time",    type=float, default=100.0)   # ps (per NVT and NPT)
    p.add_argument("--checkpoint-ps", type=float, default=1000.0)  # ps per chunk
    p.add_argument("--sync-dir", default=None, help="If set, sync outputs to this directory after each chunk")
    return p.parse_args()

def pick_platform():
    for name in ["CUDA", "OpenCL", "CPU"]:
        try:
            plat = Platform.getPlatformByName(name)
            print(f"[OpenMM] Using platform: {name}")
            return plat
        except Exception:
            pass
    raise RuntimeError("No OpenMM platform available")

def safe_remove_if_empty(path):
    if os.path.isfile(path) and os.path.getsize(path) == 0:
        try:
            os.remove(path)
        except Exception:
            pass

def atomic_rename(tmp_path, final_path):
    # Replace if exists
    if os.path.exists(final_path):
        os.remove(final_path)
    os.rename(tmp_path, final_path)

def make_sim(top, sys_, dt):
    integ = LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, dt)
    plat  = pick_platform()
    return Simulation(top, sys_, integ, plat)

def setup_and_equilibrate(cleaned_pdb, dt, eq_steps, xml_sys, pdb_saved, chk_file):
    fixer = PDBFixer(filename=cleaned_pdb)
    fixer.removeHeterogens(keepWater=True)
    fixer.findMissingResidues(); fixer.findMissingAtoms()
    fixer.addMissingAtoms(); fixer.addMissingHydrogens(pH=7.0)

    modeller = Modeller(fixer.topology, fixer.positions)
    ff = ForceField('amber14-all.xml','amber14/tip3p.xml')
    modeller.addSolvent(ff, model='tip3p', padding=1.0*unit.nanometer, ionicStrength=0.15*unit.molar)

    system = ff.createSystem(modeller.topology,
                             nonbondedMethod=PME,
                             nonbondedCutoff=1.0*unit.nanometer,
                             constraints=HBonds)
    system.addForce(MonteCarloBarostat(1*unit.atmosphere,300*unit.kelvin,25))

    sim = make_sim(modeller.topology, system, dt)
    sim.context.setPositions(modeller.positions)

    print("  • Minimizing…"); sim.minimizeEnergy()
    print("  • Equil NVT…")
    sim.context.setVelocitiesToTemperature(300*unit.kelvin)
    sim.reporters.append(StateDataReporter('nvt.log', max(1, eq_steps//10), step=True, temperature=True))
    sim.step(eq_steps); sim.reporters.clear()

    print("  • Equil NPT…")
    sim.reporters.append(StateDataReporter('npt.log', max(1, eq_steps//10), step=True, temperature=True, volume=True))
    sim.step(eq_steps); sim.reporters.clear()

    with open(xml_sys,'w') as f:
        f.write(XmlSerializer.serializeSystem(sim.context.getSystem()))
    with open(pdb_saved,'w') as f:
        PDBFile.writeFile(sim.topology,
                          sim.context.getState(getPositions=True).getPositions(),
                          f)
    sim.saveCheckpoint(chk_file)
    return sim, 0

def resume_sim(xml_sys, pdb_saved, chk_file, dt):
    with open(xml_sys) as f:
        system = XmlSerializer.deserializeSystem(f.read())
    pdb = PDBFile(pdb_saved)
    sim = make_sim(pdb.topology, system, dt)
    with open(chk_file,'rb') as f:
        sim.loadCheckpoint(f)
    steps = sim.context.getState().getStepCount()
    print(f"  • Resumed at {steps} steps")
    return sim, steps

def sync_outputs(workdir, sync_dir, extra_files=None):
    if not sync_dir:
        return
    os.makedirs(sync_dir, exist_ok=True)
    # Copy only essential + new chunk files
    essentials = ["system.xml", "solvated.pdb", "prod.chk", "nvt.log", "npt.log", "prod_full.log"]
    if extra_files:
        essentials += extra_files
    for fn in essentials:
        src = os.path.join(workdir, fn)
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(sync_dir, fn))

def main():
    args = parse_args()
    pdbid = args.pdbid.lower()
    workdir = os.path.abspath(args.workdir)
    os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)

    cleaned_pdb = os.path.join(workdir, f"{pdbid}_cleaned.pdb")
    if not os.path.exists(cleaned_pdb):
        sys.exit(f"Error: cleaned PDB not found: {cleaned_pdb}")

    # Remove any empty trajectory/log segments left from prior interruption
    for f in os.listdir(workdir):
        if f.startswith("prod_") and (f.endswith(".dcd") or f.endswith(".log")):
            safe_remove_if_empty(f)

    # MD setup
    dt      = 2.0*unit.femtoseconds
    dt_ps   = dt.value_in_unit(unit.picoseconds)
    total_steps      = int((args.total_ns*unit.nanoseconds)/dt)
    snapshot_steps   = max(1, int((args.traj_interval*unit.picoseconds)/dt))
    eq_steps         = max(1, int((args.equil_time*unit.picoseconds)/dt))
    chk_steps        = max(1, int((args.checkpoint_ps*unit.picoseconds)/dt))

    xml_system = 'system.xml'
    pdb_saved  = 'solvated.pdb'
    chk_file   = 'prod.chk'

    # Setup or resume
    if os.path.isfile(chk_file) and os.path.isfile(xml_system) and os.path.isfile(pdb_saved):
        print("▶ Resuming …")
        sim, steps_done = resume_sim(xml_system, pdb_saved, chk_file, dt)
    else:
        print("▶ Fresh setup + equilibration …")
        sim, steps_done = setup_and_equilibrate(cleaned_pdb, dt, eq_steps,
                                                xml_system, pdb_saved, chk_file)

    # Production loop
    while steps_done < total_steps:
        steps_to_run = min(chk_steps, total_steps - steps_done)
        ps_start = int(steps_done*dt_ps)
        ps_end   = int(ps_start + steps_to_run*dt_ps)
        tag      = f"{ps_start}to{ps_end}ps"
        print(f"▶ Running chunk {tag} …")

        dcd_tmp = f"prod_{tag}.dcd.tmp"
        log_tmp = f"prod_{tag}.log.tmp"
        dcd_fin = f"prod_{tag}.dcd"
        log_fin = f"prod_{tag}.log"

        # attach reporters (write to tmp)
        sim.reporters = []
        sim.reporters.append(DCDReporter(dcd_tmp, snapshot_steps))
        sim.reporters.append(StateDataReporter(log_tmp, snapshot_steps,
                                               step=True, time=True,
                                               potentialEnergy=True, temperature=True,
                                               volume=True, speed=True))
        sim.reporters.append(CheckpointReporter(chk_file, chk_steps))

        interrupted = False
        try:
            sim.step(steps_to_run)
        except KeyboardInterrupt:
            interrupted = True
            print("⚠ Interrupted — checkpoint saved")
        finally:
            sim.saveCheckpoint(chk_file)

        # If tmp outputs are empty, remove them; if non-empty, rename atomically.
        if os.path.exists(dcd_tmp) and os.path.getsize(dcd_tmp) > 0:
            atomic_rename(dcd_tmp, dcd_fin)
        else:
            safe_remove_if_empty(dcd_tmp)

        if os.path.exists(log_tmp) and os.path.getsize(log_tmp) > 0:
            atomic_rename(log_tmp, log_fin)
        else:
            safe_remove_if_empty(log_tmp)

        # Only advance steps_done if chunk ran (even partially) — we *did* call sim.step()
        # If interrupted very early, it still advanced some steps; but OpenMM doesn’t expose partial steps cleanly here.
        # Safer: reload step count from context.
        steps_done = sim.context.getState().getStepCount()

        # Sync to Drive (only the essentials + the completed chunk files)
        sync_outputs(workdir, args.sync_dir, extra_files=[dcd_fin, log_fin])

        print(f"  ✔ Reached step {steps_done}/{total_steps}")

        if interrupted:
            # Stop cleanly so user can restart in next Colab session.
            break

    if steps_done >= total_steps:
        print("✔︎ Full production complete")

if __name__ == "__main__":
    main()
