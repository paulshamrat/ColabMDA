#!/usr/bin/env python3
"""
openmm_proteinwater.py

Automatically runs successive production chunks of length --checkpoint-ps
until --total-ns is reached. Safe to interrupt and resume.
"""

import os, sys, argparse
from pdbfixer import PDBFixer
from openmm.app import (Modeller, ForceField, Simulation,
                        DCDReporter, StateDataReporter, CheckpointReporter,
                        PME, HBonds, PDBFile)
from openmm import XmlSerializer, unit, Platform, MonteCarloBarostat, LangevinMiddleIntegrator

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("pdbid")
    p.add_argument("--total-ns",      type=float, default=1.0)
    p.add_argument("--traj-interval", type=float, default=1.0)
    p.add_argument("--equil-time",    type=float, default=10.0)
    p.add_argument("--checkpoint-ps", type=float, default=10.0)
    return p.parse_args()

def make_sim(top, sys_, dt):
    integ = LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, dt)
    plat  = Platform.getPlatformByName('CUDA')
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
    sim.reporters.append(StateDataReporter('nvt.log', eq_steps//10, step=True, temperature=True))
    sim.step(eq_steps); sim.reporters.clear()

    print("  • Equil NPT…")
    sim.reporters.append(StateDataReporter('npt.log', eq_steps//10, step=True, temperature=True, volume=True))
    sim.step(eq_steps); sim.reporters.clear()

    with open(xml_sys,'w') as f:
        f.write(XmlSerializer.serializeSystem(sim.context.getSystem()))
    PDBFile.writeFile(sim.topology,
                      sim.context.getState(getPositions=True).getPositions(),
                      open(pdb_saved,'w'))
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

def main():
    args = parse_args()
    pdbid = args.pdbid.lower()
    workdir = os.path.join(os.getcwd(), pdbid)
    cleaned_pdb = os.path.join(workdir, f"{pdbid}_cleaned.pdb")
    if not os.path.exists(cleaned_pdb):
        sys.exit("Error: cleaned PDB not found")
    os.chdir(workdir)

    # MD setup
    dt      = 2.0*unit.femtoseconds
    dt_ps   = dt.value_in_unit(unit.picoseconds)
    total_steps      = int((args.total_ns*unit.nanoseconds)/dt)
    snapshot_steps   = int((args.traj_interval*unit.picoseconds)/dt)
    eq_steps         = int((args.equil_time*unit.picoseconds)/dt)
    chk_steps        = int((args.checkpoint_ps*unit.picoseconds)/dt)

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

        # attach reporters
        sim.reporters = []
        sim.reporters.append(DCDReporter(f'prod_{tag}.dcd', snapshot_steps))
        sim.reporters.append(StateDataReporter(f'prod_{tag}.log', snapshot_steps,
                                               step=True, time=True,
                                               potentialEnergy=True, temperature=True,
                                               volume=True, speed=True))
        sim.reporters.append(CheckpointReporter(chk_file, chk_steps))

        try:
            sim.step(steps_to_run)
        except KeyboardInterrupt:
            print("⚠ Interrupted — checkpoint saved")
        finally:
            sim.saveCheckpoint(chk_file)

        # **CRUCIAL** update steps_done
        steps_done += steps_to_run

        print(f"  ✔ Completed {tag}")

    print("✔︎ Full production complete")

if __name__ == "__main__":
    main()
