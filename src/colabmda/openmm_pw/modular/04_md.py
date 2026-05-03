import os, sys, argparse
from openmm.app import DCDReporter, StateDataReporter, CheckpointReporter, PDBFile
from openmm import XmlSerializer, unit
from colabmda.openmm_pw.modular.utils import make_sim, atomic_rename, safe_remove_if_empty, sync_outputs

def run_md(workdir, pdbid, total_ns=1.0, traj_interval_ps=10.0, checkpoint_ps=1000.0, sync_dir=None):
    os.chdir(workdir)
    
    xml_system = 'system.xml'
    pdb_saved  = 'solvated.pdb'
    chk_in     = 'npt.chk'
    chk_file   = 'prod.chk'

    if not all(os.path.exists(f) for f in [xml_system, pdb_saved]):
        print("Error: system.xml or solvated.pdb not found.")
        return False

    dt = 2.0*unit.femtoseconds
    dt_ps = dt.value_in_unit(unit.picoseconds)
    # Use round to avoid precision issues (e.g. 10.0/0.002 => 4999.999... -> 4999)
    total_steps = int(round((total_ns*unit.nanoseconds)/dt))
    snapshot_steps = int(round((traj_interval_ps*unit.picoseconds)/dt))
    chk_steps = int(round((checkpoint_ps*unit.picoseconds)/dt))

    with open(xml_system) as f:
        system = XmlSerializer.deserializeSystem(f.read())
    pdb = PDBFile(pdb_saved)
    sim = make_sim(pdb.topology, system, dt)

    # Resume from prod.chk if exists, else start fresh from npt.chk
    if os.path.exists(chk_file):
        print("▶ Resuming Production from prod.chk …")
        with open(chk_file, 'rb') as f:
            sim.loadCheckpoint(f)
        steps_done = sim.context.getState().getStepCount()
    else:
        print("▶ Starting Production from npt.chk …")
        if not os.path.exists(chk_in):
            print("Error: npt.chk not found. Run NPT first.")
            return False
        with open(chk_in, 'rb') as f:
            sim.loadCheckpoint(f)
        
        # Reset step count and time so production starts at 0.0
        sim.context.setStepCount(0)
        sim.context.setTime(0.0)
        steps_done = 0

    print(f"  • Progress: {steps_done*dt_ps:.1f} / {total_ns*1000:.1f} ps")

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

        sim.reporters = []
        sim.reporters.append(DCDReporter(dcd_tmp, snapshot_steps))
        sim.reporters.append(StateDataReporter(log_tmp, snapshot_steps,
                                               step=True, time=True,
                                               potentialEnergy=True, temperature=True,
                                               volume=True, speed=True))
        sim.reporters.append(CheckpointReporter(chk_file, chk_steps))

        try:
            sim.step(steps_to_run)
            interrupted = False
        except KeyboardInterrupt:
            interrupted = True
            print("⚠ Interrupted")
            break
        except Exception as e:
            print(f"⚠ Step failed: {e}")
            return False

        # Completion handling
        sim.saveCheckpoint(chk_file)
        new_steps = sim.context.getState().getStepCount()
        
        if not interrupted:
            atomic_rename(dcd_tmp, dcd_fin)
            atomic_rename(log_tmp, log_fin)
            sync_outputs(workdir, sync_dir, extra_files=[dcd_fin, log_fin])
        
        steps_done = new_steps
        print(f"  ✔ Reached step {steps_done}/{total_steps}")

    if steps_done >= total_steps:
        print("✔ Full production complete")
    return True

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("workdir")
    p.add_argument("--pdbid", required=True)
    p.add_argument("--total-ns", type=float, default=1.0)
    p.add_argument("--traj-interval", type=float, default=10.0)
    p.add_argument("--checkpoint-ps", type=float, default=1000.0)
    p.add_argument("--sync-dir", default=None)
    args = p.parse_args()
    run_md(args.workdir, args.pdbid, args.total_ns, args.traj_interval, args.checkpoint_ps, args.sync_dir)
