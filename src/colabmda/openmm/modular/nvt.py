import argparse
import os
import random

from openmm import XmlSerializer, unit
from openmm.app import PDBFile, StateDataReporter

from colabmda.openmm.modular.utils import make_sim


def run_nvt(workdir, pdbid, equil_time_ps=100.0, seed=None):
    workdir = os.path.abspath(workdir)
    os.chdir(workdir)

    xml_system = os.path.join(workdir, "system.xml")
    pdb_saved = os.path.join(workdir, "solvated.pdb")
    chk_in = os.path.join(workdir, "em.chk")
    chk_out = os.path.join(workdir, "nvt.chk")
    log_file = os.path.join(workdir, "nvt.log")

    print(f"[DEBUG] Checking for NVT files in: {workdir}")
    print(f"[DEBUG] xml: {xml_system}, pdb: {pdb_saved}, in: {chk_in}, out: {chk_out}")
    if os.path.exists(chk_out):
        print(f"✔ NVT already complete ({chk_out} found). Skipping...")
        return True

    if not all(os.path.exists(f) for f in [xml_system, pdb_saved, chk_in]):
        print("Error: Missing em.chk or system files. Run EM first.")
        return False

    print(f"▶ Starting NVT Equilibration ({equil_time_ps} ps) …")
    with open(xml_system) as f:
        system = XmlSerializer.deserializeSystem(f.read())
    pdb = PDBFile(pdb_saved)

    dt = 2.0 * unit.femtoseconds
    sim = make_sim(pdb.topology, system, dt)

    with open(chk_in, "rb") as f:
        sim.loadCheckpoint(f)

    # Assign velocities
    if seed is None:
        seed = random.randint(1, 1000000)
    print(f"  • Setting velocities (seed={seed})")
    sim.context.setVelocitiesToTemperature(300 * unit.kelvin, seed)

    eq_steps = max(1, int((equil_time_ps * unit.picoseconds) / dt))
    sim.reporters.append(
        StateDataReporter(
            log_file,
            max(1, eq_steps // 50),
            step=True,
            time=True,
            temperature=True,
            potentialEnergy=True,
        )
    )

    sim.step(eq_steps)
    sim.saveCheckpoint(chk_out)
    print("✔ NVT complete (nvt.chk saved)")
    return True


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("workdir")
    p.add_argument("--pdbid", required=True)
    p.add_argument("--equil-time", type=float, default=100.0)
    p.add_argument("--seed", type=int, default=None)
    args = p.parse_args()
    run_nvt(args.workdir, args.pdbid, args.equil_time, args.seed)
