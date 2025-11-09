#!/usr/bin/env python3
"""
openmm_proteinwater.py

Runs minimization, equilibration, and production MD directly in the
same directory where <pdbid>_cleaned.pdb lives.

Usage:
    python3 openmm_proteinwater.py <pdbid> [options]

Positional arguments:
  pdbid                4-letter PDB identifier (e.g. 4ldj)

Optional arguments:
  -h, --help           show this help message and exit
  --total-ns FLOAT     total production time in ns (default: 1.0)
  --traj-interval FLOAT
                       trajectory write interval in ps (default: 10.0)
  --equil-time FLOAT   each equilibration phase duration in ps (default: 10.0)
  --checkpoint-ps FLOAT
                       checkpoint interval in ps (default: 10.0)
"""

import os, sys
import argparse
from pdbfixer import PDBFixer
from openmm.app import (
    Modeller, ForceField, Simulation,
    DCDReporter, StateDataReporter, CheckpointReporter,
    PME, HBonds, PDBFile
)
from openmm import XmlSerializer, unit, Platform, MonteCarloBarostat, LangevinMiddleIntegrator

def main():
    p = argparse.ArgumentParser(
        description="Run MD in-place in the <pdbid>/ folder",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument("pdbid", help="4-letter PDB identifier (e.g. 4ldj)")
    p.add_argument("--total-ns",      type=float, default=1.0,  help="Total production time (ns)")
    p.add_argument("--traj-interval", type=float, default=10.0, help="Trajectory interval (ps)")
    p.add_argument("--equil-time",    type=float, default=10.0, help="Equilibration phase length (ps)")
    p.add_argument("--checkpoint-ps", type=float, default=10.0, help="Checkpoint interval (ps)")
    args = p.parse_args()

    pdb_id = args.pdbid.lower()
    workdir = os.path.join(os.getcwd(), pdb_id)
    cleaned_pdb = os.path.join(workdir, f"{pdb_id}_cleaned.pdb")
    if not os.path.exists(cleaned_pdb):
        print(f"Error: cleaned PDB not found: {cleaned_pdb}", file=sys.stderr)
        sys.exit(1)

    # Enter the folder
    os.chdir(workdir)
    print(f"→ Running MD in: {workdir}")

    # User‐tweakable (now from args)
    total_ns      = args.total_ns
    traj_interval = args.traj_interval
    equil_time    = args.equil_time
    checkpoint_ps = args.checkpoint_ps

    # Derived
    dt               = 2.0 * unit.femtoseconds
    steps_per_ps     = int((1*unit.picoseconds) / dt)
    steps_snapshot   = int((traj_interval*unit.picoseconds) / dt)
    n_eq_steps       = int((equil_time*unit.picoseconds) / dt)
    checkpoint_steps = int(checkpoint_ps * steps_per_ps)
    total_steps      = int((total_ns*unit.nanoseconds) / dt)

    def make_sim(topology, system):
        integ = LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, dt)
        plat  = Platform.getPlatformByName('CUDA')
        return Simulation(topology, system, integ, plat)

    # Filenames
    xml_system = 'system.xml'
    pdb_saved  = 'solvated.pdb'
    chk_file   = 'prod.chk'

    # Resume?
    if os.path.isfile(chk_file) and os.path.isfile(xml_system) and os.path.isfile(pdb_saved):
        print("→ Resuming from checkpoint")
        with open(xml_system) as f:
            system = XmlSerializer.deserializeSystem(f.read())
        pdb = PDBFile(pdb_saved)
        sim = make_sim(pdb.topology, system)
        with open(chk_file,'rb') as f:
            sim.loadCheckpoint(f)
        state = sim.context.getState()
        steps_done = state.getStepCount()
        print(f"Resumed at {steps_done} steps ({(steps_done*dt)/unit.nanoseconds:.3f} ns)")
    else:
        print("→ Fresh setup + equilibration")
        fixer = PDBFixer(filename=cleaned_pdb)
        fixer.removeHeterogens(keepWater=True)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(pH=7.0)

        modeller = Modeller(fixer.topology, fixer.positions)
        ff = ForceField('amber14-all.xml','amber14/tip3p.xml')
        modeller.addSolvent(ff, model='tip3p',
                           padding=1.0*unit.nanometer,
                           ionicStrength=0.15*unit.molar)

        system = ff.createSystem(modeller.topology,
                                 nonbondedMethod=PME,
                                 nonbondedCutoff=1.0*unit.nanometer,
                                 constraints=HBonds)
        system.addForce(MonteCarloBarostat(1*unit.atmosphere,300*unit.kelvin,25))

        sim = make_sim(modeller.topology, system)
        sim.context.setPositions(modeller.positions)

        print("Minimizing…"); sim.minimizeEnergy()

        print("Equil NVT…")
        sim.context.setVelocitiesToTemperature(300*unit.kelvin)
        sim.reporters.append(StateDataReporter('nvt.log', n_eq_steps//10,
                                               step=True, temperature=True))
        sim.step(n_eq_steps); sim.reporters.clear()

        print("Equil NPT…")
        sim.reporters.append(StateDataReporter('npt.log', n_eq_steps//10,
                                               step=True, temperature=True, volume=True))
        sim.step(n_eq_steps); sim.reporters.clear()

        with open(xml_system,'w') as f:
            f.write(XmlSerializer.serializeSystem(sim.context.getSystem()))
        PDBFile.writeFile(sim.topology,
                          sim.context.getState(getPositions=True).getPositions(),
                          open(pdb_saved,'w'))
        sim.saveCheckpoint(chk_file)
        steps_done = 0

    # Production
    print(f"Starting production: {steps_done}/{total_steps} steps")
    sim.reporters.append(DCDReporter(f'prod_{int(total_ns*1000)}ps.dcd', steps_snapshot))
    sim.reporters.append(StateDataReporter(f'prod_{int(total_ns*1000)}ps.log',
                                           steps_snapshot,
                                           step=True, time=True,
                                           potentialEnergy=True,
                                           temperature=True, volume=True,
                                           speed=True))
    sim.reporters.append(CheckpointReporter(chk_file, checkpoint_steps))

    while steps_done < total_steps:
        to_run = min(checkpoint_steps, total_steps-steps_done)
        sim.step(to_run)
        steps_done += to_run
        print(f"  ↪︎ {(steps_done*dt)/unit.nanoseconds:.3f} ns — checkpoint saved", flush=True)

    print("✔︎ Production complete")
    print("Final files:", os.listdir('.'))

if __name__ == "__main__":
    main()
