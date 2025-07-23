#!/usr/bin/env python3
import os, shutil, time
from pdbfixer import PDBFixer
from openmm.app import (
    Modeller, ForceField, Simulation,
    DCDReporter, StateDataReporter, CheckpointReporter,
    PME, HBonds, PDBFile
)
from openmm import XmlSerializer, unit, Platform, MonteCarloBarostat, LangevinMiddleIntegrator

# ----- User‐tweakable parameters -----
total_ns          = 5.0    # total simulation length in ns
snapshot_interval = 100.0  # snapshot every X ps
equil_ps          = 100.0  # equilibration length in ps
# ----------------------------------

# Paths & filenames
workdir    = '/content/drive/MyDrive/openmm/1aki'
system_xml = os.path.join(workdir, 'system.xml')
eq_state   = os.path.join(workdir, 'eq_state.xml')
chk_file   = os.path.join(workdir, 'prod.chk')
pdb_saved  = os.path.join(workdir, 'solvated.pdb')

# Derived quantities
dt               = 2.0 * unit.femtoseconds
steps_per_ps     = int((1*unit.picoseconds) / dt)
steps_per_snap   = int((snapshot_interval*unit.picoseconds) / dt)
total_steps      = int((total_ns*unit.nanoseconds) / dt)
n_eq_steps       = int((equil_ps*unit.picoseconds) / dt)
checkpoint_steps = steps_per_snap * 10

# Fresh-run guard: only delete if no eq_state and no checkpoint
if not (os.path.isfile(eq_state) or os.path.isfile(chk_file)):
    if os.path.isdir(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
else:
    print("eq_state or checkpoint found; preserving directory")

os.chdir(workdir)
print("Working directory:", workdir)

# Decide whether to rebuild or load
if os.path.isfile(eq_state) and os.path.isfile(system_xml) and os.path.isfile(pdb_saved):
    # --- Resume path ---
    print("Loading prepared system and topology...")
    with open(system_xml) as f:
        system = XmlSerializer.deserializeSystem(f.read())
    pdb = PDBFile('solvated.pdb')
    integrator = LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, dt)
    platform   = Platform.getPlatformByName('CUDA')
    sim = Simulation(pdb.topology, system, integrator, platform)
    sim.loadState(eq_state)
    steps_done = sim.context.getState().getStepCount()
    print(f"Resuming production at {steps_done} steps ({(steps_done*dt)/unit.nanoseconds:.3f} ns)")
else:
    # --- Fresh setup ---
    # Download PDB
    if not os.path.exists('1AKI.pdb'):
        os.system('wget -q -O 1AKI.pdb https://files.rcsb.org/download/1AKI.pdb')
    # Fix PDB
    fixer = PDBFixer('1AKI.pdb')
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.0)
    modeller = Modeller(fixer.topology, fixer.positions)
    ff = ForceField('amber14-all.xml', 'amber14/tip3p.xml')
    modeller.addSolvent(ff, model='tip3p', padding=1.0*unit.nanometer, ionicStrength=0.15*unit.molar)
    # Build system
    system = ff.createSystem(modeller.topology, nonbondedMethod=PME,
                             nonbondedCutoff=1.0*unit.nanometer, constraints=HBonds)
    system.addForce(MonteCarloBarostat(1*unit.atmosphere, 300*unit.kelvin, 25))
    # Simulation setup
    integrator = LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, dt)
    platform   = Platform.getPlatformByName('CUDA')
    sim = Simulation(modeller.topology, system, integrator, platform)
    sim.context.setPositions(modeller.positions)
    # Minimize & equilibrate
    print("Minimizing…"); sim.minimizeEnergy()
    print("Equilibrating NVT…")
    sim.context.setVelocitiesToTemperature(300*unit.kelvin)
    sim.reporters.append(StateDataReporter('nvt.log', max(1, n_eq_steps//10), step=True, temperature=True))
    sim.step(n_eq_steps); sim.reporters.clear()
    print("Equilibrating NPT…")
    sim.reporters.append(StateDataReporter('npt.log', max(1, n_eq_steps//10), step=True, temperature=True, volume=True))
    sim.step(n_eq_steps); sim.reporters.clear()
    print("✔︎ Equilibration complete")
    # Serialize system + state + PDB
    with open(system_xml, 'w') as f:
        f.write(XmlSerializer.serializeSystem(sim.context.getSystem()))
    sim.saveState(eq_state)
    PDBFile.writeFile(sim.topology,
                      sim.context.getState(getPositions=True).getPositions(),
                      open('solvated.pdb','w'))
    steps_done = 0

# --- Production ---
print(f"Running production from {steps_done} to {total_steps} steps")
dcd = f'prod_{int(total_ns*1000)}ps.dcd'
log = f'prod_{int(total_ns*1000)}ps.log'
sim.reporters.append(DCDReporter(dcd, steps_per_snap))
sim.reporters.append(StateDataReporter(log, steps_per_snap,
    step=True, time=True, potentialEnergy=True,
    temperature=True, volume=True, speed=True))
sim.reporters.append(CheckpointReporter(chk_file, checkpoint_steps))

while steps_done < total_steps:
    to_run = min(checkpoint_steps, total_steps - steps_done)
    sim.step(to_run)
    steps_done += to_run
    sim.saveCheckpoint(chk_file)
    print(f"  ✔︎ Checkpoint at {(steps_done*dt)/unit.nanoseconds:.3f} ns")

print("✔︎ Production complete")
print("Outputs:", os.listdir('.'))
