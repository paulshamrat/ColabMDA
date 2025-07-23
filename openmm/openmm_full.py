#!/usr/bin/env python3
import os
import shutil
import time
from pdbfixer import PDBFixer
from openmm.app import (
    Modeller, ForceField, Simulation,
    DCDReporter, StateDataReporter, CheckpointReporter,
    PME, HBonds, PDBFile
)
from openmm import unit, Platform, MonteCarloBarostat, LangevinMiddleIntegrator

# ----- User‐tweakable parameters -----
total_ns          = 2.0    # total simulation length in ns
snapshot_interval = 10.0   # snapshot every X ps
equil_ps          = 10.0   # equilibration length in ps for both NVT and NPT
# ------------------------------------

# Derived parameters
dt               = 2.0 * unit.femtoseconds
steps_per_ps     = int((1*unit.picoseconds) / dt)
steps_per_snap   = int((snapshot_interval*unit.picoseconds) / dt)
total_steps      = int((total_ns*unit.nanoseconds) / dt)
n_eq_steps       = int((equil_ps*unit.picoseconds) / dt)
checkpoint_steps = steps_per_snap * 10  # checkpoint every 10 snapshots

# 1) Determine script directory and create fresh "1aki" subfolder
script_dir = os.path.dirname(os.path.abspath(__file__))
workdir    = os.path.join(script_dir, '1aki')
if os.path.isdir(workdir):
    shutil.rmtree(workdir)
os.makedirs(workdir)
os.chdir(workdir)
print("Working directory:", os.getcwd())

# 2) Download PDB if missing
if not os.path.exists('1AKI.pdb'):
    os.system('wget -q -O 1AKI.pdb https://files.rcsb.org/download/1AKI.pdb')

# 3) Load & fix PDB
fixer = PDBFixer('1AKI.pdb')
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(pH=7.0)

# 4) Solvate & ionize
modeller = Modeller(fixer.topology, fixer.positions)
ff       = ForceField('amber14-all.xml', 'amber14/tip3p.xml')
modeller.addSolvent(ff,
                    model='tip3p',
                    padding=1.0*unit.nanometer,
                    ionicStrength=0.15*unit.molar)

# 5) Build system + barostat
system = ff.createSystem(modeller.topology,
                         nonbondedMethod=PME,
                         nonbondedCutoff=1.0*unit.nanometer,
                         constraints=HBonds)
system.addForce(MonteCarloBarostat(1*unit.atmosphere, 300*unit.kelvin, 25))

# 6) Integrator & platform
integrator = LangevinMiddleIntegrator(300*unit.kelvin,
                                      1/unit.picosecond,
                                      dt)
platform   = Platform.getPlatformByName('CUDA')
platform.setPropertyDefaultValue('DeviceIndex','0')
platform.setPropertyDefaultValue('Precision','mixed')

# 7) Create Simulation & minimize
sim = Simulation(modeller.topology, system, integrator, platform)
sim.context.setPositions(modeller.positions)
print("Minimizing energy…")
sim.minimizeEnergy()
print("✔︎ Minimization complete")

# 8) Equilibration: NVT then NPT
sim.context.setVelocitiesToTemperature(300*unit.kelvin)
sim.reporters.append(StateDataReporter('nvt.log', max(1, n_eq_steps//10),
                                       step=True, temperature=True))
sim.step(n_eq_steps)
sim.reporters.clear()
sim.reporters.append(StateDataReporter('npt.log', max(1, n_eq_steps//10),
                                       step=True, temperature=True, volume=True))
sim.step(n_eq_steps)
sim.reporters.clear()
print("✔︎ Equilibration complete")

# 9) Production: chunked run with DCD, log, and checkpoint
dcd_file   = f'prod_{int(total_ns*1000)}ps.dcd'
log_file   = f'prod_{int(total_ns*1000)}ps.log'
chk_file   = 'prod.chk'

sim.reporters.append(DCDReporter(dcd_file, steps_per_snap))
sim.reporters.append(StateDataReporter(log_file,
    steps_per_snap,
    step=True, time=True, potentialEnergy=True,
    temperature=True, volume=True, speed=True
))
sim.reporters.append(CheckpointReporter(chk_file, checkpoint_steps))

print(f"Running production: {total_steps} steps (~{total_ns} ns) in chunks of {checkpoint_steps} steps")

steps_done = 0
while steps_done < total_steps:
    to_run = min(checkpoint_steps, total_steps - steps_done)
    sim.step(to_run)
    steps_done += to_run
    # also save a final checkpoint after each chunk
    sim.saveCheckpoint(chk_file)
    completed_ns = (steps_done * dt) / unit.nanoseconds
    print(f"  ✔︎ Checkpoint at {completed_ns:.3f} ns → {chk_file}")

print("✔︎ Production complete")

# 10) Save final solvated PDB
state = sim.context.getState(getPositions=True)
with open('solvated.pdb', 'w') as f:
    PDBFile.writeFile(sim.topology, state.getPositions(), f)

# 11) List outputs
print("Outputs:", os.listdir('.'))
