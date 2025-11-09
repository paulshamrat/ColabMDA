#!/usr/bin/env python3
import os, shutil
from pdbfixer import PDBFixer
from openmm.app import (
    Modeller, ForceField, Simulation,
    DCDReporter, StateDataReporter, CheckpointReporter,
    PME, HBonds, PDBFile
)
from openmm import XmlSerializer, unit, Platform, MonteCarloBarostat, LangevinMiddleIntegrator

# ----- User‐tweakable parameters -----
total_ns          = 1.0     # total production length in ns
traj_interval_ps  = 10.0   # write DCD every X ps
equil_ps          = 10.0   # equilibration length in ps
checkpoint_ps     = 10.0   # write checkpoint every X ps (same as traj)
# ----------------------------------

# Paths & filenames
workdir    = '/content/drive/MyDrive/openmm/1aki'
system_xml = os.path.join(workdir, 'system.xml')
pdb_saved  = os.path.join(workdir, 'solvated.pdb')
chk_file   = os.path.join(workdir, 'prod.chk')

# Derived quantities
dt               = 2.0 * unit.femtoseconds
steps_per_ps     = int((1*unit.picoseconds) / dt)
steps_per_snap   = int((traj_interval_ps*unit.picoseconds) / dt)
n_eq_steps       = int((equil_ps*unit.picoseconds) / dt)
checkpoint_steps = int(checkpoint_ps * steps_per_ps)
total_steps      = int((total_ns*unit.nanoseconds) / dt)

# Ensure workdir exists (unless resuming)
if not os.path.isdir(workdir):
    os.makedirs(workdir)

os.chdir(workdir)
print("Working directory:", workdir)

# Platform & integrator factory
def make_simulation(topology, system):
    integrator = LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, dt)
    platform   = Platform.getPlatformByName('CUDA')
    sim = Simulation(topology, system, integrator, platform)
    return sim

# Resume from checkpoint if available
if os.path.isfile(chk_file) and os.path.isfile(system_xml) and os.path.isfile(pdb_saved):
    print("→ Found checkpoint; resuming production!")
    # Reconstruct System & Simulation exactly as before
    with open(system_xml) as f:
        system = XmlSerializer.deserializeSystem(f.read())
    pdb = PDBFile('solvated.pdb')
    sim = make_simulation(pdb.topology, system)
    # load binary checkpoint (restores positions, velocities, integrator state)
    with open(chk_file, 'rb') as f:
        sim.loadCheckpoint(f)
    state = sim.context.getState()
    steps_done = state.getStepCount()
    print(f"Resuming at step {steps_done} ({(steps_done*dt)/unit.nanoseconds:.3f} ns)")
else:
    # Fresh setup + equilibration
    print("→ No checkpoint: doing full setup + equilibration")
    if not os.path.exists('1AKI.pdb'):
        os.system('wget -q -O 1AKI.pdb https://files.rcsb.org/download/1AKI.pdb')
    fixer = PDBFixer('1AKI.pdb')
    fixer.findMissingResidues(); fixer.findMissingAtoms()
    fixer.addMissingAtoms(); fixer.addMissingHydrogens(pH=7.0)
    modeller = Modeller(fixer.topology, fixer.positions)
    ff = ForceField('amber14-all.xml', 'amber14/tip3p.xml')
    modeller.addSolvent(ff, model='tip3p',
                       padding=1.0*unit.nanometer, ionicStrength=0.15*unit.molar)
    system = ff.createSystem(modeller.topology, nonbondedMethod=PME,
                             nonbondedCutoff=1.0*unit.nanometer, constraints=HBonds)
    system.addForce(MonteCarloBarostat(1*unit.atmosphere, 300*unit.kelvin, 25))
    sim = make_simulation(modeller.topology, system)
    sim.context.setPositions(modeller.positions)

    print("Minimizing…"); sim.minimizeEnergy()
    print("Equilibrating NVT…")
    sim.context.setVelocitiesToTemperature(300*unit.kelvin)
    sim.reporters.append(StateDataReporter('nvt.log', n_eq_steps//10, step=True, temperature=True))
    sim.step(n_eq_steps); sim.reporters.clear()

    print("Equilibrating NPT…")
    sim.reporters.append(StateDataReporter('npt.log', n_eq_steps//10,
                                          step=True, temperature=True, volume=True))
    sim.step(n_eq_steps); sim.reporters.clear()
    print("✔︎ Equilibration complete")

    # Save system topology and coordinates
    with open(system_xml, 'w') as f:
        f.write(XmlSerializer.serializeSystem(sim.context.getSystem()))
    PDBFile.writeFile(sim.topology,
                      sim.context.getState(getPositions=True).getPositions(),
                      open('solvated.pdb','w'))

    # Write an initial checkpoint so even equilibration can be resumed
    sim.saveCheckpoint(chk_file)
    steps_done = 0

# --- Production ---
print(f"→ Starting production: {steps_done} → {total_steps} steps")
# Reporters
sim.reporters.append(DCDReporter(f'prod_{int(total_ns*1000)}ps.dcd', steps_per_snap))
sim.reporters.append(StateDataReporter(f'prod_{int(total_ns*1000)}ps.log',
    steps_per_snap, step=True, time=True, potentialEnergy=True,
    temperature=True, volume=True, speed=True))
sim.reporters.append(CheckpointReporter(chk_file, checkpoint_steps))

# Main loop
while steps_done < total_steps:
    to_run = min(checkpoint_steps, total_steps - steps_done)
    sim.step(to_run)
    steps_done += to_run
    # CheckpointReporter already wrote the file for us
    print(f"  ✔︎ Checkpoint at {(steps_done*dt)/unit.nanoseconds:.3f} ns")

print("✔︎ Production complete")
print("Final files:", os.listdir('.'))
