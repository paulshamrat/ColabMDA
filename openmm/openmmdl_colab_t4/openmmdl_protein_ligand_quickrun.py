from openmmdl.openmmdl_simulation.scripts.forcefield_water import (
    ff_selection,
    water_forcefield_selection,
    water_model_selection,
    generate_forcefield,
)
from openmmdl.openmmdl_simulation.scripts.protein_ligand_prep import (
    prepare_ligand,
    rdkit_to_openmm,
    merge_protein_and_ligand,
    water_padding_solvent_builder,
)

import pdbfixer
import simtk.openmm.app as app
from simtk.openmm.app import PDBFile, PDBReporter, StateDataReporter, DCDReporter
from simtk.openmm import unit, Platform, LangevinMiddleIntegrator, MonteCarloBarostat
import sys

# Input files copied into the simulation working folder by `openmmdl simulation`
protein = "5wyz-moe-processed_openMMDL.pdb"
ligand = "5VF.sdf"
ligand_name = "UNK"

# Forcefield and solvent settings
ff = "AMBER14"
water = "TIP3P"
water_padding_distance = 1.0
water_ionicstrength = 0.15
small_molecule_ff = "gaff"
small_molecule_ff_version = "gaff-2.11"

# Short smoke-test simulation settings for Colab
steps = 2000
equilibration_steps = 200
dt = 0.002 * unit.picoseconds
temperature = 300 * unit.kelvin
friction = 1.0 / unit.picosecond
pressure = 1.0 * unit.atmospheres

print("Preparing ligand and protein...")
ligand_prepared = prepare_ligand(ligand, minimize_molecule=True)
omm_ligand = rdkit_to_openmm(ligand_prepared, ligand_name)
protein_pdb = pdbfixer.PDBFixer(str(protein))

forcefield_selected = ff_selection(ff)
water_selected = water_forcefield_selection(water=water, forcefield_selection=forcefield_selected)
model_water = water_model_selection(water=water, forcefield_selection=forcefield_selected)

forcefield = generate_forcefield(
    protein_ff=forcefield_selected,
    solvent_ff=water_selected,
    add_membrane=False,
    smallMoleculeForceField=small_molecule_ff,
    smallMoleculeForceFieldVersion=small_molecule_ff_version,
    rdkit_mol=ligand_prepared,
)

complex_topology, complex_positions = merge_protein_and_ligand(protein_pdb, omm_ligand)
modeller = app.Modeller(complex_topology, complex_positions)
water_padding_solvent_builder(
    model_water,
    forcefield,
    water_padding_distance,
    protein_pdb,
    modeller,
    "Na+",
    "Cl-",
    water_ionicstrength,
    protein,
)

topology = modeller.topology
positions = modeller.positions

print("Building system...")
system = forcefield.createSystem(
    topology,
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0 * unit.nanometers,
    constraints=app.HBonds,
    rigidWater=True,
    ewaldErrorTolerance=0.0005,
)
system.addForce(MonteCarloBarostat(pressure, temperature, 25))

integrator = LangevinMiddleIntegrator(temperature, friction, dt)
integrator.setConstraintTolerance(1e-6)

available = [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())]
platform_name = "CUDA" if "CUDA" in available else "CPU"
print("Available OpenMM platforms:", available)
print("Selected platform:", platform_name)

if platform_name == "CUDA":
    platform = Platform.getPlatformByName("CUDA")
    simulation = app.Simulation(topology, system, integrator, platform, {"Precision": "single"})
else:
    platform = Platform.getPlatformByName("CPU")
    simulation = app.Simulation(topology, system, integrator, platform)

simulation.context.setPositions(positions)

print("Minimizing and equilibrating...")
simulation.minimizeEnergy()
simulation.context.setVelocitiesToTemperature(temperature)
simulation.step(equilibration_steps)

print("Running short production...")
simulation.reporters.append(PDBReporter("output_5wyz_quickrun.pdb", 500))
simulation.reporters.append(DCDReporter("trajectory_5wyz_quickrun.dcd", 500))
simulation.reporters.append(
    StateDataReporter(
        "log_5wyz_quickrun.txt",
        100,
        totalSteps=steps,
        step=True,
        potentialEnergy=True,
        temperature=True,
        speed=True,
        progress=True,
        separator="\t",
    )
)
simulation.reporters.append(
    StateDataReporter(
        sys.stdout,
        100,
        totalSteps=steps,
        step=True,
        potentialEnergy=True,
        temperature=True,
        separator="\t",
    )
)
simulation.step(steps)

with open("final_state_5wyz_quickrun.pdb", "w") as handle:
    PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), handle)

print("Quick protein-ligand simulation completed.")
