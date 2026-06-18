import json
import os
import random
import shutil

from openmm import CustomExternalForce, LangevinMiddleIntegrator, Platform, XmlSerializer, unit
from openmm.app import Simulation


def pick_platform():
    requested = os.environ.get("COLABMDA_PLATFORM", "auto").strip().upper()
    names = [requested] if requested != "AUTO" else ["CUDA", "OPENCL", "CPU"]
    canonical = {"CUDA": "CUDA", "OPENCL": "OpenCL", "CPU": "CPU"}
    for name in names:
        try:
            platform = Platform.getPlatformByName(canonical.get(name, name))
        except Exception:
            if requested != "AUTO":
                raise
            pass
        else:
            if platform.getName() == "CPU" and os.environ.get("COLABMDA_REQUIRE_GPU") == "1":
                raise RuntimeError(
                    "COLABMDA_REQUIRE_GPU=1 but CUDA/OpenCL is unavailable; refusing CPU fallback"
                )
            return platform
    raise RuntimeError("No OpenMM platform available")


def fresh_seed(seed=None):
    """Return a positive OpenMM seed; zero has special non-reproducible semantics."""
    value = int(seed) if seed is not None else random.SystemRandom().randint(1, 2**31 - 1)
    if not 1 <= value <= 2**31 - 1:
        raise ValueError("OpenMM random seeds must be between 1 and 2147483647")
    return value


def offset_seed(seed, offset=1):
    """Derive a distinct valid OpenMM seed without overflowing its signed range."""
    return ((fresh_seed(seed) - 1 + int(offset)) % (2**31 - 1)) + 1


def configure_system(system, barostat_seed=None):
    """Apply the VarMDyn-compatible nonbonded and stochastic settings."""
    for force in system.getForces():
        name = force.__class__.__name__
        if name == "NonbondedForce":
            force.setUseSwitchingFunction(True)
            force.setSwitchingDistance(0.8 * unit.nanometer)
        elif name == "CMMotionRemover":
            force.setFrequency(100)
        elif name == "MonteCarloBarostat":
            force.setRandomNumberSeed(fresh_seed(barostat_seed))
    return system


def make_sim(top, sys_, dt, temp=300, seed=None, barostat_seed=None):
    configure_system(sys_, barostat_seed=barostat_seed)
    integ = LangevinMiddleIntegrator(temp * unit.kelvin, 1 / unit.picosecond, dt)
    integ.setRandomNumberSeed(fresh_seed(seed))
    plat = pick_platform()
    properties = {}
    if plat.getName() in {"CUDA", "OpenCL"}:
        precision = os.environ.get("COLABMDA_PRECISION", "mixed")
        if "Precision" in plat.getPropertyNames():
            properties["Precision"] = precision
    detail = f", precision={properties['Precision']}" if "Precision" in properties else ""
    print(f"[OpenMM] Platform: {plat.getName()}{detail}")
    return Simulation(top, sys_, integ, plat, properties)


def save_stage_state(simulation, path):
    """Save portable coordinates/velocities/box for transfer between ensembles."""
    state = simulation.context.getState(
        getPositions=True,
        getVelocities=True,
        enforcePeriodicBox=True,
    )
    temp_path = path + ".tmp"
    with open(temp_path, "w") as handle:
        handle.write(XmlSerializer.serialize(state))
    os.replace(temp_path, path)


def load_stage_state(simulation, path, load_velocities=True):
    """Load physical state without importing Context/RNG internals."""
    with open(path) as handle:
        state = XmlSerializer.deserialize(handle.read())
    simulation.context.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())
    simulation.context.setPositions(state.getPositions())
    if load_velocities:
        simulation.context.setVelocities(state.getVelocities())


def add_solute_heavy_atom_restraints(system, topology, positions, parameter_name):
    """Add VarMDyn-style Cartesian restraints to non-solvent heavy atoms."""
    force = CustomExternalForce(f"{parameter_name}*periodicdistance(x,y,z,x0,y0,z0)^2")
    force.addGlobalParameter(
        parameter_name,
        0.0 * unit.kilocalories_per_mole / unit.angstroms**2,
    )
    for coordinate in ("x0", "y0", "z0"):
        force.addPerParticleParameter(coordinate)
    solvent_or_bulk_ions = {"HOH", "WAT", "OPC", "NA", "CL", "K"}
    count = 0
    for atom, position in zip(topology.atoms(), positions):
        if atom.residue.name.upper() in solvent_or_bulk_ions:
            continue
        if atom.element is not None and atom.element.symbol == "H":
            continue
        force.addParticle(atom.index, position.value_in_unit(unit.nanometer))
        count += 1
    system.addForce(force)
    return parameter_name, count


def set_restraint_strength(context, parameter_name, kcal_per_mol_a2):
    value = (kcal_per_mol_a2 * unit.kilocalories_per_mole / unit.angstroms**2).value_in_unit(
        unit.kilojoules_per_mole / unit.nanometers**2
    )
    context.setParameter(parameter_name, value)


def save_stage_artifacts(simulation, workdir, stage, description):
    """Write a compact, auditable artifact pair for one equilibration stage."""
    stage_dir = os.path.join(workdir, "equilibration")
    os.makedirs(stage_dir, exist_ok=True)
    stem = os.path.join(stage_dir, f"stage_{int(stage):02d}")
    save_stage_state(simulation, stem + ".state.xml")
    state = simulation.context.getState(getEnergy=True)
    record = {
        "stage": int(stage),
        "description": description,
        "step": int(state.getStepCount()),
        "time_ps": state.getTime().value_in_unit(unit.picoseconds),
        "potential_energy_kj_mol": state.getPotentialEnergy().value_in_unit(
            unit.kilojoules_per_mole
        ),
    }
    summary_path = stem + ".json"
    with open(summary_path + ".tmp", "w") as handle:
        json.dump(record, handle, indent=2)
        handle.write("\n")
    os.replace(summary_path + ".tmp", summary_path)


def atomic_rename(tmp_path, final_path):
    if os.path.exists(tmp_path):
        os.replace(tmp_path, final_path)


def safe_remove_if_empty(path):
    if os.path.isfile(path) and os.path.getsize(path) == 0:
        try:
            os.remove(path)
        except Exception:
            pass


def sync_outputs(workdir, sync_dir, extra_files=None):
    if not sync_dir:
        return
    os.makedirs(sync_dir, exist_ok=True)
    essentials = [
        "system.xml",
        "solvated.pdb",
        "em.chk",
        "em.state.xml",
        "nvt.chk",
        "nvt.state.xml",
        "nvt_protocol.json",
        "npt.chk",
        "npt.state.xml",
        "npt_protocol.json",
        "parameterization.json",
        "prod.chk",
        "nvt.log",
        "npt.log",
        "prod_full.log",
    ]
    if extra_files:
        essentials += extra_files
    for fn in essentials:
        src = os.path.join(workdir, fn)
        if os.path.exists(src):
            try:
                shutil.copy2(src, os.path.join(sync_dir, fn))
            except Exception:
                pass
    stage_dir = os.path.join(workdir, "equilibration")
    if os.path.isdir(stage_dir):
        shutil.copytree(stage_dir, os.path.join(sync_dir, "equilibration"), dirs_exist_ok=True)
