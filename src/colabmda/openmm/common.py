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


def derive_replica_seed(seed=None, replica_name=None):
    """Derive a unique random seed for a replica (e.g. r1 -> offset 1, r2 -> offset 2)."""
    base = fresh_seed(seed)
    if not replica_name:
        return base
    import re

    m = re.search(r"(\d+)", str(replica_name))
    offset = int(m.group(1)) if m else sum(ord(c) for c in str(replica_name))
    return offset_seed(base, offset=offset)


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
    simulation.context.setTime(state.getTime())
    simulation.context.setStepCount(state.getStepCount())
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
        residue_name = atom.residue.name.upper().rstrip("+-0123456789")
        if residue_name in solvent_or_bulk_ions:
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
    """Write an OpenMM physical restart state for one equilibration stage."""
    stage_dir = os.path.join(workdir, "stages")
    os.makedirs(stage_dir, exist_ok=True)
    state_path = os.path.join(stage_dir, f"state_{int(stage):02d}.xml")
    save_stage_state(simulation, state_path)


def atomic_rename(tmp_path, final_path):
    if os.path.exists(tmp_path):
        os.replace(tmp_path, final_path)


def safe_remove_if_empty(path):
    if os.path.isfile(path) and os.path.getsize(path) == 0:
        try:
            os.remove(path)
        except Exception:
            pass


def find_highest_completed_stage(
    workdir: str, start_stage: int = 1, end_stage: int = 24
) -> tuple[int, str | None]:
    """Scan workdir/stages or workdir/equilibration to find the highest completed stage.

    Returns:
        (highest_stage_number, path_to_highest_state_xml)
        e.g., (23, ".../stages/state_23.xml")
        If no stage completed, returns (0, None).
    """
    stages_dir = os.path.join(workdir, "stages")
    legacy_dir = os.path.join(workdir, "equilibration")
    eq_dir = stages_dir if os.path.exists(stages_dir) else legacy_dir
    if not os.path.exists(eq_dir):
        return 0, None

    highest, latest = 0, None
    gap = False
    for stage_num in range(start_stage, end_stage + 1):
        candidates = [
            os.path.join(eq_dir, f"state_{stage_num:02d}.xml"),
            os.path.join(eq_dir, f"stage_{stage_num:02d}.state.xml"),
        ]
        state_file = next(
            (p for p in candidates if os.path.isfile(p) and os.path.getsize(p) > 0), None
        )
        json_file = os.path.join(eq_dir, f"stage_{stage_num:02d}.json")
        complete = state_file is not None
        if not complete:
            gap = True
            continue
        if gap:
            raise RuntimeError(
                "Non-contiguous equilibration stages; inspect artifacts before recovery"
            )
        if os.path.isfile(json_file):
            try:
                with open(json_file) as handle:
                    if json.load(handle).get("stage") != stage_num:
                        raise RuntimeError(f"Incorrect equilibration stage metadata: {json_file}")
            except (json.JSONDecodeError, OSError):
                pass
        highest, latest = stage_num, state_file
    return highest, latest


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
        "eq.pdb",
        "equilibrated.pdb",
        "state.xml",
        "md.dcd",
        "md.log",
        "qc.png",
        "qc.txt",
        "qc.json",
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
    for folder in ["stages", "equilibration", "chunks"]:
        src_folder = os.path.join(workdir, folder)
        if os.path.isdir(src_folder):
            shutil.copytree(src_folder, os.path.join(sync_dir, folder), dirs_exist_ok=True)
