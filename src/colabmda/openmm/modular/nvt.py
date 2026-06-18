import argparse
import json
import os

from openmm import XmlSerializer, unit
from openmm.app import PDBFile, StateDataReporter

from colabmda.openmm.modular.utils import (
    add_solute_heavy_atom_restraints,
    fresh_seed,
    load_stage_state,
    make_sim,
    offset_seed,
    save_stage_artifacts,
    save_stage_state,
    set_restraint_strength,
)


def run_nvt(workdir, pdbid, equil_time_ps=100.0, seed=None, protocol="varmdyn"):
    workdir = os.path.abspath(workdir)
    os.chdir(workdir)

    xml_system = os.path.join(workdir, "system.xml")
    pdb_saved = os.path.join(workdir, "solvated.pdb")
    state_in = os.path.join(workdir, "em.state.xml")
    chk_out = os.path.join(workdir, "nvt.chk")
    state_out = os.path.join(workdir, "nvt.state.xml")
    protocol_file = os.path.join(workdir, "nvt_protocol.json")
    log_file = os.path.join(workdir, "nvt.log")

    print(f"[DEBUG] Checking for NVT files in: {workdir}")
    print(f"[DEBUG] xml: {xml_system}, pdb: {pdb_saved}, in: {state_in}, out: {state_out}")
    if os.path.exists(state_out):
        if os.path.exists(protocol_file):
            with open(protocol_file) as handle:
                recorded = json.load(handle)
            if recorded.get("protocol") != protocol:
                raise RuntimeError(
                    f"Existing NVT state used protocol={recorded.get('protocol')!r}; "
                    f"requested {protocol!r}. Use a fresh workdir."
                )
        print(f"✔ NVT already complete ({state_out} found). Skipping...")
        return True

    if not all(os.path.exists(f) for f in [xml_system, pdb_saved, state_in]):
        print("Error: Missing em.state.xml or system files. Run EM with the current version first.")
        return False

    print(f"▶ Starting NVT Equilibration ({equil_time_ps} ps) …")
    with open(xml_system) as f:
        system = XmlSerializer.deserializeSystem(f.read())
    pdb = PDBFile(pdb_saved)

    # Remove MonteCarloBarostat for true NVT constant-volume equilibration
    barostat_indices = [
        i
        for i in range(system.getNumForces())
        if system.getForce(i).__class__.__name__ == "MonteCarloBarostat"
    ]
    for idx in reversed(barostat_indices):
        system.removeForce(idx)

    restraint_parameter = None
    if protocol == "varmdyn":
        restraint_parameter, restrained_atoms = add_solute_heavy_atom_restraints(
            system, pdb.topology, pdb.positions, "nvt_restraint_k"
        )
        print(f"  • Restrained solute heavy atoms: {restrained_atoms}")

    seed = fresh_seed(seed)
    dt = (1.0 if protocol == "varmdyn" else 2.0) * unit.femtoseconds
    sim = make_sim(pdb.topology, system, dt, seed=seed, temp=50 if protocol == "varmdyn" else 300)
    load_stage_state(sim, state_in, load_velocities=False)

    # Assign velocities
    print(f"  • Setting velocities (seed={seed})")
    sim.context.setVelocitiesToTemperature(
        (50 if protocol == "varmdyn" else 300) * unit.kelvin, seed
    )

    if protocol == "varmdyn":
        for label, weight in [("01", 100), ("02", 50), ("03", 25)]:
            print(f"  • Stage {label}: restrained minimization, k={weight} kcal/mol/Å²")
            set_restraint_strength(sim.context, restraint_parameter, weight)
            sim.minimizeEnergy(maxIterations=5000)
            save_stage_artifacts(
                sim, workdir, label, f"restrained minimization; k={weight} kcal/mol/A^2"
            )
        for index, temperature in enumerate((50, 100, 150, 200, 250, 300), start=4):
            print(f"  • Stage {index:02d}: NVT heating, 10 ps at {temperature} K, 1 fs")
            sim.integrator.setTemperature(temperature * unit.kelvin)
            sim.step(10000)
            save_stage_artifacts(
                sim, workdir, index, f"NVT heating; 10 ps; {temperature} K; 1 fs; k=25"
            )
        for label, weight, iterations in [
            ("10", 25, 5000),
            ("11", 10, 1000),
            ("12", 5, 1000),
            ("13", 2, 1000),
            ("14", 1, 1000),
        ]:
            print(f"  • Stage {label}: restrained minimization, k={weight} kcal/mol/Å²")
            set_restraint_strength(sim.context, restraint_parameter, weight)
            sim.minimizeEnergy(maxIterations=iterations)
            save_stage_artifacts(
                sim, workdir, label, f"post-heating minimization; k={weight} kcal/mol/A^2"
            )
        # Rebuild at 2 fs for constrained 300 K equilibration, transferring only physical state.
        save_stage_state(sim, "nvt_heating.state.xml")
        system_2fs = XmlSerializer.deserializeSystem(XmlSerializer.serializeSystem(system))
        sim = make_sim(pdb.topology, system_2fs, 2.0 * unit.femtoseconds, seed=offset_seed(seed))
        load_stage_state(sim, "nvt_heating.state.xml")
        restraint_parameter = "nvt_restraint_k"
        eq_steps = 50000
    else:
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

    if protocol == "varmdyn":
        for label, weight in [("15", 10), ("16", 5), ("17", 2), ("18", 1)]:
            print(f"  • Stage {label}: NVT, 100 ps at 300 K, k={weight} kcal/mol/Å²")
            set_restraint_strength(sim.context, restraint_parameter, weight)
            sim.step(eq_steps)
            save_stage_artifacts(
                sim, workdir, label, f"NVT; 100 ps; 300 K; 2 fs; k={weight} kcal/mol/A^2"
            )
    else:
        sim.step(eq_steps)
    sim.saveCheckpoint(chk_out)
    save_stage_state(sim, state_out)
    with open(protocol_file, "w") as handle:
        json.dump({"protocol": protocol}, handle, indent=2)
        handle.write("\n")
    print("✔ NVT complete (nvt.state.xml saved; nvt.chk is same-stage recovery only)")
    return True


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("workdir")
    p.add_argument("--pdbid", required=True)
    p.add_argument("--equil-time", type=float, default=100.0)
    p.add_argument("--seed", type=int, default=None)
    p.add_argument("--protocol", choices=["varmdyn", "quick"], default="varmdyn")
    args = p.parse_args()
    run_nvt(args.workdir, args.pdbid, args.equil_time, args.seed, args.protocol)
