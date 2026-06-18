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


def run_npt(workdir, pdbid, equil_time_ps=100.0, seed=None, protocol="varmdyn"):
    workdir = os.path.abspath(workdir)
    os.chdir(workdir)

    xml_system = os.path.join(workdir, "system.xml")
    pdb_saved = os.path.join(workdir, "solvated.pdb")
    state_in = os.path.join(workdir, "nvt.state.xml")
    chk_out = os.path.join(workdir, "npt.chk")
    state_out = os.path.join(workdir, "npt.state.xml")
    protocol_file = os.path.join(workdir, "npt_protocol.json")
    log_file = os.path.join(workdir, "npt.log")

    print(f"[DEBUG] Checking for NPT files in: {workdir}")
    if os.path.exists(state_out):
        if os.path.exists(protocol_file):
            with open(protocol_file) as handle:
                recorded = json.load(handle)
            if recorded.get("protocol") != protocol:
                raise RuntimeError(
                    f"Existing NPT state used protocol={recorded.get('protocol')!r}; "
                    f"requested {protocol!r}. Use a fresh workdir."
                )
        print(f"✔ NPT already complete ({state_out} found). Skipping...")
        return True

    if not all(os.path.exists(f) for f in [xml_system, pdb_saved, state_in]):
        print("Error: Missing nvt.state.xml or system files. Run NVT first.")
        return False

    print(f"▶ Starting NPT Equilibration ({equil_time_ps} ps) …")
    with open(xml_system) as f:
        system = XmlSerializer.deserializeSystem(f.read())
    pdb = PDBFile(pdb_saved)

    restraint_parameter = None
    if protocol == "varmdyn":
        restraint_parameter, restrained_atoms = add_solute_heavy_atom_restraints(
            system, pdb.topology, pdb.positions, "npt_restraint_k"
        )
        print(f"  • Restrained solute heavy atoms: {restrained_atoms}")

    dt = 2.0 * unit.femtoseconds
    seed = fresh_seed(seed)
    sim = make_sim(pdb.topology, system, dt, seed=seed, barostat_seed=offset_seed(seed))
    load_stage_state(sim, state_in)

    eq_steps = (
        50000 if protocol == "varmdyn" else max(1, int((equil_time_ps * unit.picoseconds) / dt))
    )
    sim.reporters.append(
        StateDataReporter(
            log_file,
            max(1, eq_steps // 50),
            step=True,
            time=True,
            temperature=True,
            potentialEnergy=True,
            volume=True,
            density=True,
        )
    )

    if protocol == "varmdyn":
        for label, weight in [("19", 10), ("20", 5), ("21", 2), ("22", 1)]:
            print(f"  • Stage {label}: NPT, 100 ps at 300 K/1 bar, k={weight} kcal/mol/Å²")
            set_restraint_strength(sim.context, restraint_parameter, weight)
            sim.step(eq_steps)
            save_stage_artifacts(
                sim,
                workdir,
                label,
                f"NPT; 100 ps; 300 K; 1 bar; 2 fs; k={weight} kcal/mol/A^2",
            )
        print("  • Stage 23: unrestrained NPT, 100 ps at 300 K/1 bar")
        set_restraint_strength(sim.context, restraint_parameter, 0)
        sim.step(eq_steps)
        save_stage_artifacts(sim, workdir, 23, "unrestrained NPT; 100 ps; 300 K; 1 bar; 2 fs")
        print("  • Stage 24: unrestrained NPT, 1000 ps at 300 K/1 bar")
        sim.step(500000)
        save_stage_artifacts(sim, workdir, 24, "unrestrained NPT; 1000 ps; 300 K; 1 bar; 2 fs")
    else:
        sim.step(eq_steps)

    # Save final equilibrated PDB
    with open("equilibrated.pdb", "w") as f:
        PDBFile.writeFile(sim.topology, sim.context.getState(getPositions=True).getPositions(), f)

    sim.saveCheckpoint(chk_out)
    save_stage_state(sim, state_out)
    with open(protocol_file, "w") as handle:
        json.dump({"protocol": protocol}, handle, indent=2)
        handle.write("\n")
    print("✔ NPT complete (npt.state.xml saved; npt.chk is same-stage recovery only)")
    return True


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("workdir")
    p.add_argument("--pdbid", required=True)
    p.add_argument("--equil-time", type=float, default=100.0)
    p.add_argument("--seed", type=int, default=None)
    p.add_argument("--protocol", choices=["varmdyn", "quick"], default="varmdyn")
    args = p.parse_args()
    run_npt(args.workdir, args.pdbid, args.equil_time, args.seed, args.protocol)
