import argparse
import json
import os

from openmm import XmlSerializer, unit
from openmm.app import PDBFile, StateDataReporter

from colabmda.openmm.common import (
    add_solute_heavy_atom_restraints,
    find_highest_completed_stage,
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
        if protocol == "varmdyn":
            completed, _ = find_highest_completed_stage(workdir, 19, 24)
            if completed != 24:
                raise RuntimeError(
                    "Final NPT state lacks complete staged artifacts; inspect before recovery"
                )
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

    # Check for granular stage recovery
    highest_stage, latest_state_path = find_highest_completed_stage(
        workdir, start_stage=19, end_stage=24
    )
    # Ensure required input files exist
    if highest_stage >= 19 and latest_state_path:
        start_stage = highest_stage + 1
        starting_state = latest_state_path
        print(f"[RECOVERY] Found existing progress up to Stage {highest_stage}!")
        print(
            f"[RECOVERY] Resuming directly from {os.path.basename(latest_state_path)} (starting Stage {start_stage})..."
        )
        if not all(os.path.exists(f) for f in [xml_system, pdb_saved]):
            print("Error: Missing system files. Cannot resume NPT.")
            return False
    else:
        start_stage = 19
        starting_state = state_in
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
    load_stage_state(sim, starting_state)

    test_mode = os.environ.get("COLABMDA_TEST_1_24", "").strip() in {"1", "true", "True"}
    eq_steps = (
        (10 if test_mode else 50000)
        if protocol == "varmdyn"
        else max(1, int((equil_time_ps * unit.picoseconds) / dt))
    )
    report_interval = 2 if test_mode else max(1, eq_steps // 50)
    sim.reporters.append(
        StateDataReporter(
            log_file,
            report_interval,
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
            stage_num = int(label)
            if stage_num < start_stage:
                print(f"  • Stage {label}: already complete. Fast-forwarding...")
                continue
            print(f"  • Stage {label}: NPT, 100 ps at 300 K/1 bar, k={weight} kcal/mol/Å²")
            set_restraint_strength(sim.context, restraint_parameter, weight)
            sim.step(eq_steps)
            save_stage_artifacts(
                sim,
                workdir,
                label,
                f"NPT; 100 ps; 300 K; 1 bar; 2 fs; k={weight} kcal/mol/A^2",
            )

        if 23 >= start_stage:
            print("  • Stage 23: unrestrained NPT, 100 ps at 300 K/1 bar")
            set_restraint_strength(sim.context, restraint_parameter, 0)
            sim.step(eq_steps)
            save_stage_artifacts(sim, workdir, 23, "unrestrained NPT; 100 ps; 300 K; 1 bar; 2 fs")
        else:
            print("  • Stage 23: already complete. Fast-forwarding...")

        if 24 >= start_stage:
            print("  • Stage 24: unrestrained NPT, 1000 ps at 300 K/1 bar")
            stage24_steps = 10 if test_mode else 500000
            set_restraint_strength(sim.context, restraint_parameter, 0)
            sim.step(stage24_steps)
            save_stage_artifacts(sim, workdir, 24, "unrestrained NPT; 1000 ps; 300 K; 1 bar; 2 fs")
        else:
            print("  • Stage 24: already complete. Fast-forwarding...")
    else:
        sim.step(eq_steps)

    # Save final equilibrated PDB (standard OpenMM convention: eq.pdb)
    final_positions = sim.context.getState(getPositions=True).getPositions()
    with open(os.path.join(workdir, "eq.pdb"), "w") as f:
        PDBFile.writeFile(sim.topology, final_positions, f)

    sim.saveCheckpoint(chk_out)
    save_stage_state(sim, state_out)
    print("✔ NPT complete (npt.state.xml and eq.pdb saved)")

    # Generate automated QC report if possible
    try:
        from colabmda.openmm.qc import analyze_logs

        analyze_logs(workdir, strict=False, protocol=protocol)
    except Exception:
        pass
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
