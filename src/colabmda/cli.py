#!/usr/bin/env python3
import argparse
import os
import shutil
from pathlib import Path

DEFAULT_DRIVE_ROOT = "/content/drive/MyDrive/research/projects"
ENV_ROOT = "COLABMDA_ROOT"
ENV_PROJECT = "COLABMDA_PROJECT_DIR"
ROOT_HELP = "Override project root (defaults to current working directory)"


def _normalize_equil_protocol(protocol: str) -> str:
    """Map user-facing equilibration labels to internal implementation names."""
    if protocol in {"default", "full", "staged"}:
        return "varmdyn"
    if protocol in {"varmdyn", "quick"}:
        return protocol
    raise SystemExit(
        "ERROR: --equil-protocol must be 'default' for the full staged protocol "
        "or 'quick' for a short smoke-test equilibration."
    )


def _resolve_root(use_drive: bool, root: str | None) -> str:
    # 1. Use explicit root if provided
    if root:
        return str(Path(root).resolve())

    # 2. Use environment variable if provided
    for env_var in [ENV_PROJECT, ENV_ROOT]:
        val = os.environ.get(env_var)
        if val:
            return str(Path(val).resolve())

    cwd = os.getcwd()
    # 3. If on Colab but in the temporary /content folder, check Drive
    if cwd == "/content" and os.path.exists("/content/drive/MyDrive"):
        if os.path.exists("/content/drive/MyDrive/research/projects"):
            return "/content/drive/MyDrive/research/projects"
        if os.path.exists("/content/drive/MyDrive/ColabMDA"):
            return "/content/drive/MyDrive/ColabMDA"
        return "/content/drive/MyDrive"

    # 4. Otherwise, use the Current Working Directory (standard for HPC/Laptop/Colab projects)
    return cwd


def _ensure_dir(path: str):
    Path(path).mkdir(parents=True, exist_ok=True)


def _prepare_run_inputs(root: str, pdbid: str, name: str):
    prep_dir = Path(root) / pdbid / "prep"
    run_dir = Path(root) / pdbid / "run"
    run_dir.mkdir(parents=True, exist_ok=True)

    prep_clean = prep_dir / f"{name}_cleaned.pdb"
    run_clean = run_dir / f"{name}_cleaned.pdb"
    if not run_clean.exists() and prep_clean.exists():
        shutil.copy2(prep_clean, run_clean)

    prep_raw = prep_dir / f"{name}.pdb"
    run_raw = run_dir / f"{name}.pdb"
    if not run_raw.exists() and prep_raw.exists():
        shutil.copy2(prep_raw, run_raw)

    if not run_clean.exists():
        raise SystemExit(
            f"ERROR: cleaned PDB missing for run.\n"
            f"  Expected at: {run_clean}\n"
            f"Run prep first, e.g.:\n"
            f"  colabmda openmm prep --pdb-id {pdbid}\n"
        )


def _guess_pdbid_from_workdir(workdir: str) -> str | None:
    candidates = list(Path(workdir).glob("*_cleaned.pdb"))
    if not candidates:
        return None
    if len(candidates) > 1:
        names = ", ".join([c.name for c in candidates])
        raise SystemExit(
            f"ERROR: multiple *_cleaned.pdb files found in {workdir}: {names}\n"
            f"Please specify --name."
        )
    return candidates[0].name.replace("_cleaned.pdb", "")


def _default_project_root() -> str:
    return _resolve_root(use_drive=True, root=None) or DEFAULT_DRIVE_ROOT


def _simulation_workdir(root: str, name: str) -> Path:
    """Resolve a named simulation under the current and legacy layouts."""
    root_path = Path(root).resolve()
    candidates = [
        root_path / "data" / "simulation" / name,
        root_path / "data" / "sim" / name,
        root_path / "simulation" / name,
        root_path / "sim" / name,
        root_path / "simulations" / name,
    ]
    return next((path for path in candidates if path.is_dir()), candidates[0])


def _analysis_outdir(root: str, name: str) -> Path:
    """Resolve a named analysis directory under the current and legacy layouts."""
    root_path = Path(root).resolve()
    candidates = [
        root_path / "data" / "analysis" / name,
        root_path / "analysis" / name,
    ]
    return next((path for path in candidates if path.is_dir()), candidates[0])


def main():
    p = argparse.ArgumentParser(prog="colabmda")
    sub = p.add_subparsers(dest="tool", required=True)

    # ---------------- OpenMM ----------------
    p_openmm = sub.add_parser("openmm", help="OpenMM protein-water workflow")
    sub_openmm = p_openmm.add_subparsers(dest="cmd", required=True)

    # prep
    p_prep = sub_openmm.add_parser("prep", help="Prepare/clean PDB")
    g = p_prep.add_mutually_exclusive_group(required=True)
    g.add_argument("--pdb-id", help="4-letter PDB id (downloads from RCSB)")
    g.add_argument("--pdb-file", help="Local PDB file path (no download)")
    p_prep.add_argument(
        "--outdir", default=None, help="Output directory for --pdb-file (default: ./<name>)"
    )
    p_prep.add_argument(
        "--name", default=None, help="Prefix name (default: from --pdb-id or file stem)"
    )
    p_prep.add_argument("--ph", type=float, default=7.4, help="Hydrogen pH (default: 7.4)")
    p_prep.add_argument(
        "--sync-dir", default=None, help="Optional: copy cleaned prep outputs to this directory"
    )
    p_prep.add_argument(
        "--drive", action="store_true", help="(compat) Use Drive root (default behavior)"
    )
    p_prep.add_argument(
        "--root",
        "--project",
        dest="root",
        default=None,
        help=ROOT_HELP,
    )

    # run (colab-safe runner)
    p_run = sub_openmm.add_parser("run", help="Run/resume chunked MD (colab-safe)")
    g = p_run.add_mutually_exclusive_group(required=False)
    g.add_argument("--pdb-id", help="Use ./<pdb-id> as workdir (pdb-id download workflow)")
    g.add_argument("--workdir", help="Folder containing <name>_cleaned.pdb")
    p_run.add_argument(
        "--name", default=None, help="Prefix name (default: --pdb-id or inferred from workdir)"
    )
    p_run.add_argument("--total-ns", type=float, default=100.0)
    p_run.add_argument(
        "--traj-interval",
        type=float,
        default=10.0,
        help="ps between saved frames (default: 10.0 ps)",
    )
    p_run.add_argument(
        "--equil-time",
        type=float,
        default=100.0,
        help="NVT and NPT duration in ps for --equil-protocol quick",
    )
    p_run.add_argument("--checkpoint-ps", type=float, default=1000.0, help="ps per chunk")
    p_run.add_argument("--sync-dir", default=None, help="Optional: sync outputs to this directory")
    p_run.add_argument("--replica", default=None, help="Optional: replica subfolder (e.g. r1, r2, r3)")
    p_run.add_argument("--ligand", default=None, help="Optional: Path to ligand SDF/MOL2 file")
    p_run.add_argument(
        "--small-molecule-ff",
        default="gaff-2.11",
        help="Small-molecule force field for ligand setup (default: gaff-2.11; e.g. gaff-2.2.20 if installed)",
    )
    p_run.add_argument(
        "--keep-mg", action="store_true", help="Optional: Keep Mg2+ ion from raw structure"
    )
    p_run.add_argument(
        "--padding-nm",
        default="auto",
        help="Solvent padding in nm, or auto; 1.4 is conservative and 1.0 is compact",
    )
    p_run.add_argument(
        "--protein-ff",
        choices=["ff19SB", "ff14SB"],
        default="ff19SB",
        help="Protein force field for native OpenMM setup (default: ff19SB)",
    )
    p_run.add_argument(
        "--water-model",
        choices=["opc", "tip3p", "tip3pfb", "tip4pew"],
        default="opc",
        help="Water/ion model for native OpenMM setup (default: opc)",
    )
    p_run.add_argument(
        "--equil-protocol",
        default="default",
        metavar="{default,quick}",
        help="Full staged equilibration or short smoke-test equilibration",
    )
    p_run.add_argument("--amber-prmtop", default=None, help="Pre-parameterized AMBER topology")
    p_run.add_argument("--amber-inpcrd", default=None, help="Matching AMBER coordinates")
    p_run.add_argument(
        "--seed", type=int, default=None, help="Optional: random seed for velocity assignment"
    )
    p_run.add_argument(
        "--equil-only",
        action="store_true",
        help="Run Energy Minimization, NVT, and NPT equilibration, then exit before production MD",
    )
    p_run.add_argument(
        "--drive", action="store_true", help="(compat) Use Drive root (default behavior)"
    )
    p_run.add_argument(
        "--root",
        "--project",
        dest="root",
        default=None,
        help=ROOT_HELP,
    )

    # merge
    p_merge = sub_openmm.add_parser("merge", help="Merge chunk DCDs/logs")
    g = p_merge.add_mutually_exclusive_group(required=False)
    g.add_argument("--pdb-id", help="Use ./<pdb-id> as simulation directory")
    g.add_argument("--pdb-dir", help="Simulation directory (e.g. 4ldj_wt)")
    g.add_argument("--name", help="Simulation system name (e.g. 4ldj_wt)")
    p_merge.add_argument("--replica", default="r1", help="Replica subfolder (default: r1)")
    p_merge.add_argument(
        "--topology", default=None, help="Topology PDB (default: <dir>/solvated.pdb)"
    )
    p_merge.add_argument("--out-traj", default="prod_full.dcd")
    p_merge.add_argument("--out-log", default="prod_full.log")
    p_merge.add_argument(
        "--stride", type=int, default=1, help="Keep every Nth frame while merging (default: 1)"
    )
    p_merge.add_argument("--center", action="store_true", help="Center protein in the box")
    p_merge.add_argument(
        "--wrap", action="store_true", help="Wrap solvent molecules (image_molecules)"
    )
    p_merge.add_argument(
        "--mda", action="store_true", help="Use MDAnalysis for merging and PBC correction"
    )
    p_merge.add_argument(
        "--selection", default=None, help="MDAnalysis selection string (implies --mda)"
    )
    p_merge.add_argument(
        "--ca-only", action="store_true", help="Extract Cα atoms only (implies --mda)"
    )
    p_merge.add_argument(
        "--protein-only", action="store_true", help="Extract protein-only atoms (implies --mda)"
    )
    p_merge.add_argument(
        "--drive", action="store_true", help="(compat) Use Drive root (default behavior)"
    )
    p_merge.add_argument(
        "--root",
        "--project",
        dest="root",
        default=None,
        help=ROOT_HELP,
    )

    # analysis
    p_ana = sub_openmm.add_parser("analysis", help="RMSD/Rg/RMSF analysis")
    g = p_ana.add_mutually_exclusive_group(required=False)
    g.add_argument("--pdb-id", help="Use ./<pdb-id> as simulation directory")
    g.add_argument("--pdb-dir", help="Simulation directory (e.g. 4ldj_wt)")
    g.add_argument("--name", help="Simulation system name (e.g. 4ldj_wt)")
    p_ana.add_argument("--replica", default=None, help="Replica subfolder (default: all replicas if present, or r1)")
    p_ana.add_argument("--topology", default=None)
    p_ana.add_argument("--trajectory", default=None)
    p_ana.add_argument(
        "--interval", type=float, default=None, help="ps per frame (if not auto-detected)"
    )
    p_ana.add_argument("--outdir", default=None)
    p_ana.add_argument(
        "--drive", action="store_true", help="(compat) Use Drive root (default behavior)"
    )
    p_ana.add_argument(
        "--root",
        "--project",
        dest="root",
        default=None,
        help=ROOT_HELP,
    )

    # status
    p_stat = sub_openmm.add_parser("status", help="Sanity-check frames/time/resume readiness")
    g = p_stat.add_mutually_exclusive_group(required=False)
    g.add_argument("--pdb-id", help="Use ./<pdb-id> as simulation directory")
    g.add_argument("--pdb-dir", help="Simulation directory (e.g. 4ldj_wt)")
    g.add_argument("--name", help="Simulation system name (e.g. 4ldj_wt)")
    p_stat.add_argument("--replica", default="r1", help="Replica subfolder (default: r1)")
    p_stat.add_argument(
        "--drive", action="store_true", help="(compat) Use Drive root (default behavior)"
    )
    p_stat.add_argument(
        "--root",
        "--project",
        dest="root",
        default=None,
        help=ROOT_HELP,
    )

    # compare
    p_comp = sub_openmm.add_parser("compare", help="Compare multiple systems (aggregate replicas)")
    p_comp.add_argument(
        "--series",
        action="append",
        required=True,
        help="LABEL=DIR1,DIR2 (e.g. WT=data/analysis/wt/r1,data/analysis/wt/r2)",
    )
    p_comp.add_argument("--outdir", required=True, help="Output directory for plots")

    # view
    p_view = sub_openmm.add_parser("view", help="View trajectory in PyMOL")
    p_view.add_argument(
        "--pdb-dir",
        default=None,
        help="Directory containing the merged files (default: current directory)",
    )
    p_view.add_argument(
        "-t", "--topology", default=None, help="Topology file (default: prod_full.pdb)"
    )
    p_view.add_argument(
        "-x", "--trajectory", default=None, help="Trajectory file (default: prod_full.dcd)"
    )
    p_view.add_argument(
        "-r", "--resi", type=int, default=12, help="Highlight residue index (default: 12)"
    )

    # snapshots
    p_snap = sub_openmm.add_parser(
        "snapshots", help="Generate comparative snapshot grid (e.g. 3x8 WT vs Mutants)"
    )
    p_snap.add_argument(
        "-c", "--config", default=None, help="Path to custom JSON configuration file"
    )
    p_snap.add_argument(
        "-o",
        "--output",
        default="data/analysis/snapshots/master_3x8_snapshots.png",
        help="Path to output image grid file",
    )
    p_snap.add_argument(
        "--temp-dir", default="scratch/master_grid_3x8", help="Temporary folder for panel rendering"
    )
    p_snap.add_argument("--font-path", default=None, help="Custom TrueType font path")

    # stage
    p_stage = sub_openmm.add_parser(
        "stage", help="Stage a WT/mutant structure into simulations/<name>"
    )
    p_stage.add_argument(
        "--pdb-file", required=True, help="Input structure PDB file (typically from data/structure/)"
    )
    p_stage.add_argument("--name", required=True, help="Simulation name (e.g. 4ldj_wt, 4ldj_G12C)")
    p_stage.add_argument(
        "--replica", default=None, help="Optional: create nested replica subfolder (e.g. r1, r2, r3)"
    )
    p_stage.add_argument("--ph", type=float, default=7.4, help="Hydrogen pH (default: 7.4)")
    p_stage.add_argument("--root", "--project", dest="root", default=None, help=ROOT_HELP)

    # em/nvt/npt/check-equil/md (individual modular steps)
    p_em = sub_openmm.add_parser("em", help="Modular: Minimization")
    p_em.add_argument("--name", required=True)
    p_em.add_argument("--workdir", default=None)
    p_em.add_argument("--root", "--project", dest="root", default=None)
    p_em.add_argument("--ligand", default=None, help="Optional: Path to ligand SDF/MOL2 file")
    p_em.add_argument(
        "--small-molecule-ff",
        default="gaff-2.11",
        help="Small-molecule force field for ligand setup",
    )
    p_em.add_argument(
        "--keep-mg", action="store_true", help="Optional: Keep Mg2+ ion from raw structure"
    )
    p_em.add_argument("--padding-nm", default="auto")
    p_em.add_argument("--protein-ff", choices=["ff19SB", "ff14SB"], default="ff19SB")
    p_em.add_argument(
        "--water-model", choices=["opc", "tip3p", "tip3pfb", "tip4pew"], default="opc"
    )
    p_em.add_argument("--amber-prmtop", default=None)
    p_em.add_argument("--amber-inpcrd", default=None)

    p_nvt = sub_openmm.add_parser("nvt", help="Modular: NVT Equilibration")
    p_nvt.add_argument("--name", required=True)
    p_nvt.add_argument("--equil-time", type=float, default=100.0)
    p_nvt.add_argument("--seed", type=int, default=None)
    p_nvt.add_argument("--equil-protocol", default="default", metavar="{default,quick}")
    p_nvt.add_argument("--workdir", default=None)
    p_nvt.add_argument("--root", "--project", dest="root", default=None)

    p_npt = sub_openmm.add_parser("npt", help="Modular: NPT Equilibration")
    p_npt.add_argument("--name", required=True)
    p_npt.add_argument("--equil-time", type=float, default=100.0)
    p_npt.add_argument("--seed", type=int, default=None)
    p_npt.add_argument("--equil-protocol", default="default", metavar="{default,quick}")
    p_npt.add_argument("--workdir", default=None)
    p_npt.add_argument("--root", "--project", dest="root", default=None)

    p_qc = sub_openmm.add_parser("qc", aliases=["check-equil"], help="Modular: Stability Check & QC Plots")
    p_qc.add_argument("--name", required=True)
    p_qc.add_argument("--workdir", default=None)
    p_qc.add_argument("--root", "--project", dest="root", default=None)
    p_qc.add_argument(
        "--warn-only", action="store_true", help="Report failed QC without a nonzero exit"
    )

    p_md = sub_openmm.add_parser("md", help="Modular: Production MD")
    p_md.add_argument("--name", required=True)
    p_md.add_argument("--total-ns", type=float, default=100.0)
    p_md.add_argument("--traj-interval", type=float, default=10.0)
    p_md.add_argument("--checkpoint-ps", type=float, default=1000.0)
    p_md.add_argument("--sync-dir", default=None)
    p_md.add_argument("--workdir", default=None)
    p_md.add_argument("--root", "--project", dest="root", default=None)
    p_md.add_argument("--replica", default="r1", help="Production replica subfolder (default: r1)")
    p_md.add_argument("--seed", type=int, default=None, help="Random seed for velocity assignment")

    # ---------------- Modeller ----------------
    p_mod = sub.add_parser("modeller", help="Modeller workflows")
    sub_mod = p_mod.add_subparsers(dest="cmd", required=True)

    # build
    p_build = sub_mod.add_parser("build", help="Build homology model")
    p_build.add_argument("--pdb-id", required=True, help="Template PDB ID (e.g. 4bgq)")
    p_build.add_argument("--uniprot-id", required=True, help="UniProt ID (e.g. O76039)")
    p_build.add_argument("--chain", default="A", help="Chain ID (default: A)")
    p_build.add_argument("--range", nargs=2, type=int, metavar=("START", "END"))
    p_build.add_argument("--truncate", action="store_true")
    p_build.add_argument("--uniprot-numbering", action="store_true")
    p_build.add_argument("--mut", default=None)
    p_build.add_argument("--list", default=None, help="File with one mutation per line")
    p_build.add_argument("--outdir", default=None)
    p_build.add_argument("--outdir-mut", default=None)
    p_build.add_argument("--seed", type=int, default=None)
    p_build.add_argument("--logfile", default=None)
    p_build.add_argument("--verbose", action="store_true")

    # mutate
    p_mut = sub_mod.add_parser("mutate", help="Mutate an existing PDB")
    p_mut.add_argument("--pdb-in", required=True, help="Input PDB to mutate")
    p_mut.add_argument("--chain", default="A", help="Chain ID (default: A)")
    p_mut.add_argument("--mut", default=None)
    p_mut.add_argument("--list", default=None, help="File with one mutation per line")
    p_mut.add_argument("--outdir-mut", default=None)
    p_mut.add_argument("--seed", type=int, default=None)
    p_mut.add_argument("--logfile", default=None)
    p_mut.add_argument("--verbose", action="store_true")

    args = p.parse_args()

    if args.tool == "openmm":
        from colabmda.openmm.commands import (
            openmm_analysis,
            openmm_check_equil,
            openmm_compare,
            openmm_em,
            openmm_md,
            openmm_merge,
            openmm_npt,
            openmm_nvt,
            openmm_prep_from_file,
            openmm_prep_from_pdbid,
            openmm_qc,
            openmm_snapshots,
            openmm_status,
            openmm_view,
        )

        if args.cmd == "prep":
            if args.pdb_id:
                root = _resolve_root(args.drive, args.root)
                if root:
                    _ensure_dir(root)
                openmm_prep_from_pdbid(args.pdb_id, root_dir=root, ph=args.ph, sync_dir=args.sync_dir)
            else:
                name = args.name or Path(args.pdb_file).stem
                root = _resolve_root(args.drive, args.root)
                if root:
                    _ensure_dir(root)
                outdir = args.outdir or (str(Path(root) / name / "prep") if root else name)
                openmm_prep_from_file(
                    args.pdb_file, outdir, pdbid=name, ph=args.ph, sync_dir=args.sync_dir
                )

        elif args.cmd == "stage":
            name = args.name or Path(args.pdb_file).stem
            root_path = Path(
                _resolve_root(getattr(args, "drive", False), getattr(args, "root", None))
            ).resolve()
            _ensure_dir(str(root_path))

            # Priority: data/simulation > data/sim > simulation > sim > simulations (legacy)
            if (root_path / "data" / "simulation").exists():
                sim_base = root_path / "data" / "simulation"
            elif (root_path / "data" / "sim").exists():
                sim_base = root_path / "data" / "sim"
            elif (root_path / "simulation").exists():
                sim_base = root_path / "simulation"
            elif (root_path / "sim").exists():
                sim_base = root_path / "sim"
            elif (root_path / "simulations").exists():
                sim_base = root_path / "simulations"
            else:
                sim_base = root_path / "data" / "simulation"

            sim_dir = sim_base / name
            sim_dir.mkdir(parents=True, exist_ok=True)

            target_clean = sim_dir / f"{name}_cleaned.pdb"
            in_pdb = Path(args.pdb_file).resolve()

            if in_pdb.exists():
                shutil.copy2(in_pdb, target_clean)
                summary_json = in_pdb.parent / "protonation_summary.json"
                if summary_json.exists():
                    shutil.copy2(summary_json, sim_dir / "protonation_summary.json")
            else:
                openmm_prep_from_file(str(in_pdb), str(sim_dir), pdbid=name, ph=args.ph)

            # Always stage the shared equilibration hub 'eq'
            eq_dir = sim_dir / "eq"
            eq_dir.mkdir(parents=True, exist_ok=True)
            if target_clean.exists():
                shutil.copy2(target_clean, eq_dir / f"{name}_cleaned.pdb")
            if (sim_dir / "protonation_summary.json").exists():
                shutil.copy2(sim_dir / "protonation_summary.json", eq_dir / "protonation_summary.json")

            if args.replica:
                rep_dir = sim_dir / args.replica
                rep_dir.mkdir(parents=True, exist_ok=True)
                if target_clean.exists():
                    shutil.copy2(target_clean, rep_dir / f"{name}_cleaned.pdb")
                if (sim_dir / "protonation_summary.json").exists():
                    shutil.copy2(sim_dir / "protonation_summary.json", rep_dir / "protonation_summary.json")

            print(f"[INFO] Staged simulation folder: {sim_dir}")

        elif args.cmd == "run":
            if args.pdb_id:
                root = _resolve_root(args.drive, args.root)
                if root:
                    _ensure_dir(root)
                workdir = (
                    str(Path(root) / args.pdb_id / "run")
                    if root
                    else str(Path(args.pdb_id) / "run")
                )
                name = args.name or args.pdb_id
                if root:
                    _prepare_run_inputs(root=root, pdbid=args.pdb_id, name=name)
            elif args.workdir:
                workdir = str(Path(args.workdir).resolve())
                name = args.name or _guess_pdbid_from_workdir(workdir)
            else:
                root = Path(_resolve_root(args.drive, args.root)).resolve()
                name = args.name
                cwd = Path.cwd().resolve()

                print(f"[DEBUG] root: {root}")
                print(f"[DEBUG] cwd:  {cwd}")
                print(f"[DEBUG] name: {name}")

                # Search strategy: prioritize data/simulation > data/sim > simulation > sim > simulations
                search_paths = [
                    cwd / "data" / "simulation" / name if name else None,
                    cwd / "data" / "sim" / name if name else None,
                    cwd / "simulation" / name if name else None,
                    cwd / "sim" / name if name else None,
                    cwd / "simulations" / name if name else None,
                    cwd / name if name else None,
                    cwd if (name and name in str(cwd)) else None,  # Check if we are already inside
                    root / "data" / "simulation" / name if name else None,
                    root / "data" / "sim" / name if name else None,
                    root / "simulation" / name if name else None,
                    root / "sim" / name if name else None,
                    root / "simulations" / name if name else None,
                    root / name if name else None,
                ]

                workdir = None
                for p in search_paths:
                    if p:
                        exists = p.exists()
                        is_dir = p.is_dir() if exists else False
                        print(f"[DEBUG] Checking: {p} (exists={exists}, is_dir={is_dir})")
                        if exists and is_dir:
                            # If we found the replica folder itself, use its parent as the simulation base
                            if p.name == args.replica:
                                workdir = p.parent
                            else:
                                workdir = p
                            break

                if workdir is None:
                    # If we are in a folder named 'name', use it
                    if name and cwd.name == name:
                        workdir = cwd
                    else:
                        print(f"[DEBUG] No folder found, defaulting to cwd: {cwd}")
                        workdir = cwd

                if not name:
                    name = _guess_pdbid_from_workdir(str(workdir))

            workdir = Path(workdir).resolve()
            if args.replica:
                # If we are already inside the replica folder, don't double-append
                if workdir.name != args.replica:
                    workdir = workdir / args.replica

                _ensure_dir(str(workdir))

                # Ensure cleaned PDB is available in the replica folder
                local_clean = workdir / f"{name}_cleaned.pdb"
                parent_clean = workdir.parent / f"{name}_cleaned.pdb"
                if not local_clean.exists() and parent_clean.exists():
                    shutil.copy2(parent_clean, local_clean)
                elif not local_clean.exists():
                    # Fallback: find any cleaned PDB in parent to infer name
                    guesses = list(workdir.parent.glob("*_cleaned.pdb"))
                    if guesses and not name:
                        name = guesses[0].name.replace("_cleaned.pdb", "")
                        shutil.copy2(guesses[0], local_clean)

            workdir_str = str(workdir.resolve())
            if not name:
                name = _guess_pdbid_from_workdir(workdir_str)

            if not name:
                raise SystemExit("ERROR: could not infer pdbid; please specify --name.")

            # Auto-discover cleaned starting structure if missing from workdir
            workdir_path = Path(workdir_str)
            clean_target = workdir_path / f"{name}_cleaned.pdb"
            if not clean_target.exists():
                candidates = [
                    workdir_path / "r1" / f"{name}_cleaned.pdb",
                    workdir_path / "eq" / f"{name}_cleaned.pdb",
                    workdir_path / "equil" / f"{name}_cleaned.pdb",
                    workdir_path / "prep" / f"{name}_cleaned.pdb",
                    workdir_path.parent / f"{name}_cleaned.pdb",
                    workdir_path.parent / "eq" / f"{name}_cleaned.pdb",
                    workdir_path.parent / "equil" / f"{name}_cleaned.pdb",
                ]
                glob_cands = list(workdir_path.glob("*/*_cleaned.pdb")) + list(workdir_path.glob("*_cleaned.pdb"))
                str_cands = (
                    list(Path.cwd().glob(f"data/structure/**/{name}*/*_cleaned.pdb"))
                    + list(Path.cwd().glob(f"data/structures/**/{name}*/*_cleaned.pdb"))
                    + list(Path.cwd().glob(f"data/str/**/{name}*/*_cleaned.pdb"))
                    + list(Path.cwd().glob(f"structure/**/{name}*/*_cleaned.pdb"))
                    + list(Path.cwd().glob(f"str/**/{name}*/*_cleaned.pdb"))
                )
                for c in glob_cands + str_cands:
                    candidates.append(c)

                for c in candidates:
                    if c.exists() and c.is_file():
                        print(f"[INFO] Auto-discovered cleaned starting structure: {c}")
                        shutil.copy2(c, clean_target)
                        eq_dir = workdir_path / "eq"
                        if eq_dir.exists() and eq_dir.is_dir():
                            eq_clean = eq_dir / f"{name}_cleaned.pdb"
                            if not eq_clean.exists():
                                shutil.copy2(c, eq_clean)
                        r1_dir = workdir_path / "r1"
                        if r1_dir.exists() and r1_dir.is_dir():
                            r1_clean = r1_dir / f"{name}_cleaned.pdb"
                            if not r1_clean.exists():
                                shutil.copy2(c, r1_clean)
                        break

            # Copy ligand file to workdir if provided
            if args.ligand:
                ligand_path = Path(args.ligand)
                local_ligand = Path(workdir_str) / ligand_path.name
                if ligand_path.exists() and ligand_path.resolve() != local_ligand.resolve():
                    shutil.copy2(ligand_path, local_ligand)
                args.ligand = str(local_ligand.resolve())

            # MODULAR RUN: EM -> NVT -> NPT -> Check -> MD
            from colabmda.openmm.modular.utils import derive_replica_seed
            replica_seed = derive_replica_seed(args.seed, args.replica)

            # Skip equilibration for non-r1 replicas (e.g. r2, r3) as they inherit from r1/parent
            is_non_r1_replica = False
            if args.replica and args.replica != "r1":
                is_non_r1_replica = True

            if is_non_r1_replica:
                print(
                    f"[INFO] Replica '{args.replica}' detected (seed={replica_seed}). Skipping EM/NVT/NPT and inheriting equilibrated structure."
                )
            else:
                equil_protocol = _normalize_equil_protocol(args.equil_protocol)
                openmm_em(
                    workdir_str,
                    name,
                    ligand=args.ligand,
                    keep_mg=args.keep_mg,
                    amber_prmtop=args.amber_prmtop,
                    amber_inpcrd=args.amber_inpcrd,
                    padding_nm=args.padding_nm,
                    protein_ff=args.protein_ff,
                    water_model=args.water_model,
                    small_molecule_ff=args.small_molecule_ff,
                )
                openmm_nvt(
                    workdir_str, name, args.equil_time, seed=replica_seed, protocol=equil_protocol
                )
                openmm_npt(
                    workdir_str, name, args.equil_time, seed=replica_seed, protocol=equil_protocol
                )
                openmm_check_equil(workdir_str, strict=equil_protocol == "varmdyn")

            if getattr(args, "equil_only", False):
                print("✔ Shared equilibration complete. Stopping before production MD as requested by --equil-only.")
                return

            openmm_md(
                workdir=workdir_str,
                pdbid=name,
                total_ns=args.total_ns,
                traj_interval=args.traj_interval,
                checkpoint_ps=args.checkpoint_ps,
                sync_dir=args.sync_dir,
                seed=replica_seed,
            )

        elif args.cmd == "merge":
            root = _resolve_root(args.drive, args.root)
            name = getattr(args, "name", None)
            replica = getattr(args, "replica", "r1") or "r1"
            if name:
                sim_base = _simulation_workdir(root or ".", name)
                pdbid_dir = str((sim_base / replica).resolve())
            elif args.pdb_id and root:
                pdbid_dir = str(Path(root) / args.pdb_id / "run")
            elif args.pdb_id:
                pdbid_dir = str(Path(args.pdb_id) / "run")
            elif args.pdb_dir:
                pdb_path = Path(args.pdb_dir)
                if not pdb_path.is_absolute() and root and (Path(root) / pdb_path).exists():
                    pdbid_dir = str((Path(root) / pdb_path).resolve())
                else:
                    pdbid_dir = args.pdb_dir
            else:
                pdbid_dir = "."

            # Resolve output paths: bare filenames go inside pdbid_dir, not CWD
            _pdbid_path = Path(pdbid_dir).resolve()
            out_traj = (
                args.out_traj
                if Path(args.out_traj).is_absolute()
                else str(_pdbid_path / args.out_traj)
            )
            out_log = (
                args.out_log
                if Path(args.out_log).is_absolute()
                else str(_pdbid_path / args.out_log)
            )

            openmm_merge(
                pdbid_dir,
                args.topology,
                out_traj,
                out_log,
                stride=args.stride,
                center=args.center,
                wrap=args.wrap,
                mda=args.mda,
                selection=args.selection,
                ca_only=args.ca_only,
                protein_only=args.protein_only,
            )

        elif args.cmd == "analysis":
            root = _resolve_root(args.drive, args.root)
            name = getattr(args, "name", None)
            replica = getattr(args, "replica", None)
            if name:
                sim_base = _simulation_workdir(root or ".", name)
                sim_dir = str((sim_base / replica).resolve()) if replica else str(sim_base.resolve())
                if args.outdir:
                    out_base = Path(args.outdir)
                elif root:
                    out_base = _analysis_outdir(root, name)
                    if replica:
                        out_base = out_base / replica
                else:
                    out_base = Path("data/analysis") / name
                    if replica:
                        out_base = out_base / replica
            elif args.pdb_id and root:
                sim_path = Path(root) / "data" / "sim" / args.pdb_id
                if not sim_path.exists():
                    sim_path = Path(root) / "sim" / args.pdb_id
                if not sim_path.exists():
                    sim_path = Path(root) / "simulations" / args.pdb_id
                sim_dir = str(sim_path)
                out_base = Path(root) / "data" / "analysis" / args.pdb_id
            elif args.pdb_id:
                sim_path = Path("data") / "sim" / args.pdb_id
                if not sim_path.exists():
                    sim_path = Path("sim") / args.pdb_id
                if not sim_path.exists():
                    sim_path = Path("simulations") / args.pdb_id
                sim_dir = str(sim_path)
                out_base = Path("data") / "analysis" / args.pdb_id
            elif args.pdb_dir:
                pdb_path = Path(args.pdb_dir)
                if not pdb_path.is_absolute() and root and (Path(root) / pdb_path).exists():
                    sim_dir = str((Path(root) / pdb_path).resolve())
                else:
                    sim_dir = args.pdb_dir
                parent_name = Path(sim_dir).parent.name if Path(sim_dir).name.startswith("r") else Path(sim_dir).name
                if args.outdir:
                    out_base = Path(args.outdir)
                elif root:
                    out_base = _analysis_outdir(root, parent_name)
                else:
                    out_base = Path("data/analysis") / parent_name
            else:
                sim_dir = "."
                out_base = Path(args.outdir) if args.outdir else Path("data/analysis/current")

            openmm_analysis(sim_dir, args.topology, args.trajectory, args.interval, str(out_base))

        elif args.cmd == "compare":
            openmm_compare(args.series, args.outdir)

        elif args.cmd == "view":
            openmm_view(args.pdb_dir, args.topology, args.trajectory, resi=args.resi)

        elif args.cmd == "snapshots":
            openmm_snapshots(
                config_path=args.config,
                output=args.output,
                temp_dir=args.temp_dir,
                font_path=args.font_path,
            )

        elif args.cmd == "status":
            root = _resolve_root(args.drive, args.root)
            name = getattr(args, "name", None)
            replica = getattr(args, "replica", "r1") or "r1"
            if name:
                sim_base = _simulation_workdir(root or ".", name)
                pdbid_dir = str((sim_base / replica).resolve())
            elif args.pdb_id and root:
                pdbid_dir = str(Path(root) / args.pdb_id / "run")
            elif args.pdb_id:
                pdbid_dir = str(Path(args.pdb_id) / "run")
            elif args.pdb_dir:
                pdb_path = Path(args.pdb_dir)
                if not pdb_path.is_absolute() and root and (Path(root) / pdb_path).exists():
                    pdbid_dir = str((Path(root) / pdb_path).resolve())
                else:
                    pdbid_dir = args.pdb_dir
            else:
                pdbid_dir = "."
            openmm_status(pdbid_dir)

        elif args.cmd in ["em", "nvt", "npt", "qc", "check-equil", "md"]:
            name = args.name
            if args.workdir:
                workdir = str(Path(args.workdir).resolve())
            else:
                root_val = args.root or _resolve_root(getattr(args, "drive", False), None)
                sim_base = _simulation_workdir(root_val, name)
                if args.cmd in ["em", "nvt", "npt", "qc", "check-equil"]:
                    if (sim_base / "eq").exists():
                        workdir_path = sim_base / "eq"
                    elif (sim_base / "equil").exists():
                        workdir_path = sim_base / "equil"
                    else:
                        workdir_path = sim_base / "eq"
                elif args.cmd == "md":
                    rep = getattr(args, "replica", "r1") or "r1"
                    workdir_path = sim_base / rep
                else:
                    workdir_path = sim_base
                workdir_path.mkdir(parents=True, exist_ok=True)
                workdir = str(workdir_path.resolve())

                # If cleaned PDB exists in sim_base, eq, or equil, ensure it's available in workdir
                workdir_clean = workdir_path / f"{name}_cleaned.pdb"
                if not workdir_clean.exists():
                    for src in [
                        sim_base / f"{name}_cleaned.pdb",
                        sim_base / "eq" / f"{name}_cleaned.pdb",
                        sim_base / "equil" / f"{name}_cleaned.pdb",
                    ]:
                        if src.exists():
                            shutil.copy2(src, workdir_clean)
                            break
            if args.cmd == "em":
                if args.ligand:
                    ligand_path = Path(args.ligand)
                    local_ligand = Path(workdir) / ligand_path.name
                    if ligand_path.exists() and ligand_path.resolve() != local_ligand.resolve():
                        shutil.copy2(ligand_path, local_ligand)
                    args.ligand = str(local_ligand.resolve())
                openmm_em(
                    workdir,
                    name,
                    ligand=args.ligand,
                    keep_mg=args.keep_mg,
                    amber_prmtop=args.amber_prmtop,
                    amber_inpcrd=args.amber_inpcrd,
                    padding_nm=args.padding_nm,
                    protein_ff=args.protein_ff,
                    water_model=args.water_model,
                    small_molecule_ff=args.small_molecule_ff,
                )
            elif args.cmd == "nvt":
                openmm_nvt(
                    workdir,
                    name,
                    args.equil_time,
                    args.seed,
                    _normalize_equil_protocol(args.equil_protocol),
                )
            elif args.cmd == "npt":
                openmm_npt(
                    workdir,
                    name,
                    args.equil_time,
                    seed=args.seed,
                    protocol=_normalize_equil_protocol(args.equil_protocol),
                )
            elif args.cmd in ["qc", "check-equil"]:
                openmm_qc(workdir, strict=not args.warn_only)
            elif args.cmd == "md":
                openmm_md(
                    workdir,
                    name,
                    args.total_ns,
                    args.traj_interval,
                    args.checkpoint_ps,
                    args.sync_dir,
                    seed=args.seed,
                )

    elif args.tool == "modeller":
        from colabmda.modeller.commands import (
            modeller_build,
            modeller_mutate,
        )

        if args.cmd == "build":
            modeller_build(args)
        elif args.cmd == "mutate":
            if not (args.mut or args.list):
                raise SystemExit("ERROR: mutate requires --mut or --list")
            modeller_mutate(args)


if __name__ == "__main__":
    main()
