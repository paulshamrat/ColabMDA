#!/usr/bin/env python3
import argparse
import os
import shutil
import datetime
from pathlib import Path

from colabmda.openmm_pw.commands import (
    openmm_prep_from_pdbid,
    openmm_prep_from_file,
    openmm_run_colab,
    openmm_merge,
    openmm_analysis,
    openmm_status,
    openmm_em,
    openmm_nvt,
    openmm_npt,
    openmm_check_equil,
    openmm_md,
)
from colabmda.modeller.commands import (
    modeller_build,
    modeller_mutate,
)

DEFAULT_DRIVE_ROOT = "/content/drive/MyDrive/openmm"
ENV_ROOT = "COLABMDA_ROOT"

def _resolve_root(use_drive: bool, root: str | None) -> str | None:
    if root:
        return root
    env_root = os.environ.get(ENV_ROOT)
    if env_root:
        return env_root
    # Default to Drive-first layout for persistent Colab workflows.
    return DEFAULT_DRIVE_ROOT

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
        raise SystemExit(f"ERROR: multiple *_cleaned.pdb files found in {workdir}: {names}\n"
                         f"Please specify --name.")
    return candidates[0].name.replace("_cleaned.pdb", "")

def _default_project_root() -> str:
    return _resolve_root(use_drive=True, root=None) or DEFAULT_DRIVE_ROOT

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
    p_prep.add_argument("--outdir", default=None, help="Output directory for --pdb-file (default: ./<name>)")
    p_prep.add_argument("--name", default=None, help="Prefix name (default: from --pdb-id or file stem)")
    p_prep.add_argument("--ph", type=float, default=7.0, help="Hydrogen pH (default: 7.0)")
    p_prep.add_argument("--sync-dir", default=None, help="Optional: copy cleaned prep outputs to this directory")
    p_prep.add_argument("--drive", action="store_true", help="(compat) Use Drive root (default behavior)")
    p_prep.add_argument("--root", default=None, help=f"Override base directory (default: ${ENV_ROOT} if set, else {DEFAULT_DRIVE_ROOT} when --drive)")

    # run (colab-safe runner)
    p_run = sub_openmm.add_parser("run", help="Run/resume chunked MD (colab-safe)")
    g = p_run.add_mutually_exclusive_group(required=False)
    g.add_argument("--pdb-id", help="Use ./<pdb-id> as workdir (pdb-id download workflow)")
    g.add_argument("--workdir", help="Folder containing <name>_cleaned.pdb")
    p_run.add_argument("--name", default=None, help="Prefix name (default: --pdb-id or inferred from workdir)")
    p_run.add_argument("--total-ns", type=float, default=100.0)
    p_run.add_argument("--traj-interval", type=float, default=100.0, help="ps between saved frames")
    p_run.add_argument("--equil-time", type=float, default=100.0, help="ps for NVT and ps for NPT")
    p_run.add_argument("--checkpoint-ps", type=float, default=1000.0, help="ps per chunk")
    p_run.add_argument("--sync-dir", default=None, help="Optional: sync outputs to this directory")
    p_run.add_argument("--replica", default=None, help="Optional: replica subfolder (e.g. r1, r2)")
    p_run.add_argument("--seed", type=int, default=None, help="Optional: random seed for velocity assignment")
    p_run.add_argument("--drive", action="store_true", help="(compat) Use Drive root (default behavior)")
    p_run.add_argument("--root", default=None, help=f"Override base directory (default: ${ENV_ROOT} if set, else {DEFAULT_DRIVE_ROOT} when --drive)")

    # merge
    p_merge = sub_openmm.add_parser("merge", help="Merge chunk DCDs/logs")
    g = p_merge.add_mutually_exclusive_group(required=False)
    g.add_argument("--pdb-id", help="Use ./<pdb-id> as simulation directory")
    g.add_argument("--pdb-dir", help="Simulation directory (e.g. 4ldj_wt)")
    p_merge.add_argument("--topology", default=None, help="Topology PDB (default: <dir>/solvated.pdb)")
    p_merge.add_argument("--out-traj", default="prod_full.dcd")
    p_merge.add_argument("--out-log", default="prod_full.log")
    p_merge.add_argument("--stride", type=int, default=1, help="Keep every Nth frame while merging (default: 1)")
    p_merge.add_argument("--center", action="store_true", help="Center protein in the box")
    p_merge.add_argument("--wrap", action="store_true", help="Wrap solvent molecules (image_molecules)")
    p_merge.add_argument("--drive", action="store_true", help="(compat) Use Drive root (default behavior)")
    p_merge.add_argument("--root", default=None, help=f"Override base directory (default: ${ENV_ROOT} if set, else {DEFAULT_DRIVE_ROOT} when --drive)")

    # analysis
    p_ana = sub_openmm.add_parser("analysis", help="RMSD/Rg/RMSF analysis")
    g = p_ana.add_mutually_exclusive_group(required=False)
    g.add_argument("--pdb-id", help="Use ./<pdb-id> as simulation directory")
    g.add_argument("--pdb-dir", help="Simulation directory (e.g. 4ldj_wt)")
    p_ana.add_argument("--topology", default=None)
    p_ana.add_argument("--trajectory", default=None)
    p_ana.add_argument("--interval", type=float, default=None, help="ps per frame (if not auto-detected)")
    p_ana.add_argument("--outdir", default=None)
    p_ana.add_argument("--drive", action="store_true", help="(compat) Use Drive root (default behavior)")
    p_ana.add_argument("--root", default=None, help=f"Override base directory (default: ${ENV_ROOT} if set, else {DEFAULT_DRIVE_ROOT} when --drive)")

    # status
    p_stat = sub_openmm.add_parser("status", help="Sanity-check frames/time/resume readiness")
    g = p_stat.add_mutually_exclusive_group(required=False)
    g.add_argument("--pdb-id", help="Use ./<pdb-id> as simulation directory")
    g.add_argument("--pdb-dir", help="Simulation directory (e.g. 4ldj_wt)")
    p_stat.add_argument("--drive", action="store_true", help="(compat) Use Drive root (default behavior)")
    p_stat.add_argument("--root", default=None, help=f"Override base directory (default: ${ENV_ROOT} if set, else {DEFAULT_DRIVE_ROOT} when --drive)")

    # stage
    p_stage = sub_openmm.add_parser("stage", help="Stage a WT/mutant structure into simulations/<name>")
    p_stage.add_argument("--pdb-file", required=True, help="Input structure PDB file (typically from structures/)")
    p_stage.add_argument("--name", required=True, help="Simulation name (e.g. 4ldj_wt, 4ldj_G12C)")
    p_stage.add_argument("--replica", default=None, help="Optional: create nested replica subfolder (e.g. r1, r2)")
    p_stage.add_argument("--ph", type=float, default=7.0, help="Hydrogen pH (default: 7.0)")
    p_stage.add_argument("--root", default=None, help=f"Project root (default: ${ENV_ROOT} or {DEFAULT_DRIVE_ROOT})")

    # em/nvt/npt/check-equil/md (individual modular steps)
    p_em = sub_openmm.add_parser("em", help="Modular: Minimization")
    p_em.add_argument("--name", required=True)
    p_em.add_argument("--workdir", default=None)
    p_em.add_argument("--root", default=None)

    p_nvt = sub_openmm.add_parser("nvt", help="Modular: NVT Equilibration")
    p_nvt.add_argument("--name", required=True)
    p_nvt.add_argument("--equil-time", type=float, default=100.0)
    p_nvt.add_argument("--seed", type=int, default=None)
    p_nvt.add_argument("--workdir", default=None)
    p_nvt.add_argument("--root", default=None)

    p_npt = sub_openmm.add_parser("npt", help="Modular: NPT Equilibration")
    p_npt.add_argument("--name", required=True)
    p_npt.add_argument("--equil-time", type=float, default=100.0)
    p_npt.add_argument("--workdir", default=None)
    p_npt.add_argument("--root", default=None)

    p_chk = sub_openmm.add_parser("check-equil", help="Modular: Stability Check & QC Plots")
    p_chk.add_argument("--name", required=True)
    p_chk.add_argument("--workdir", default=None)
    p_chk.add_argument("--root", default=None)

    p_md = sub_openmm.add_parser("md", help="Modular: Production MD")
    p_md.add_argument("--name", required=True)
    p_md.add_argument("--total-ns", type=float, default=100.0)
    p_md.add_argument("--traj-interval", type=float, default=100.0)
    p_md.add_argument("--checkpoint-ps", type=float, default=1000.0)
    p_md.add_argument("--sync-dir", default=None)
    p_md.add_argument("--workdir", default=None)
    p_md.add_argument("--root", default=None)

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
        if args.cmd == "prep":
            if args.pdb_id:
                root = _resolve_root(args.drive, args.root)
                if root:
                    _ensure_dir(root)
                openmm_prep_from_pdbid(args.pdb_id, root_dir=root, sync_dir=args.sync_dir)
            else:
                name = args.name or Path(args.pdb_file).stem
                root = _resolve_root(args.drive, args.root)
                if root:
                    _ensure_dir(root)
                outdir = args.outdir or (str(Path(root) / name / "prep") if root else name)
                openmm_prep_from_file(args.pdb_file, outdir, pdbid=name, ph=args.ph, sync_dir=args.sync_dir)

        elif args.cmd == "run":
            if args.pdb_id:
                root = _resolve_root(args.drive, args.root)
                if root:
                    _ensure_dir(root)
                workdir = str(Path(root) / args.pdb_id / "run") if root else str(Path(args.pdb_id) / "run")
                name = args.name or args.pdb_id
                if root:
                    _prepare_run_inputs(root=root, pdbid=args.pdb_id, name=name)
            elif args.workdir:
                workdir = str(Path(args.workdir).resolve())
                name = args.name or _guess_pdbid_from_workdir(workdir)
            else:
                root = _resolve_root(args.drive, args.root)
                name = args.name
                if name and root and (Path(root) / "simulations" / name).exists():
                    workdir = str(Path(root) / "simulations" / name)
                elif name and (Path.cwd() / "simulations" / name).exists():
                    workdir = str(Path.cwd() / "simulations" / name)
                else:
                    workdir = os.getcwd()
                
                if not name:
                    name = _guess_pdbid_from_workdir(workdir)
            
            if args.replica:
                # If we are already inside the replica folder, don't append it again
                if Path(workdir).name == args.replica:
                    pass
                else:
                    workdir = str(Path(workdir) / args.replica)
                _ensure_dir(workdir)
                # If we are in a replica folder, we need the *_cleaned.pdb from the parent or here
                local_clean = Path(workdir) / f"{name}_cleaned.pdb"
                parent_clean = Path(workdir).parent / f"{name}_cleaned.pdb"
                if not local_clean.exists() and parent_clean.exists():
                    shutil.copy2(parent_clean, local_clean)
                elif not local_clean.exists():
                    # Check if ANY _cleaned.pdb exists here or parent to infer name
                    guesses = list(Path(workdir).glob("*_cleaned.pdb")) or list(Path(workdir).parent.glob("*_cleaned.pdb"))
                    if guesses and not name:
                        name = guesses[0].name.replace("_cleaned.pdb", "")
                        shutil.copy2(guesses[0], Path(workdir) / guesses[0].name)

            if not name:
                raise SystemExit("ERROR: could not infer pdbid; please specify --name.")

            # MODULAR RUN: EM -> NVT -> NPT -> Check -> MD
            openmm_em(workdir, name)
            openmm_nvt(workdir, name, args.equil_time, seed=args.seed)
            openmm_npt(workdir, name, args.equil_time)
            openmm_check_equil(workdir)
            openmm_md(
                workdir=workdir,
                pdbid=name,
                total_ns=args.total_ns,
                traj_interval=args.traj_interval,
                checkpoint_ps=args.checkpoint_ps,
                sync_dir=args.sync_dir,
            )

        elif args.cmd == "merge":
            root = _resolve_root(args.drive, args.root)
            if args.pdb_id and root:
                pdbid_dir = str(Path(root) / args.pdb_id / "run")
            elif args.pdb_id:
                pdbid_dir = str(Path(args.pdb_id) / "run")
            elif args.pdb_dir:
                pdbid_dir = args.pdb_id or args.pdb_dir
            else:
                pdbid_dir = "."
            openmm_merge(pdbid_dir, args.topology, args.out_traj, args.out_log, stride=args.stride, center=args.center, wrap=args.wrap)

        elif args.cmd == "analysis":
            root = _resolve_root(args.drive, args.root)
            if args.pdb_id and root:
                sim_dir = str(Path(root) / "simulations" / args.pdb_id)
                out_base = Path(root) / "analysis" / "single" / args.pdb_id
            elif args.pdb_id:
                sim_dir = str(Path("simulations") / args.pdb_id)
                out_base = Path("analysis") / "single" / args.pdb_id
            elif args.pdb_dir:
                sim_dir = args.pdb_dir
                out_base = Path(args.outdir) if args.outdir else Path("analysis/single") / Path(sim_dir).name
            else:
                sim_dir = "."
                out_base = Path(args.outdir) if args.outdir else Path("analysis/single/current")

            openmm_analysis(sim_dir, args.topology, args.trajectory, args.interval, str(out_base))

        elif args.cmd == "status":
            root = _resolve_root(args.drive, args.root)
            if args.pdb_id and root:
                pdbid_dir = str(Path(root) / args.pdb_id / "run")
            elif args.pdb_id:
                pdbid_dir = str(Path(args.pdb_id) / "run")
            elif args.pdb_dir:
                pdbid_dir = args.pdb_id or args.pdb_dir
            else:
                pdbid_dir = "."
            openmm_status(pdbid_dir)

        elif args.cmd == "stage":
            root = args.root or _default_project_root()
            _ensure_dir(root)
            outdir = Path(root) / "simulations" / args.name
            if args.replica:
                outdir = outdir / args.replica
            _ensure_dir(str(outdir))
            openmm_prep_from_file(args.pdb_file, str(outdir), pdbid=args.name, ph=args.ph)
            print(f"[INFO] Staged simulation folder: {outdir}")

        elif args.cmd in ["em", "nvt", "npt", "check-equil", "md"]:
            root = _resolve_root(False, args.root)
            workdir = args.workdir or os.getcwd()
            name = args.name
            if args.cmd == "em":
                openmm_em(workdir, name)
            elif args.cmd == "nvt":
                openmm_nvt(workdir, name, args.equil_time, args.seed)
            elif args.cmd == "npt":
                openmm_npt(workdir, name, args.equil_time)
            elif args.cmd == "check-equil":
                openmm_check_equil(workdir)
            elif args.cmd == "md":
                openmm_md(workdir, name, args.total_ns, args.traj_interval, args.checkpoint_ps, args.sync_dir)

    elif args.tool == "modeller":
        if args.cmd == "build":
            modeller_build(args)
        elif args.cmd == "mutate":
            if not (args.mut or args.list):
                raise SystemExit("ERROR: mutate requires --mut or --list")
            modeller_mutate(args)

if __name__ == "__main__":
    main()
