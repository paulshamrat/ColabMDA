#!/usr/bin/env python3
import argparse
from pathlib import Path

from colabmda.openmm_pw.commands import (
    openmm_prep_from_pdbid,
    openmm_prep_from_file,
    openmm_run_colab,
    openmm_merge,
    openmm_analysis,
)
from colabmda.modeller.commands import (
    modeller_build,
    modeller_mutate,
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

    # run (colab-safe runner)
    p_run = sub_openmm.add_parser("run", help="Run/resume chunked MD (colab-safe)")
    g = p_run.add_mutually_exclusive_group(required=True)
    g.add_argument("--pdb-id", help="Use ./<pdb-id> as workdir (pdb-id download workflow)")
    g.add_argument("--workdir", help="Folder containing <name>_cleaned.pdb")
    p_run.add_argument("--name", default=None, help="Prefix name (default: --pdb-id or inferred from workdir)")
    p_run.add_argument("--total-ns", type=float, default=100.0)
    p_run.add_argument("--traj-interval", type=float, default=100.0, help="ps between saved frames")
    p_run.add_argument("--equil-time", type=float, default=100.0, help="ps for NVT and ps for NPT")
    p_run.add_argument("--checkpoint-ps", type=float, default=1000.0, help="ps per chunk")
    p_run.add_argument("--sync-dir", default=None, help="Optional: sync outputs to this directory")

    # merge
    p_merge = sub_openmm.add_parser("merge", help="Merge chunk DCDs/logs")
    g = p_merge.add_mutually_exclusive_group(required=True)
    g.add_argument("--pdb-id", help="Use ./<pdb-id> as simulation directory")
    g.add_argument("--pdb-dir", help="Simulation directory (e.g. 4ldj_wt)")
    p_merge.add_argument("--topology", default=None, help="Topology PDB (default: <dir>/solvated.pdb)")
    p_merge.add_argument("--out-traj", default="prod_full.dcd")
    p_merge.add_argument("--out-log", default="prod_full.log")

    # analysis
    p_ana = sub_openmm.add_parser("analysis", help="RMSD/Rg/RMSF analysis")
    g = p_ana.add_mutually_exclusive_group(required=True)
    g.add_argument("--pdb-id", help="Use ./<pdb-id> as simulation directory")
    g.add_argument("--pdb-dir", help="Simulation directory (e.g. 4ldj_wt)")
    p_ana.add_argument("--topology", default=None)
    p_ana.add_argument("--trajectory", default=None)
    p_ana.add_argument("--interval", type=float, default=None, help="ps per frame (if not auto-detected)")
    p_ana.add_argument("--outdir", default=None)

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
                openmm_prep_from_pdbid(args.pdb_id)
            else:
                name = args.name or Path(args.pdb_file).stem
                outdir = args.outdir or name
                openmm_prep_from_file(args.pdb_file, outdir, pdbid=name, ph=args.ph)

        elif args.cmd == "run":
            if args.pdb_id:
                workdir = args.pdb_id
                name = args.name or args.pdb_id
            else:
                workdir = args.workdir
                name = args.name or _guess_pdbid_from_workdir(workdir)
                if not name:
                    raise SystemExit("ERROR: could not infer pdbid from workdir; please specify --name.")
            openmm_run_colab(
                workdir=workdir,
                pdbid=name,
                total_ns=args.total_ns,
                traj_interval=args.traj_interval,
                equil_time=args.equil_time,
                checkpoint_ps=args.checkpoint_ps,
                sync_dir=args.sync_dir,
            )

        elif args.cmd == "merge":
            pdbid_dir = args.pdb_id or args.pdb_dir
            openmm_merge(pdbid_dir, args.topology, args.out_traj, args.out_log)

        elif args.cmd == "analysis":
            pdbid_dir = args.pdb_id or args.pdb_dir
            openmm_analysis(pdbid_dir, args.topology, args.trajectory, args.interval, args.outdir)

    elif args.tool == "modeller":
        if args.cmd == "build":
            modeller_build(args)
        elif args.cmd == "mutate":
            if not (args.mut or args.list):
                raise SystemExit("ERROR: mutate requires --mut or --list")
            modeller_mutate(args)

if __name__ == "__main__":
    main()
