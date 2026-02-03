#!/usr/bin/env python3
import sys
import subprocess
from pathlib import Path
from importlib import resources

def _py():
    return sys.executable

def _run(argv: list[str]):
    script = resources.files("colabmda.legacy.modeller").joinpath("modeller6.py")
    if not script.exists():
        raise SystemExit(
            "ERROR: bundled modeller6.py not found.\n"
            "This indicates an incomplete install. Reinstall with:\n"
            "  pip install -e .\n"
        )
    cmd = [_py(), str(script)] + argv
    print("\n[RUN]", " ".join(cmd), "\n")
    raise SystemExit(subprocess.call(cmd))

def modeller_build(args):
    argv = [args.pdb_id, args.uniprot_id]
    if args.chain: argv += ["--chain", args.chain]
    if args.range: argv += ["--range", str(args.range[0]), str(args.range[1])]
    if args.truncate: argv += ["--truncate"]
    if args.mut: argv += ["--mut", args.mut]
    if args.list: argv += ["--list", args.list]
    if args.outdir: argv += ["--outdir", args.outdir]
    if args.outdir_mut: argv += ["--outdir-mut", args.outdir_mut]
    if args.seed is not None: argv += ["--seed", str(args.seed)]
    if args.logfile: argv += ["--logfile", args.logfile]
    if args.verbose: argv += ["--verbose"]
    _run(argv)

def modeller_mutate(args):
    argv = ["--pdb-in", args.pdb_in, "--chain", args.chain]
    if args.mut: argv += ["--mut", args.mut]
    if args.list: argv += ["--list", args.list]
    if args.outdir_mut: argv += ["--outdir-mut", args.outdir_mut]
    if args.seed is not None: argv += ["--seed", str(args.seed)]
    if args.logfile: argv += ["--logfile", args.logfile]
    if args.verbose: argv += ["--verbose"]
    _run(argv)
