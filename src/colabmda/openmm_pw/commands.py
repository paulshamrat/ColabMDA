#!/usr/bin/env python3
import os
import sys
import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]   # .../src/colabmda/openmm_pw/commands.py -> repo root

SCRIPTS = {
    # "old" workflow scripts (pdb-id download)
    "clean_by_pdbid": REPO_ROOT / "openmm" / "openmm_proteinwater" / "pdbfixer_cleaning.py",
    "run_basic":      REPO_ROOT / "openmm" / "openmm_proteinwater" / "openmm_proteinwater.py",
    "merge":          REPO_ROOT / "openmm" / "openmm_proteinwater" / "openmm_trajmerge.py",
    "analysis":       REPO_ROOT / "openmm" / "openmm_proteinwater" / "openmm_trajanalysis.py",

    # "colab-safe" workflow scripts (local pdb file cleaning + robust resume)
    "clean_from_file": REPO_ROOT / "openmm" / "openmm_proteinwater_260203" / "pdbfixer_clean_fromfile.py",
    "run_colab":       REPO_ROOT / "openmm" / "openmm_proteinwater_260203" / "openmm_proteinwater_colab.py",
}

def _py():
    return sys.executable

def _run(script: Path, argv: list[str]):
    if not script.exists():
        raise SystemExit(f"ERROR: expected script not found:\n  {script}\n"
                         f"Are you running from the ColabMDA repo root?")
    cmd = [_py(), str(script)] + argv
    print("\n[RUN]", " ".join(cmd), "\n")
    raise SystemExit(subprocess.call(cmd))

def openmm_prep_from_pdbid(pdbid: str):
    _run(SCRIPTS["clean_by_pdbid"], [pdbid])

def openmm_prep_from_file(pdb_file: str, outdir: str, pdbid: str = "4ldj", ph: float = 7.0):
    _run(SCRIPTS["clean_from_file"], ["--in", pdb_file, "--outdir", outdir, "--pdbid", pdbid, "--ph", str(ph)])

def openmm_run_colab(workdir: str, pdbid: str, total_ns: float, traj_interval: float,
                     equil_time: float, checkpoint_ps: float, sync_dir: str | None):
    argv = [
        workdir,
        "--pdbid", pdbid,
        "--total-ns", str(total_ns),
        "--traj-interval", str(traj_interval),
        "--equil-time", str(equil_time),
        "--checkpoint-ps", str(checkpoint_ps),
    ]
    if sync_dir:
        argv += ["--sync-dir", sync_dir]
    _run(SCRIPTS["run_colab"], argv)

def openmm_merge(pdbid_dir: str, topology: str | None, out_traj: str, out_log: str):
    argv = [pdbid_dir, "--out-traj", out_traj, "--out-log", out_log]
    if topology:
        argv += ["--topology", topology]
    _run(SCRIPTS["merge"], argv)

def openmm_analysis(pdbid_dir: str, topology: str | None, trajectory: str | None, interval_ps: float | None, outdir: str | None):
    argv = [pdbid_dir]
    if topology:
        argv += ["--topology", topology]
    if trajectory:
        argv += ["--trajectory", trajectory]
    if interval_ps is not None:
        argv += ["--interval", str(interval_ps)]
    if outdir:
        argv += ["--outdir", outdir]
    _run(SCRIPTS["analysis"], argv)
