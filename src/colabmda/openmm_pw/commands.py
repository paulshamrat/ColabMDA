#!/usr/bin/env python3
import os
import sys
import subprocess
from pathlib import Path
from importlib import resources
import re

SCRIPTS = {
    # Bundled workflow scripts (pdb-id download)
    "clean_by_pdbid": ("colabmda.legacy.openmm_proteinwater", "pdbfixer_cleaning.py"),
    "run_basic":      ("colabmda.legacy.openmm_proteinwater", "openmm_proteinwater.py"),
    "merge":          ("colabmda.legacy.openmm_proteinwater", "openmm_trajmerge.py"),
    "analysis":       ("colabmda.legacy.openmm_proteinwater", "openmm_trajanalysis.py"),

    # Bundled colab-safe workflow scripts (local pdb file cleaning + robust resume)
    "clean_from_file": ("colabmda.legacy.openmm_proteinwater_260203", "pdbfixer_clean_fromfile.py"),
    "run_colab":       ("colabmda.legacy.openmm_proteinwater_260203", "openmm_proteinwater_colab.py"),
}

def _py():
    return sys.executable

def _script_path(pkg: str, name: str) -> Path:
    return resources.files(pkg).joinpath(name)

def _iter_data_lines(path: Path):
    if not path.exists():
        return
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s[0].isdigit():
                yield s

def _parse_last_step_time(log_path: Path):
    last = None
    for s in _iter_data_lines(log_path):
        last = s
    if not last:
        return None, None
    parts = re.split(r"\s+", last)
    if len(parts) < 2:
        return None, None
    try:
        step = int(float(parts[0]))
        time_ps = float(parts[1])
        return step, time_ps
    except Exception:
        return None, None

def _parse_chunk_ranges(workdir: Path):
    ranges = []
    rx = re.compile(r"^prod_(\d+)to(\d+)ps\.(?:dcd|log)$")
    for p in workdir.iterdir():
        m = rx.match(p.name)
        if m:
            start = int(m.group(1))
            end = int(m.group(2))
            ranges.append((start, end))
    # de-dup and sort
    ranges = sorted(set(ranges))
    return ranges

def _range_gaps(ranges):
    gaps = []
    if not ranges:
        return gaps
    prev_end = ranges[0][1]
    for start, end in ranges[1:]:
        if start > prev_end:
            gaps.append((prev_end, start))
        prev_end = max(prev_end, end)
    return gaps

def openmm_status(pdbid_dir: str):
    workdir = Path(pdbid_dir).resolve()
    if not workdir.exists():
        raise SystemExit(f"ERROR: directory not found: {workdir}")

    logs = sorted(workdir.glob("prod_*to*ps.log"))
    dcds = sorted(workdir.glob("prod_*to*ps.dcd"))
    merged_dcd = workdir / "prod_full.dcd"
    merged_log = workdir / "prod_full.log"

    ranges = _parse_chunk_ranges(workdir)
    gaps = _range_gaps(ranges)

    total_frames = 0
    max_step = None
    max_time_ps = None
    for lp in logs:
        # Count frames by data lines
        frames = sum(1 for _ in _iter_data_lines(lp))
        total_frames += frames
        step, time_ps = _parse_last_step_time(lp)
        if step is not None and (max_step is None or step > max_step):
            max_step = step
        if time_ps is not None and (max_time_ps is None or time_ps > max_time_ps):
            max_time_ps = time_ps

    chk = workdir / "prod.chk"
    xml = workdir / "system.xml"
    solv = workdir / "solvated.pdb"
    can_resume = chk.exists() and xml.exists() and solv.exists()

    print("\n[STATUS]")
    print(f"  Workdir          : {workdir}")
    print(f"  Chunks (DCD/log) : {len(dcds)} / {len(logs)}")
    if total_frames:
        print(f"  Frames (from logs): {total_frames}")
    else:
        print("  Frames           : (no chunk logs found)")

    if ranges:
        max_end_ps = max(end for _, end in ranges)
        ns = max_end_ps / 1000.0
        print(f"  Sim time (ps/ns) : {max_end_ps:.2f} ps / {ns:.4f} ns (from chunk names)")
    elif max_time_ps is not None:
        ns = max_time_ps / 1000.0
        print(f"  Sim time (ps/ns) : {max_time_ps:.2f} ps / {ns:.4f} ns (from logs)")
    else:
        print("  Sim time         : (unknown)")

    if gaps:
        gap_str = ", ".join([f"{a}to{b}ps" for a, b in gaps])
        print(f"  Gaps detected    : {gap_str}")

    if max_step is not None:
        print(f"  Last step        : {max_step}")

    print(f"  Resume-ready     : {'YES' if can_resume else 'NO'} (needs prod.chk, system.xml, solvated.pdb)")
    print(f"  Merged DCD       : {'YES' if merged_dcd.exists() else 'NO'}")
    print(f"  Merged log       : {'YES' if merged_log.exists() else 'NO'}")
def _run(script: Path, argv: list[str], cwd: str | None = None):
    if not script.exists():
        raise SystemExit(
            "ERROR: expected bundled script not found:\n"
            f"  {script}\n"
            "This indicates an incomplete install. Reinstall with:\n"
            "  pip install -e .\n"
        )
    cmd = [_py(), str(script)] + argv
    print("\n[RUN]", " ".join(cmd), "\n")
    raise SystemExit(subprocess.call(cmd, cwd=cwd))

def openmm_prep_from_pdbid(pdbid: str, root_dir: str | None = None):
    pkg, name = SCRIPTS["clean_by_pdbid"]
    base = Path(root_dir or os.getcwd()).resolve()
    outdir = base / pdbid
    print(f"[INFO] Prep output will be written to: {outdir}")
    _run(_script_path(pkg, name), [pdbid], cwd=root_dir)

def openmm_prep_from_file(pdb_file: str, outdir: str, pdbid: str = "4ldj", ph: float = 7.0):
    pkg, name = SCRIPTS["clean_from_file"]
    print(f"[INFO] Prep output will be written to: {Path(outdir).resolve()}")
    _run(_script_path(pkg, name), ["--in", pdb_file, "--outdir", outdir, "--pdbid", pdbid, "--ph", str(ph)])

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
    pkg, name = SCRIPTS["run_colab"]
    _run(_script_path(pkg, name), argv)

def openmm_merge(pdbid_dir: str, topology: str | None, out_traj: str, out_log: str):
    argv = [pdbid_dir, "--out-traj", out_traj, "--out-log", out_log]
    if topology:
        argv += ["--topology", topology]
    pkg, name = SCRIPTS["merge"]
    _run(_script_path(pkg, name), argv)

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
    pkg, name = SCRIPTS["analysis"]
    _run(_script_path(pkg, name), argv)
