#!/usr/bin/env python3
import os
import re
import shutil
import subprocess
import sys
from importlib import resources
from pathlib import Path

# Imports relocated internally to functions to avoid module-level import errors
# when running visualization subcommands in environments without simulation packages.

SCRIPTS = {
    "clean_by_pdbid": ("colabmda.openmm", "prepare.py"),
    "clean_from_file": ("colabmda.openmm", "prepare.py"),
    "run_colab": ("colabmda.openmm", "run.py"),
    "merge": ("colabmda.openmm", "merge.py"),
    "analysis": ("colabmda.openmm", "analysis.py"),
    **{stage: ("colabmda.openmm", f"{stage}.py") for stage in ("em", "nvt", "npt", "qc", "md")},
    "check_equil": ("colabmda.openmm", "qc.py"),
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
    parts = last.split(",") if "," in last else re.split(r"\s+", last)
    if len(parts) < 2:
        return None, None
    try:
        step = int(float(parts[0].strip()))
        time_ps = float(parts[1].strip())
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


def _read_dcd_frames(dcd_path: Path) -> int | None:
    try:
        with open(dcd_path, "rb") as f:
            header = f.read(12)
            if len(header) < 12:
                return None
            import struct

            block_size, signature, nset = struct.unpack("i4sI", header[:12])
            if block_size != 84:
                block_size, signature, nset = struct.unpack(">i4sI", header[:12])
            if signature == b"CORD":
                return nset
    except Exception:
        pass
    return None


def _read_topology_info(pdb_path: Path) -> dict | None:
    if not pdb_path.is_file():
        return None
    try:
        atoms = 0
        residues = set()
        waters = 0
        ions = 0
        ion_names = {"NA", "CL", "SOD", "CLA", "K", "MG", "CA", "ZN"}
        water_names = {"HOH", "WAT", "SOL", "TIP3"}

        with open(pdb_path) as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    atoms += 1
                    res_name = line[17:20].strip()
                    res_seq = line[22:26].strip()
                    chain_id = line[21].strip()
                    res_key = (chain_id, res_seq, res_name)
                    if res_key not in residues:
                        residues.add(res_key)
                        if res_name in water_names:
                            waters += 1
                        elif res_name.upper() in ion_names:
                            ions += 1

        protein_res = len(residues) - waters - ions
        return {
            "atoms": atoms,
            "residues": len(residues),
            "protein_residues": protein_res,
            "waters": waters,
            "ions": ions,
        }
    except Exception:
        pass
    return None


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

    topo_pdb = (
        solv
        if solv.exists()
        else (merged_dcd.with_suffix(".pdb") if merged_dcd.with_suffix(".pdb").exists() else None)
    )

    print("\n[STATUS]")
    print(f"  Workdir          : {workdir}")
    print(f"  Chunks (DCD/log) : {len(dcds)} / {len(logs)}")
    if topo_pdb:
        t_info = _read_topology_info(topo_pdb)
        if t_info:
            print(f"  Topology File    : {topo_pdb.name}")
            print(f"                     └─ {t_info['atoms']} atoms, {t_info['residues']} residues")
            print(
                f"                     └─ ({t_info['protein_residues']} protein, {t_info['waters']} water, {t_info['ions']} ions)"
            )
    else:
        print("  Topology File    : NOT FOUND")

    merged_dcd_str = "NOT FOUND"
    if merged_dcd.exists():
        frames = _read_dcd_frames(merged_dcd)
        if frames is not None:
            merged_dcd_str = f"{merged_dcd.name} ({frames} frames)"
        else:
            merged_dcd_str = f"{merged_dcd.name} (exists)"

    merged_log_str = "NOT FOUND"
    if merged_log.exists():
        try:
            with open(merged_log) as f:
                lines = sum(1 for line in f if line.strip())
                frames = max(0, lines - 1)
                merged_log_str = f"{merged_log.name} ({frames} frames)"
        except Exception:
            merged_log_str = f"{merged_log.name} (exists)"

    print(f"  Trajectory File  : {merged_dcd_str}")
    print(f"  Log File         : {merged_log_str}")

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

    print(
        f"  Resume-ready     : {'YES' if can_resume else 'NO'} (needs prod.chk, system.xml, solvated.pdb)"
    )


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
    rc = subprocess.call(cmd, cwd=cwd)
    if rc != 0:
        raise SystemExit(rc)


def _sync_tree(src_dir: Path, dst_dir: Path):
    if not src_dir.exists():
        print(f"[WARN] Sync source does not exist: {src_dir}")
        return
    if src_dir.resolve() == dst_dir.resolve():
        print(f"[INFO] Sync skipped (same path): {src_dir}")
        return
    dst_dir.mkdir(parents=True, exist_ok=True)
    copied = 0
    for p in src_dir.iterdir():
        if p.is_file():
            shutil.copy2(p, dst_dir / p.name)
            copied += 1
    print(f"[INFO] Synced {copied} files: {src_dir} -> {dst_dir}")


def openmm_prep_from_pdbid(
    pdbid: str, root_dir: str | None = None, ph: float = 7.4, sync_dir: str | None = None
):
    from colabmda.openmm.prepare import run_clean_by_pdbid

    base = Path(root_dir or os.getcwd()).resolve()
    outdir = base / pdbid / "prep"
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"[INFO] Prep output will be written to: {outdir}")
    run_clean_by_pdbid(pdbid, outdir=str(outdir), ph=ph)
    if sync_dir:
        _sync_tree(outdir, Path(sync_dir).resolve())


def openmm_prep_from_file(
    pdb_file: str, outdir: str, pdbid: str = "4ldj", ph: float = 7.0, sync_dir: str | None = None
):
    from colabmda.openmm.prepare import run_clean_from_file

    outdir_path = Path(outdir).resolve()
    print(f"[INFO] Prep output will be written to: {outdir_path}")
    run_clean_from_file(pdb_file, outdir, pdbid=pdbid, ph=ph)
    if sync_dir:
        _sync_tree(outdir_path, Path(sync_dir).resolve())


def openmm_run_colab(
    workdir: str,
    pdbid: str,
    total_ns: float,
    traj_interval: float,
    equil_time: float,
    checkpoint_ps: float,
    sync_dir: str | None,
    equil_only: bool = False,
):
    from colabmda.openmm.run import run_colab_md

    run_colab_md(
        workdir=workdir,
        pdbid=pdbid,
        total_ns=total_ns,
        traj_interval=traj_interval,
        equil_time=equil_time,
        checkpoint_ps=checkpoint_ps,
        sync_dir=sync_dir,
        equil_only=equil_only,
    )


def openmm_merge(
    pdbid_dir: str,
    topology: str | None,
    out_traj: str,
    out_log: str,
    stride: int = 1,
    center: bool = False,
    wrap: bool = False,
    mda: bool = False,
    selection: str | None = None,
    ca_only: bool = False,
    protein_only: bool = False,
):
    from colabmda.openmm.merge import (
        merge_logs,
        merge_trajectories,
        merge_trajectories_mda,
    )

    # Automatically trigger MDAnalysis mode if any MDAnalysis specific flags are active
    mda_active = mda or ca_only or protein_only or (selection is not None)

    if mda_active:
        merge_trajectories_mda(
            pdbid=pdbid_dir,
            topology=topology,
            out_traj=out_traj,
            stride=stride,
            center=center,
            wrap=wrap,
            selection=selection,
            ca_only=ca_only,
            protein_only=protein_only,
        )
    else:
        merge_trajectories(
            pdbid=pdbid_dir,
            topology=topology,
            out_traj=out_traj,
            stride=stride,
            center=center,
            wrap=wrap,
        )
    merge_logs(pdbid_dir, out_log, stride=stride)


def openmm_analysis(
    pdbid_dir: str,
    topology: str | None,
    trajectory: str | None,
    interval_ps: float | None,
    outdir: str | None,
    reference: str | None = None,
):
    from colabmda.openmm.analysis import run_trajectory_analysis

    result_dir = run_trajectory_analysis(
        pdbid=pdbid_dir,
        topology=topology,
        trajectory=trajectory,
        interval=interval_ps,
        outdir=outdir,
        reference=reference,
    )

    # If outdir exists, also copy equilibration plots there for completeness
    if outdir:
        out_path = Path(outdir)
        out_path.mkdir(parents=True, exist_ok=True)
        qc_file = Path(pdbid_dir) / "equilibration_qc.png"
        if qc_file.exists() and qc_file.resolve() != (out_path / "equilibration_qc.png").resolve():
            shutil.copy2(qc_file, out_path / "equilibration_qc.png")
            print(f"[INFO] Copied equilibration QC plot to {out_path}")

    from colabmda.openmm.figures import plot_analysis

    plot_analysis(result_dir)


def openmm_em(
    workdir: str,
    pdbid: str,
    ligand: str | None = None,
    keep_mg: bool = False,
    amber_prmtop: str | None = None,
    amber_inpcrd: str | None = None,
    padding_nm: float | str = "auto",
    protein_ff: str = "ff19SB",
    water_model: str = "opc",
    small_molecule_ff: str = "gaff-2.11",
    max_iterations: int = 0,
):
    from colabmda.openmm.em import run_em

    success = run_em(
        workdir,
        pdbid,
        ligand=ligand,
        keep_mg=keep_mg,
        amber_prmtop=amber_prmtop,
        amber_inpcrd=amber_inpcrd,
        padding_nm=padding_nm,
        protein_ff=protein_ff,
        water_model=water_model,
        small_molecule_ff=small_molecule_ff,
        max_iterations=max_iterations,
    )
    if not success:
        raise SystemExit(1)


def openmm_nvt(
    workdir: str,
    pdbid: str,
    equil_time: float,
    seed: int | None = None,
    protocol: str = "varmdyn",
):
    from colabmda.openmm.nvt import run_nvt

    success = run_nvt(workdir, pdbid, equil_time_ps=equil_time, seed=seed, protocol=protocol)
    if not success:
        raise SystemExit(1)


def openmm_npt(
    workdir: str,
    pdbid: str,
    equil_time: float,
    seed: int | None = None,
    protocol: str = "varmdyn",
):
    from colabmda.openmm.npt import run_npt

    success = run_npt(workdir, pdbid, equil_time_ps=equil_time, seed=seed, protocol=protocol)
    if not success:
        raise SystemExit(1)


def openmm_qc(workdir: str, strict: bool = True, protocol: str | None = None):
    from colabmda.openmm.qc import analyze_logs

    passed = analyze_logs(workdir, strict=strict, protocol=protocol)
    from colabmda.openmm.figures import plot_equilibration

    plot_equilibration(workdir)
    if strict and not passed:
        raise SystemExit(
            "Equilibration QC failed; inspect equilibration_qc.json/PNG before production"
        )
    return passed


openmm_check_equil = openmm_qc


def openmm_md(
    workdir: str,
    pdbid: str,
    total_ns: float,
    traj_interval: float,
    checkpoint_ps: float,
    sync_dir: str | None = None,
    seed: int | None = None,
):
    from colabmda.openmm.md import run_md

    success = run_md(
        workdir=workdir,
        pdbid=pdbid,
        total_ns=total_ns,
        traj_interval_ps=traj_interval,
        checkpoint_ps=checkpoint_ps,
        sync_dir=sync_dir,
        seed=seed,
    )
    if not success:
        raise SystemExit(1)


def openmm_compare(*args, **kwargs):
    from colabmda.openmm.figures import openmm_compare as render

    return render(*args, **kwargs)


def openmm_view(*args, **kwargs):
    from colabmda.openmm.figures import openmm_view as render

    return render(*args, **kwargs)


def openmm_snapshots(*args, **kwargs):
    from colabmda.openmm.figures import openmm_snapshots as render

    return render(*args, **kwargs)
