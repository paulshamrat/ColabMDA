#!/usr/bin/env python3
import os
import re
import shutil
import subprocess
import sys
from importlib import resources
from pathlib import Path

SCRIPTS = {
    # Bundled workflow scripts (pdb-id download)
    "clean_by_pdbid": ("colabmda.openmm_pw.engine", "pdbfixer_cleaning.py"),
    "run_basic": ("colabmda.openmm_pw.engine", "openmm_proteinwater.py"),
    "merge": ("colabmda.openmm_pw.engine", "openmm_trajmerge.py"),
    "analysis": ("colabmda.openmm_pw.engine", "openmm_trajanalysis.py"),
    # Bundled colab-safe workflow scripts (local pdb file cleaning + robust resume)
    "clean_from_file": ("colabmda.openmm_pw.engine", "pdbfixer_clean_fromfile.py"),
    "run_colab": ("colabmda.openmm_pw.engine", "openmm_proteinwater_colab.py"),
    # New Modular Workflow
    "em": ("colabmda.openmm_pw.modular", "01_em.py"),
    "nvt": ("colabmda.openmm_pw.modular", "02_nvt.py"),
    "npt": ("colabmda.openmm_pw.modular", "03_npt.py"),
    "check_equil": ("colabmda.openmm_pw.modular", "check_equil.py"),
    "md": ("colabmda.openmm_pw.modular", "04_md.py"),
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


def openmm_prep_from_pdbid(pdbid: str, root_dir: str | None = None, sync_dir: str | None = None):
    pkg, name = SCRIPTS["clean_by_pdbid"]
    base = Path(root_dir or os.getcwd()).resolve()
    outdir = base / pdbid / "prep"
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"[INFO] Prep output will be written to: {outdir}")
    _run(_script_path(pkg, name), [pdbid, "--outdir", str(outdir)])
    if sync_dir:
        _sync_tree(outdir, Path(sync_dir).resolve())


def openmm_prep_from_file(
    pdb_file: str, outdir: str, pdbid: str = "4ldj", ph: float = 7.0, sync_dir: str | None = None
):
    pkg, name = SCRIPTS["clean_from_file"]
    outdir_path = Path(outdir).resolve()
    print(f"[INFO] Prep output will be written to: {outdir_path}")
    _run(
        _script_path(pkg, name),
        ["--in", pdb_file, "--outdir", outdir, "--pdbid", pdbid, "--ph", str(ph)],
    )
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
):
    argv = [
        workdir,
        "--pdbid",
        pdbid,
        "--total-ns",
        str(total_ns),
        "--traj-interval",
        str(traj_interval),
        "--equil-time",
        str(equil_time),
        "--checkpoint-ps",
        str(checkpoint_ps),
    ]
    if sync_dir:
        argv += ["--sync-dir", sync_dir]
    pkg, name = SCRIPTS["run_colab"]
    _run(_script_path(pkg, name), argv)


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
    argv = [pdbid_dir, "--out-traj", out_traj, "--out-log", out_log, "--stride", str(stride)]
    if topology:
        argv += ["--topology", topology]
    if center:
        argv += ["--center"]
    if wrap:
        argv += ["--wrap"]
    if mda:
        argv += ["--mda"]
    if selection:
        argv += ["--selection", selection]
    if ca_only:
        argv += ["--ca-only"]
    if protein_only:
        argv += ["--protein-only"]
    pkg, name = SCRIPTS["merge"]
    _run(_script_path(pkg, name), argv)


def openmm_analysis(
    pdbid_dir: str,
    topology: str | None,
    trajectory: str | None,
    interval_ps: float | None,
    outdir: str | None,
):
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

    # If outdir exists, also copy equilibration plots there for completeness
    if outdir:
        out_path = Path(outdir)
        out_path.mkdir(parents=True, exist_ok=True)
        qc_file = Path(pdbid_dir) / "equilibration_qc.png"
        if qc_file.exists():
            shutil.copy2(qc_file, out_path / "equilibration_qc.png")
            print(f"[INFO] Copied equilibration QC plot to {out_path}")


def openmm_em(workdir: str, pdbid: str):
    pkg, name = SCRIPTS["em"]
    _run(_script_path(pkg, name), [workdir, "--pdbid", pdbid])


def openmm_nvt(workdir: str, pdbid: str, equil_time: float, seed: int | None = None):
    pkg, name = SCRIPTS["nvt"]
    argv = [workdir, "--pdbid", pdbid, "--equil-time", str(equil_time)]
    if seed is not None:
        argv += ["--seed", str(seed)]
    _run(_script_path(pkg, name), argv)


def openmm_npt(workdir: str, pdbid: str, equil_time: float):
    pkg, name = SCRIPTS["npt"]
    _run(_script_path(pkg, name), [workdir, "--pdbid", pdbid, "--equil-time", str(equil_time)])


def openmm_check_equil(workdir: str):
    pkg, name = SCRIPTS["check_equil"]
    _run(_script_path(pkg, name), [workdir])


def openmm_md(
    workdir: str,
    pdbid: str,
    total_ns: float,
    traj_interval: float,
    checkpoint_ps: float,
    sync_dir: str | None = None,
):
    pkg, name = SCRIPTS["md"]
    argv = [
        workdir,
        "--pdbid",
        pdbid,
        "--total-ns",
        str(total_ns),
        "--traj-interval",
        str(traj_interval),
        "--checkpoint-ps",
        str(checkpoint_ps),
    ]
    if sync_dir:
        argv += ["--sync-dir", sync_dir]
    _run(_script_path(pkg, name), argv)


def openmm_compare(series_list, outdir):
    import matplotlib.pyplot as plt
    import pandas as pd

    outpath = Path(outdir)
    outpath.mkdir(parents=True, exist_ok=True)

    plt.style.use("seaborn-v0_8-muted")
    plt.rcParams.update({"font.size": 12, "axes.grid": True, "grid.alpha": 0.3})

    metrics = {
        "rmsd.csv": ("Time (ps)", "RMSD (Å)", "System Stability (RMSD)"),
        "rg.csv": ("Time (ps)", "Radius of Gyration (Å)", "Compactness (Rg)"),
        "rmsf.csv": ("Residue Index", "RMSF (Å)", "Flexibility (RMSF)"),
    }

    def aggregate_system(dirs, metric_file):
        dfs = []
        for d in dirs:
            p = Path(d) / metric_file
            if p.exists():
                dfs.append(pd.read_csv(p))
        if not dfs:
            return None, None, None
        combined = pd.concat(dfs)
        col_name = dfs[0].columns[1]
        grouped = combined.groupby(combined.iloc[:, 0])
        mean = grouped[col_name].mean()
        std = grouped[col_name].std()
        return mean.index, mean.values, std.values

    for filename, (xlabel, ylabel, title) in metrics.items():
        plt.figure(figsize=(10, 6))
        found_any = False
        for item in series_list:
            if "=" not in item:
                continue
            label, dirs_str = item.split("=", 1)
            dirs = [d.strip() for d in dirs_str.split(",")]
            x, mean, std = aggregate_system(dirs, filename)
            if x is None:
                continue
            found_any = True
            p = plt.plot(x, mean, label=f"{label} (avg)", lw=2)
            color = p[0].get_color()
            plt.fill_between(x, mean - std, mean + std, color=color, alpha=0.2)

        if found_any:
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.title(title)
            plt.legend()
            plt.tight_layout()
            plt.savefig(outpath / filename.replace(".csv", "_avg.png"), dpi=300)
        plt.close()
    print(f"✅ Aggregate plots saved in {outdir}")


def openmm_view(pdb_dir: str | None, topology: str | None, trajectory: str | None, resi: int = 12):
    import subprocess
    from pathlib import Path

    # 1. Resolve simulation directory
    sim_dir = Path(pdb_dir or os.getcwd()).resolve()
    if not sim_dir.is_dir():
        sys.exit(f"Error: directory '{sim_dir}' does not exist.")

    # 2. Resolve topology and trajectory filenames
    topo_name = topology
    traj_name = trajectory

    if not topo_name:
        for candidate in ["prod_full.pdb", "kras_protein.pdb", "solvated.pdb"]:
            if (sim_dir / candidate).is_file():
                topo_name = candidate
                break
        if not topo_name:
            sys.exit(f"Error: No topology PDB file found in '{sim_dir}'. Run merge first.")

    if not traj_name:
        for candidate in ["prod_full.dcd", "kras_protein.dcd"]:
            if (sim_dir / candidate).is_file():
                traj_name = candidate
                break
        if not traj_name:
            sys.exit(f"Error: No trajectory DCD file found in '{sim_dir}'. Run merge first.")

    topo_path = sim_dir / topo_name
    traj_path = sim_dir / traj_name

    if not topo_path.is_file():
        sys.exit(f"Error: topology file not found: {topo_path}")
    if not traj_path.is_file():
        sys.exit(f"Error: trajectory file not found: {traj_path}")

    # Find crystal reference PDB (e.g., 4ldj_wt.pdb) to align coordinate space
    ref_name = f"{sim_dir.parent.name}.pdb"
    ref_path = sim_dir / ref_name
    if not ref_path.is_file():
        for p in sim_dir.glob("*.pdb"):
            if (
                p.name not in [topo_name, "solvated.pdb", "equilibrated.pdb"]
                and "cleaned" not in p.name
                and "protein" not in p.name
            ):
                ref_name = p.name
                ref_path = p
                break

    align_cmds = ""
    if ref_path.is_file():
        align_cmds = f"""
# 2b. Align trajectory to the crystal reference to ensure identical viewing orientation
load {ref_name}, crystal_ref
align kras and name CA, crystal_ref and name CA, mobile_state=1, target_state=1
delete crystal_ref
"""

    # 3. Create visualize.pml file
    pml_content = f"""# PyMOL Script generated by ColabMDA
# 1. Clear out all the old overlapping objects from memory
reinitialize
bg_color white

# 2. Load the trajectory files
load {topo_name}, kras
load_traj {traj_name}, kras
{align_cmds}
# 3. Hide solvent water and ions
hide everything, all
hide nonbonded, all
hide nb_spheres, all

# 4. Freeze the backbone tumbling rotation frame-by-frame
intra_fit kras and name CA

# 5. Generate your crisp secondary structure ribbon
dss kras
cartoon automatic, kras
show cartoon, kras

# 6. Apply your smooth cyan color scheme
color cyan, kras and name C*
util.cnc("kras")

# 7. Highlight the mutated residue side chain sticks
show sticks, kras and resi {resi} and not name N+C+O+H
color yellow, kras and resi {resi} and name SG
set stick_radius, 0.25

# 8. Focus camera right onto the protein
zoom kras and polymer, buffer=4
"""
    pml_file = sim_dir / "visualize.pml"
    try:
        pml_file.write_text(pml_content)
        print(f"→ Generated PyMOL script: {pml_file}")
    except Exception as e:
        sys.exit(f"Error writing PyMOL script: {e}")

    # 4. Run PyMOL
    if not shutil.which("pymol"):
        print("[WARN] PyMOL executable ('pymol') not found in PATH.")
        print("Please install PyMOL or run manually using the generated visualize.pml file:")
        print(f"  cd {sim_dir} && pymol visualize.pml")
        return

    print(f"→ Launching PyMOL: pymol visualize.pml in {sim_dir}...")
    try:
        subprocess.run(["pymol", "visualize.pml"], cwd=str(sim_dir), check=True)
    except Exception as e:
        print(f"Error running PyMOL: {e}")
