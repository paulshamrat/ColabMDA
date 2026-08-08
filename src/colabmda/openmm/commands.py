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
    # Bundled workflow scripts (pdb-id download)
    "clean_by_pdbid": ("colabmda.openmm.engine", "pdbfixer_cleaning.py"),
    "merge": ("colabmda.openmm.engine", "openmm_trajmerge.py"),
    "analysis": ("colabmda.openmm.engine", "openmm_trajanalysis.py"),
    # Bundled colab-safe workflow scripts (local pdb file cleaning + robust resume)
    "clean_from_file": ("colabmda.openmm.engine", "pdbfixer_clean_fromfile.py"),
    "run_colab": ("colabmda.openmm.engine", "openmm_proteinwater_colab.py"),
    # New Modular Workflow
    "em": ("colabmda.openmm.modular", "em.py"),
    "nvt": ("colabmda.openmm.modular", "nvt.py"),
    "npt": ("colabmda.openmm.modular", "npt.py"),
    "check_equil": ("colabmda.openmm.modular", "check_equil.py"),
    "md": ("colabmda.openmm.modular", "md.py"),
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


def openmm_prep_from_pdbid(pdbid: str, root_dir: str | None = None, sync_dir: str | None = None):
    from colabmda.openmm.engine.pdbfixer_cleaning import run_clean_by_pdbid

    base = Path(root_dir or os.getcwd()).resolve()
    outdir = base / pdbid / "prep"
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"[INFO] Prep output will be written to: {outdir}")
    run_clean_by_pdbid(pdbid, outdir=str(outdir))
    if sync_dir:
        _sync_tree(outdir, Path(sync_dir).resolve())


def openmm_prep_from_file(
    pdb_file: str, outdir: str, pdbid: str = "4ldj", ph: float = 7.0, sync_dir: str | None = None
):
    from colabmda.openmm.engine.pdbfixer_clean_fromfile import run_clean_from_file

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
    from colabmda.openmm.engine.openmm_proteinwater_colab import run_colab_md

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
    from colabmda.openmm.engine.openmm_trajmerge import (
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
):
    from colabmda.openmm.engine.openmm_trajanalysis import run_trajectory_analysis

    run_trajectory_analysis(
        pdbid=pdbid_dir,
        topology=topology,
        trajectory=trajectory,
        interval=interval_ps,
        outdir=outdir,
    )

    # If outdir exists, also copy equilibration plots there for completeness
    if outdir:
        out_path = Path(outdir)
        out_path.mkdir(parents=True, exist_ok=True)
        qc_file = Path(pdbid_dir) / "equilibration_qc.png"
        if qc_file.exists() and qc_file.resolve() != (out_path / "equilibration_qc.png").resolve():
            shutil.copy2(qc_file, out_path / "equilibration_qc.png")
            print(f"[INFO] Copied equilibration QC plot to {out_path}")


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
):
    from colabmda.openmm.modular.em import run_em

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
    from colabmda.openmm.modular.nvt import run_nvt

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
    from colabmda.openmm.modular.npt import run_npt

    success = run_npt(workdir, pdbid, equil_time_ps=equil_time, seed=seed, protocol=protocol)
    if not success:
        raise SystemExit(1)


def openmm_check_equil(workdir: str, strict: bool = True):
    from colabmda.openmm.modular.check_equil import analyze_logs

    passed = analyze_logs(workdir)
    if strict and not passed:
        raise SystemExit(
            "Equilibration QC failed; inspect equilibration_qc.json/PNG before production"
        )
    return passed


def openmm_md(
    workdir: str,
    pdbid: str,
    total_ns: float,
    traj_interval: float,
    checkpoint_ps: float,
    sync_dir: str | None = None,
    seed: int | None = None,
):
    from colabmda.openmm.modular.md import run_md

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


def openmm_compare(series_list, outdir):
    import matplotlib.pyplot as plt
    import pandas as pd

    outpath = Path(outdir)
    outpath.mkdir(parents=True, exist_ok=True)

    try:
        plt.style.use("seaborn-v0_8-muted")
    except Exception:
        try:
            plt.style.use("seaborn-muted")
        except Exception:
            plt.style.use("ggplot")
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


def openmm_snapshots(
    config_path: str | None = None,
    output: str = "data/analysis/snapshots/master_3x8_snapshots.png",
    temp_dir: str = "scratch/master_grid_3x8",
    font_path: str | None = None,
):
    import json
    from pathlib import Path

    # 1. Load configuration
    default_config = {
        "camera_view": [
            -0.597478449,
            -0.617595792,
            -0.511463583,
            -0.502994895,
            0.785388291,
            -0.360776752,
            0.624511778,
            0.041705273,
            -0.779902637,
            0.000050262,
            -0.000109358,
            -182.850265503,
            31.930021286,
            32.076828003,
            34.52369679,
            141.553222656,
            224.141143799,
            -20.000000000,
        ],
        "align_ref_pdb": "analysis/clean_trajectories/g12c_protein.pdb",
        "stable_core_sel": "name CA and (resi 1-24 or resi 41-54 or resi 81-166)",
        "systems": {
            "WT": {
                "pdb": "analysis/clean_trajectories/wt_protein.pdb",
                "dcd": "analysis/clean_trajectories/wt_protein.dcd",
                "states": [1, 51, 91, 100, 101, 106, 111, 121],
                "times": [
                    "t = 0.00 ns",
                    "t = 0.50 ns",
                    "t = 0.90 ns",
                    "t = 0.99 ns",
                    "t = 1.00 ns",
                    "t = 1.05 ns",
                    "t = 1.10 ns",
                    "t = 1.20 ns",
                ],
                "mut_residue": 12,
            },
            "G12C": {
                "pdb": "analysis/clean_trajectories/g12c_protein.pdb",
                "dcd": "analysis/clean_trajectories/g12c_protein.dcd",
                "states": [1, 51, 91, 100, 101, 106, 111, 121],
                "times": [
                    "t = 0.00 ns",
                    "t = 0.50 ns",
                    "t = 0.90 ns",
                    "t = 0.99 ns",
                    "t = 1.00 ns",
                    "t = 1.05 ns",
                    "t = 1.10 ns",
                    "t = 1.20 ns",
                ],
                "mut_residue": 12,
            },
            "G12D": {
                "pdb": "analysis/clean_trajectories/g12d_protein.pdb",
                "dcd": "analysis/clean_trajectories/g12d_protein.dcd",
                "states": [1, 51, 91, 100, 101, 106, 111, 121],
                "times": [
                    "t = 0.00 ns",
                    "t = 0.50 ns",
                    "t = 0.90 ns",
                    "t = 0.99 ns",
                    "t = 1.00 ns",
                    "t = 1.05 ns",
                    "t = 1.10 ns",
                    "t = 1.20 ns",
                ],
                "mut_residue": 12,
            },
        },
    }

    if config_path:
        print(f"Loading custom configuration from {config_path}...")
        with open(config_path) as f:
            config = json.load(f)
    else:
        print("Using default KRAS 3x8 snapshot configuration...")
        config = default_config

    camera_view = config["camera_view"]
    align_ref_pdb = config.get("align_ref_pdb")
    stable_core_sel = config.get(
        "stable_core_sel", "name CA and (resi 1-24 or resi 41-54 or resi 81-166)"
    )
    systems = config["systems"]

    # 2. Check for PyMOL
    try:
        import pymol
        from pymol import cmd
    except ImportError:
        print("[ERROR] PyMOL python package is not installed in the active environment.")
        print("Please run this command using an environment with PyMOL (e.g. 'pymol-viz').")
        sys.exit(1)

    # Start PyMOL headlessly
    pymol.finish_launching(["pymol", "-qc"])

    # Create temp directory
    temp_path = Path(temp_dir)
    temp_path.mkdir(parents=True, exist_ok=True)

    rendered_files = {}

    for sys_name, sys_cfg in systems.items():
        pdb_path = sys_cfg["pdb"]
        dcd_path = sys_cfg["dcd"]
        states = sys_cfg["states"]
        mut_residue = sys_cfg.get("mut_residue", 12)

        if not os.path.exists(pdb_path):
            sys.exit(f"Error: PDB file not found: {pdb_path}")
        if not os.path.exists(dcd_path):
            sys.exit(f"Error: DCD file not found: {dcd_path}")

        print(f"Rendering snapshots for {sys_name}...")
        cmd.reinitialize()
        cmd.bg_color("white")

        # Load reference PDB for alignment if specified
        if align_ref_pdb and os.path.exists(align_ref_pdb):
            cmd.load(align_ref_pdb, "ref_align")

        # Load system
        cmd.load(pdb_path, "prot")
        cmd.load_traj(dcd_path, "prot")

        # Hide nonbonded
        cmd.hide("everything", "all")
        cmd.hide("nonbonded", "all")
        cmd.hide("nb_spheres", "all")

        # Step 1: Align all states of prot to state 1 of prot using stable core
        cmd.intra_fit(f"prot and {stable_core_sel}", 1)

        # Step 2: Align state 1 of prot to ref_align using stable core and apply to all states
        if align_ref_pdb and os.path.exists(align_ref_pdb):
            cmd.align(
                f"prot and {stable_core_sel}",
                f"ref_align and {stable_core_sel}",
                mobile_state=1,
                target_state=1,
            )
            cmd.delete("ref_align")

        cmd.dss("prot")
        cmd.cartoon("automatic", "prot")
        cmd.cartoon("loop", "prot and resid 57-75")
        cmd.show("cartoon", "prot")

        # Colors
        cmd.color("gray90", "prot")
        cmd.color("orange", "prot and resid 57-75")  # Switch II
        cmd.color("yellow", f"prot and resid {mut_residue}")

        # Mutation sticks
        cmd.show("sticks", f"prot and resi {mut_residue} and not name N+C+O+H")
        cmd.set("stick_radius", 0.25)

        # Transparency
        cmd.set("cartoon_transparency", 0.6, f"prot and not (resid 57-75 or resid {mut_residue})")
        cmd.set("cartoon_transparency", 0.0, f"resid 57-75 or resid {mut_residue}")

        # Ray tracing options
        cmd.set("ray_shadows", 1)
        cmd.set("ray_trace_mode", 1)
        cmd.set("ray_trace_gain", 0.5)
        cmd.set("depth_cue", 1)
        cmd.set("specular", 0.3)
        cmd.set("ray_opaque_background", 1)

        cmd.set_view(camera_view)

        rendered_files[sys_name] = []
        for state in states:
            p_out = temp_path / f"{sys_name}_state_{state}.png"
            rendered_files[sys_name].append(str(p_out))
            cmd.set("state", state)
            cmd.png(str(p_out), width=800, height=800, ray=1)

    pymol.cmd.quit()

    print("All snapshots rendered! Compiling image grid...")

    import PIL.ImageOps
    from PIL import Image, ImageDraw, ImageFont

    all_images = []
    for sys_name in systems.keys():
        for path in rendered_files[sys_name]:
            all_images.append(Image.open(path))

    # Bounding boxes for cropping
    bboxes = []
    for img in all_images:
        gray = img.convert("L")
        inverted = PIL.ImageOps.invert(gray)
        bbox = inverted.getbbox()
        if bbox:
            bboxes.append(bbox)

    orig_width, orig_height = all_images[0].size
    if bboxes:
        left_crop = min(b[0] for b in bboxes)
        top_crop = min(b[1] for b in bboxes)
        right_crop = max(b[2] for b in bboxes)
        bottom_crop = max(b[3] for b in bboxes)

        padding = 20
        left_crop = max(0, left_crop - padding)
        top_crop = max(0, top_crop - padding)
        right_crop = min(orig_width, right_crop + padding)
        bottom_crop = min(orig_height, bottom_crop + padding)

        cropped_images = [
            img.crop((left_crop, top_crop, right_crop, bottom_crop)) for img in all_images
        ]
    else:
        cropped_images = all_images

    # Grid configuration (N columns, M rows total)
    num_rows = len(systems)
    num_cols = len(list(systems.values())[0]["states"])

    panel_w, panel_h = cropped_images[0].size
    spacing_x = 12
    spacing_y = 15
    label_y_area = 120

    total_width = panel_w * num_cols + spacing_x * (num_cols - 1)
    total_height = (panel_h + label_y_area) * num_rows + spacing_y * (num_rows - 1)

    combined = Image.new("RGB", (total_width, total_height), "white")
    draw = ImageDraw.Draw(combined)

    # Font resolution
    font_path_to_use = font_path
    if not font_path_to_use:
        for path in [
            "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
            "/usr/share/fonts/truetype/freefont/FreeSansBold.ttf",
            "/usr/share/fonts/liberation/LiberationSans-Bold.ttf",
        ]:
            if os.path.exists(path):
                font_path_to_use = path
                break

    try:
        if font_path_to_use:
            font_panel_letter = ImageFont.truetype(font_path_to_use, 64)
            # Find corresponding regular font for time labels
            regular_font_path = font_path_to_use.replace("-Bold", "").replace("Bold", "")
            if not os.path.exists(regular_font_path):
                regular_font_path = font_path_to_use
            font_time = ImageFont.truetype(regular_font_path, 72)
        else:
            font_panel_letter = ImageFont.load_default()
            font_time = ImageFont.load_default()
    except OSError:
        font_panel_letter = ImageFont.load_default()
        font_time = ImageFont.load_default()

    img_idx = 0
    for r, (sys_name, sys_cfg) in enumerate(systems.items()):
        row_y = r * (panel_h + label_y_area + spacing_y)
        for c in range(num_cols):
            img = cropped_images[img_idx]
            x = c * (panel_w + spacing_x)
            combined.paste(img, (x, row_y))

            # Panel letter on the first panel of each row
            if c == 0:
                letter = chr(65 + r)  # A, B, C...
                draw.text((x + 20, row_y + 20), letter, fill="black", font=font_panel_letter)

            # Time label
            time_str = sys_cfg["times"][c]
            try:
                time_w = draw.textlength(time_str, font=font_time)
            except AttributeError:
                time_w, _ = draw.textsize(time_str, font=font_time)

            time_x = x + (panel_w - time_w) // 2
            time_y = row_y + panel_h + 20
            draw.text((time_x, time_y), time_str, fill="black", font=font_time)

            img_idx += 1

    # Draw faint dashed separator line only for G12C row (row r = 1) between columns 3 (0.99ns) and 4 (1.00ns)
    # If the system name G12C exists as the second row (standard KRAS layout)
    sys_list = list(systems.keys())
    if "G12C" in sys_list and sys_list.index("G12C") == 1 and num_cols >= 8:
        x_divider = 4 * panel_w + 3 * spacing_x + spacing_x // 2
        y_start = 1 * (panel_h + label_y_area + spacing_y)
        y_end = y_start + panel_h
        dash_len = 15
        gap_len = 10
        curr_y = y_start
        while curr_y < y_end:
            draw.line(
                [(x_divider, curr_y), (x_divider, min(curr_y + dash_len, y_end))],
                fill=(180, 180, 180),
                width=4,
            )
            curr_y += dash_len + gap_len

    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    combined.save(str(output_path), dpi=(300, 300))
    print(f"Master snapshots grid successfully saved to {output_path}!")

    # Clean up temp files
    for sys_name in systems.keys():
        for path in rendered_files[sys_name]:
            if os.path.exists(path):
                os.remove(path)
    if temp_path.exists():
        try:
            temp_path.rmdir()
        except OSError:
            pass
