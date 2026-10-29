#!/usr/bin/env python3
"""
analysis.py

Analyze trajectory (RMSD, Rg, RMSF) for a single system or multiple replicas.
Automatically detects r1, r2, r3... subfolders if present.

Usage:
    colabmda openmm analysis --name 4ldj_wt
"""

import argparse
import glob
import hashlib
import json
import os
from pathlib import Path

import MDAnalysis as mda
import mdtraj as md
import numpy as np
import pandas as pd
from MDAnalysis.analysis import rms


def _data_lines(path):
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s[0] == "#" or s[0].isalpha():
                continue
            yield s


def _infer_interval(sim_dir, n_frames):
    merged_log = (
        os.path.join(sim_dir, "md.log")
        if os.path.isfile(os.path.join(sim_dir, "md.log"))
        else os.path.join(sim_dir, "prod_full.log")
    )
    if os.path.isfile(merged_log):
        times = []
        for s in _data_lines(merged_log):
            parts = s.split(",") if "," in s else s.split()
            if len(parts) >= 2:
                try:
                    times.append(float(parts[1]))
                except Exception:
                    continue
        if len(times) >= 2:
            return (times[-1] - times[0]) / max(1, n_frames - 1)
    return None


def _residue_labels(topology, ca_indices, reference=None):
    """Map C-alpha residues back to native author labels, with an identity gate."""
    residues = [topology.atom(i).residue for i in ca_indices]
    if reference is None:
        # Amber topology indices do not encode native author numbering.
        return [r.resSeq for r in residues]
    from Bio.PDB import MMCIFParser, PDBParser

    reference = Path(reference)
    parser = (
        MMCIFParser(QUIET=True) if reference.suffix.lower() == ".cif" else PDBParser(QUIET=True)
    )
    structure = parser.get_structure("reference", str(reference))
    native = [r for r in next(iter(structure)).get_residues() if r.id[0] == " " and "CA" in r]
    aliases = {
        "HID": "HIS",
        "HIE": "HIS",
        "HIP": "HIS",
        "ASH": "ASP",
        "GLH": "GLU",
        "CYX": "CYS",
        "LYN": "LYS",
    }

    def normalize(name):
        return aliases.get(name, name)

    if len(native) != len(residues) or any(
        normalize(a.resname) != normalize(b.name) for a, b in zip(native, residues)
    ):
        raise ValueError("Reference/trajectory C-alpha residue identity mismatch")
    labels = [r.id[1] if not r.id[2].strip() else f"{r.id[1]}{r.id[2].strip()}" for r in native]
    if len(set(labels)) != len(labels):
        labels = [f"{r.parent.id}:{label}" for r, label in zip(native, labels)]
    return labels


def analyze_single(sim_dir, topo_path, traj_path, interval_user=None, reference=None):
    u = mda.Universe(topo_path, traj_path)
    n_frames = len(u.trajectory)
    interval = interval_user or _infer_interval(sim_dir, n_frames) or 10.0
    # Correct time: Frames are at 1*int, 2*int ... N*int
    times = (np.arange(n_frames) + 1) * interval

    # RMSD
    rmsd_calc = rms.RMSD(u, u, select="backbone", ref_frame=0)
    rmsd_calc.run()
    rmsd = rmsd_calc.results.rmsd[:, 2]

    # Rg
    heavy = u.select_atoms("protein and not name H*")
    rg = []
    for _ts in u.trajectory:
        coords = heavy.positions
        cog = coords.mean(axis=0)
        rg.append(np.sqrt(((coords - cog) ** 2).sum(axis=1).mean()))
    rg = np.array(rg)

    # RMSF
    t = md.load(traj_path, top=topo_path)
    ca_idx = t.topology.select("name CA")
    t.superpose(t, 0, atom_indices=ca_idx)
    rmsf = md.rmsf(t, t, atom_indices=ca_idx) * 10.0
    if reference is None and str(topo_path).endswith((".prmtop", ".parm7")):
        raise ValueError("Amber topology requires a reference structure for native residue labels")
    resids = _residue_labels(t.topology, ca_idx, reference)

    return {"times": times, "rmsd": rmsd, "rg": rg, "rmsf": rmsf, "resids": resids}


def run_trajectory_analysis(
    pdbid, topology=None, trajectory=None, interval=None, outdir=None, reference=None
):
    sim_dir = os.path.abspath(pdbid)
    outdir_path = Path(outdir or f"analysis_{os.path.basename(sim_dir)}").resolve()
    outdir_path.mkdir(parents=True, exist_ok=True)

    # Detect replicas
    replica_dirs = sorted(glob.glob(os.path.join(sim_dir, "r[0-9]*")))
    results = []
    labels = []
    res_outdirs = []

    if not replica_dirs:
        # Single mode
        topo = topology or os.path.join(sim_dir, "solvated.pdb")
        traj = trajectory or (
            os.path.join(sim_dir, "md.dcd")
            if os.path.exists(os.path.join(sim_dir, "md.dcd"))
            else os.path.join(sim_dir, "prod_full.dcd")
        )
        if os.path.exists(traj):
            results.append(analyze_single(sim_dir, topo, traj, interval, reference))
            labels.append("Sim")
            # For single mode, save directly in outdir
            res_outdirs = [outdir_path]
    else:
        # Multi-replica mode
        print(
            f"Detected {len(replica_dirs)} replicas: {[os.path.basename(r) for r in replica_dirs]}"
        )
        for rd in replica_dirs:
            topo = os.path.join(rd, "solvated.pdb")
            traj = (
                os.path.join(rd, "md.dcd")
                if os.path.exists(os.path.join(rd, "md.dcd"))
                else os.path.join(rd, "prod_full.dcd")
            )
            if os.path.exists(traj):
                results.append(analyze_single(rd, topo, traj, interval, reference))
                rep_name = os.path.basename(rd)
                labels.append(rep_name)
                # Individual replica output folder
                rep_out = outdir_path / rep_name
                rep_out.mkdir(parents=True, exist_ok=True)
                res_outdirs.append(rep_out)

                # Copy equilibration QC plot if it exists in the replica folder
                qc_src = Path(rd) / "equilibration_qc.png"
                if qc_src.exists():
                    import shutil

                    shutil.copy2(qc_src, rep_out / "equilibration_qc.png")

    if not results:
        raise ValueError("Error: No valid trajectory data found.")

    for result, target in zip(results, res_outdirs):
        pd.DataFrame({"Time (ps)": result["times"], "RMSD (Å)": result["rmsd"]}).to_csv(
            target / "rmsd.csv", index=False
        )
        pd.DataFrame({"Time (ps)": result["times"], "Rg (Å)": result["rg"]}).to_csv(
            target / "rg.csv", index=False
        )
        pd.DataFrame({"Residue": result["resids"], "RMSF (Å)": result["rmsf"]}).to_csv(
            target / "rmsf.csv", index=False
        )
    if len(results) > 1:
        agg = outdir_path / "aggregate"
        agg.mkdir(exist_ok=True)
        for key in ("rmsd", "rg", "rmsf"):
            length = min(len(r[key]) for r in results)
            if key == "rmsf" and any(r["resids"] != results[0]["resids"] for r in results):
                raise ValueError("Replica residue identities differ")
            values = np.array([r[key][:length] for r in results])
            x = results[0]["resids"] if key == "rmsf" else results[0]["times"][:length]
            pd.DataFrame(
                {
                    "Residue" if key == "rmsf" else "Time (ps)": x,
                    "Mean": values.mean(axis=0),
                    "Std": values.std(axis=0),
                }
            ).to_csv(agg / f"{key}_avg.csv", index=False)
    metadata = {"analysis_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    snapshot = Path(__file__).resolve().parents[3] / "SNAPSHOT.json"
    if snapshot.exists():
        metadata["tool"] = json.loads(snapshot.read_text())
    if reference:
        metadata["reference"] = str(Path(reference).resolve())
        metadata["reference_sha256"] = hashlib.sha256(Path(reference).read_bytes()).hexdigest()
    (outdir_path / "analysis.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"Analysis tables saved in: {outdir_path}")
    return outdir_path


def aggregate_reports(folders, outdir):
    """Aggregate saved replica metrics only when time/native residue axes match."""
    folders = [Path(p) for p in folders]
    if len(folders) < 2 or len(set(folders)) != len(folders):
        raise ValueError("Aggregate requires at least two distinct replicas")
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "AGGREGATE_DONE").unlink(missing_ok=True)
    tables = {}
    hashes = {}
    for metric in ["rmsd", "rg", "rmsf"]:
        inputs = [folder / f"{metric}.csv" for folder in folders]
        frames = [pd.read_csv(p) for p in inputs]
        axis = frames[0].iloc[:, 0]
        if axis.empty:
            raise ValueError("Cannot aggregate empty tables")
        values = []
        for path, frame in zip(inputs, frames):
            if len(frame) != len(axis) or frame.columns.tolist() != frames[0].columns.tolist():
                raise ValueError("Replica table lengths or columns differ")
            current = frame.iloc[:, 0]
            if metric == "rmsf":
                matched = current.astype(str).tolist() == axis.astype(str).tolist()
            else:
                matched = np.allclose(
                    current.to_numpy(float), axis.to_numpy(float), atol=1e-8, rtol=0
                )
            if not matched:
                raise ValueError("Replica time or native residue axes differ")
            data = frame.iloc[:, 1].to_numpy(float)
            if not np.isfinite(data).all():
                raise ValueError("Nonfinite replica metric")
            values.append(data)
            hashes[str(path.resolve())] = hashlib.sha256(path.read_bytes()).hexdigest()
        data = np.asarray(values)
        tables[metric] = pd.DataFrame(
            {axis.name: axis, "Mean": data.mean(axis=0), "Std": data.std(axis=0)}
        )
    for metric, frame in tables.items():
        frame.to_csv(outdir / f"{metric}_avg.csv", index=False)
    (outdir / "analysis.json").write_text(
        json.dumps(
            {
                "replicas": [str(p.resolve()) for p in folders],
                "input_sha256": hashes,
                "standard_deviation": "population (ddof=0)",
                "analysis_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            },
            indent=2,
        )
        + "\n"
    )
    return outdir


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("pdbid", help="Simulation directory")
    parser.add_argument("-t", "--topology")
    parser.add_argument("-x", "--trajectory")
    parser.add_argument("-i", "--interval", type=float)
    parser.add_argument("-o", "--outdir")
    parser.add_argument(
        "--reference", help="Native-numbered reference structure for Amber topology"
    )
    args = parser.parse_args()

    run_trajectory_analysis(
        args.pdbid,
        topology=args.topology,
        trajectory=args.trajectory,
        interval=args.interval,
        outdir=args.outdir,
        reference=args.reference,
    )


if __name__ == "__main__":
    main()
