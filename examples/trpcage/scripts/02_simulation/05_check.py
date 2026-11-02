"""Verify the complete example before certifying software completion."""

import hashlib
import json
import os
from pathlib import Path

import mdtraj as md
import numpy as np
import pandas as pd

root = Path(__file__).resolve().parents[2]
run_id = os.environ["COLABMDA_RUN_ID"]
if not run_id or any(
    c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-" for c in run_id
):
    raise ValueError("Invalid run ID")
run = root / "data/02_simulation" / run_id
system = run / "apo/WT"
prod = system / "03_md/r1"
analysis = system / "04_analysis"
# Keep this aligned with 03_production.sh when changing the example duration.
requests = run / "production.json"
ns = (
    float(json.loads(requests.read_text())[-1]["requested_ns"])
    if requests.exists()
    else float(os.environ.get("COLABMDA_EXAMPLE_NS", "1"))
)
expected_times = np.arange(10, ns * 1000 + 0.001, 10)
if ns <= 0 or not np.isclose(ns * 10, round(ns * 10)):
    raise ValueError("Use a positive duration divisible by 0.1 ns")
if len(list(prod.glob("prod_*to*ps.dcd"))) != round(ns * 10):
    raise RuntimeError("Unexpected production chunk count")
trajectory = md.load(str(prod / "prod_keep_ligand.dcd"), top=str(prod / "stripped.prmtop"))
ca = md.load(str(prod / "prod_ca.dcd"), top=str(prod / "ca.prmtop"))
if (
    trajectory.n_frames != len(expected_times)
    or ca.n_frames != len(expected_times)
    or ca.n_atoms != 20
):
    raise RuntimeError("Unexpected frame or C-alpha count")
if not np.isfinite(trajectory.xyz).all():
    raise RuntimeError("Nonfinite coordinates")
log = pd.read_csv(prod / "prod_full.log")
if len(log) != len(expected_times) or not np.allclose(log["Time (ps)"], expected_times):
    raise RuntimeError("Trajectory times have gaps, overlaps or wrong duration")
if not np.isfinite(log.select_dtypes(include="number").to_numpy()).all():
    raise RuntimeError("Nonfinite production log")
for metric in ["rmsd", "rg", "rmsf"]:
    values = pd.read_csv(analysis / f"{metric}.csv")
    if values.empty or not np.isfinite(values.to_numpy()).all():
        raise RuntimeError("Invalid analysis table")
if pd.read_csv(analysis / "rmsf.csv")["Residue"].tolist() != list(range(1, 21)):
    raise RuntimeError("Native author numbering was lost")
for name in ["rmsd_vs_time.png", "rg_vs_time.png", "rmsf_per_residue.png"]:
    path = analysis / "figures" / name
    if not path.is_file() or not path.stat().st_size:
        raise RuntimeError("Missing figure")
setup = json.loads((root / "data/00_setup/example.json").read_text())
qc = json.loads((system / "02_eq/equilibration_qc.json").read_text())
record = {
    "example": "trpcage",
    "pdb": "1L2Y",
    "residues": 20,
    "production_ns": ns,
    "frames": trajectory.n_frames,
    "solute_atoms": trajectory.n_atoms,
    "software_pass": True,
    "equilibration_qc_pass": qc["passed"],
    "setup": setup,
    "analysis": json.loads((analysis / "analysis.json").read_text()),
    "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
}
(run / "summary.json").write_text(json.dumps(record, indent=2) + "\n")
(run / "RUN_DONE").write_text(
    "Software checks passed; physical equilibration is reported separately in summary.json\n"
)
print("RUN_DONE", run, "physical QC:", qc["passed"])
