"""Compatibility adapter for the combined protein-in-water command."""

from pathlib import Path

from colabmda.openmm.em import run_em
from colabmda.openmm.md import run_md
from colabmda.openmm.npt import run_npt
from colabmda.openmm.nvt import run_nvt


def run_colab_md(
    workdir,
    pdbid,
    total_ns=1.0,
    traj_interval=100.0,
    equil_time=100.0,
    checkpoint_ps=1000.0,
    sync_dir=None,
    equil_only=False,
):
    path = Path(workdir).resolve()
    path.mkdir(parents=True, exist_ok=True)
    if not (path / "prod.chk").exists() and not (path / "npt.state.xml").exists():
        # Retain the combined command's explicit ff14SB/TIP3P physical model.
        if not run_em(str(path), pdbid, protein_ff="ff14SB", water_model="tip3p", padding_nm=1.0):
            raise RuntimeError("Minimization failed")
        if not run_nvt(str(path), pdbid, equil_time_ps=equil_time, protocol="quick"):
            raise RuntimeError("NVT failed")
        if not run_npt(str(path), pdbid, equil_time_ps=equil_time, protocol="quick"):
            raise RuntimeError("NPT failed")
    if equil_only:
        return True
    return run_md(str(path), pdbid, total_ns, traj_interval, checkpoint_ps, sync_dir)
