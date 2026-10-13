"""CLI invocation and subparser tests."""

import subprocess
import sys


def test_colabmda_help():
    res = subprocess.run(
        [sys.executable, "-m", "colabmda.cli", "--help"], capture_output=True, text=True
    )
    assert res.returncode == 0
    assert "openmm" in res.stdout
    assert "modeller" in res.stdout
    # Docking should be parked on the docking branch
    assert "dock" not in res.stdout


def test_openmm_help():
    res = subprocess.run(
        [sys.executable, "-m", "colabmda.cli", "openmm", "--help"], capture_output=True, text=True
    )
    assert res.returncode == 0
    assert "prep" in res.stdout
    assert "stage" in res.stdout
    assert "em" in res.stdout
    assert "nvt" in res.stdout
    assert "npt" in res.stdout
    assert "qc" in res.stdout
    assert "md" in res.stdout


def test_modeller_help():
    res = subprocess.run(
        [sys.executable, "-m", "colabmda.cli", "modeller", "--help"], capture_output=True, text=True
    )
    assert res.returncode == 0
    assert "build" in res.stdout
    assert "mutate" in res.stdout


def test_qc_alias():
    res = subprocess.run(
        [sys.executable, "-m", "colabmda.cli", "openmm", "qc", "--help"],
        capture_output=True,
        text=True,
    )
    assert res.returncode == 0
    assert "--warn-only" in res.stdout

    res_alias = subprocess.run(
        [sys.executable, "-m", "colabmda.cli", "openmm", "check-equil", "--help"],
        capture_output=True,
        text=True,
    )
    assert res_alias.returncode == 0
    assert "--warn-only" in res_alias.stdout


def test_merge_analysis_status_name_options():
    for cmd in ["merge", "analysis", "status"]:
        res = subprocess.run(
            [sys.executable, "-m", "colabmda.cli", "openmm", cmd, "--help"],
            capture_output=True,
            text=True,
        )
        assert res.returncode == 0
        assert "--name" in res.stdout
        assert "--replica" in res.stdout
        assert "--project" in res.stdout


def test_simulation_workdir_resolution(tmp_path):
    from colabmda.cli import _simulation_workdir

    # 1. Default should point to data/simulation
    res = _simulation_workdir(str(tmp_path), "test_sys")
    assert str(res) == str(tmp_path / "data" / "simulation" / "test_sys")

    # 2. If data/sim exists (legacy), it should discover it
    legacy_dir = tmp_path / "data" / "sim" / "legacy_sys"
    legacy_dir.mkdir(parents=True)
    res_legacy = _simulation_workdir(str(tmp_path), "legacy_sys")
    assert str(res_legacy) == str(legacy_dir)


def test_stage_creates_eq_and_replica(tmp_path):
    # Setup dummy structure in data/structure
    str_dir = tmp_path / "data" / "structure" / "4ldj" / "wt" / "prep"
    str_dir.mkdir(parents=True)
    dummy_pdb = str_dir / "4ldj_wt_cleaned.pdb"
    dummy_pdb.write_text("ATOM      1  N   MET A   1       0.000   0.000   0.000  1.00  0.00\nEND\n")

    res = subprocess.run(
        [
            sys.executable,
            "-m",
            "colabmda.cli",
            "openmm",
            "stage",
            "--pdb-file",
            str(dummy_pdb),
            "--name",
            "4ldj_wt",
            "--replica",
            "r1",
            "--project",
            str(tmp_path),
        ],
        capture_output=True,
        text=True,
    )
    assert res.returncode == 0
    # eq hub should always be created and populated
    assert (tmp_path / "data" / "simulation" / "4ldj_wt" / "eq" / "4ldj_wt_cleaned.pdb").exists()
    # r1 replica should be created and populated
    assert (tmp_path / "data" / "simulation" / "4ldj_wt" / "r1" / "4ldj_wt_cleaned.pdb").exists()

