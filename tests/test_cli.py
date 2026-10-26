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
    dummy_pdb.write_text(
        "ATOM      1  N   MET A   1       0.000   0.000   0.000  1.00  0.00\nEND\n"
    )

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


def test_top_level_single_keyword_help():
    res = subprocess.run(
        [sys.executable, "-m", "colabmda.cli", "--help"],
        capture_output=True,
        text=True,
    )
    assert res.returncode == 0
    for cmd in ["prep", "stage", "sim", "merge", "analyze", "compare", "colab"]:
        assert cmd in res.stdout


def test_top_level_stage_command(tmp_path, monkeypatch):
    monkeypatch.setenv("COLABMDA_RUN_ID", "cdkl5_staged_v1")
    res = subprocess.run(
        [
            sys.executable,
            "-m",
            "colabmda.cli",
            "stage",
            "--state",
            "atp_2mg_eb2",
            "--variant",
            "WT",
            "--replicas",
            "1,2,3",
            "--root",
            str(tmp_path),
        ],
        capture_output=True,
        text=True,
    )
    assert res.returncode == 0
    variant_base = tmp_path / "data/atp_2mg_eb2/simulation/cdkl5_staged_v1/WT"
    for folder in ["01_prep", "02_eq", "03_md/r1", "03_md/r2", "03_md/r3"]:
        assert (variant_base / folder).is_dir()
    assert (tmp_path / "data/atp_2mg_eb2/analysis/cdkl5_staged_v1/WT/r1").is_dir()
    assert not (variant_base / "04_analysis").exists()
    assert (tmp_path / "logs/atp_2mg_eb2/cdkl5_staged_v1/WT").is_dir()
    assert not (tmp_path / "scripts").exists()


def test_top_level_subcommands_help():
    for subcmd in ["prep", "stage", "sim", "merge", "analyze", "compare", "colab"]:
        res = subprocess.run(
            [sys.executable, "-m", "colabmda.cli", subcmd, "--help"],
            capture_output=True,
            text=True,
        )
        assert res.returncode == 0, f"Command {subcmd} --help failed: {res.stderr}"


def test_top_level_compare_command(tmp_path, monkeypatch):
    monkeypatch.setenv("COLABMDA_RUN_ID", "run01")
    for var in ["WT", "L119R"]:
        rep_dir = tmp_path / f"data/apo/analysis/run01/{var}/r1"
        rep_dir.mkdir(parents=True, exist_ok=True)
        (rep_dir / "rmsd.csv").write_text("Time (ps),RMSD (Å)\n0,1.0\n10,1.2\n")
        (rep_dir / "rg.csv").write_text("Time (ps),Radius of Gyration (Å)\n0,15.0\n10,15.1\n")
        (rep_dir / "rmsf.csv").write_text("Residue Index,RMSF (Å)\n1,0.5\n2,0.6\n")

    res = subprocess.run(
        [
            sys.executable,
            "-m",
            "colabmda.cli",
            "compare",
            "--state",
            "apo",
            "--variants",
            "WT,L119R",
            "--root",
            str(tmp_path),
        ],
        capture_output=True,
        text=True,
    )
    assert res.returncode == 0, f"Compare failed: {res.stderr}"
    compare_fig_dir = tmp_path / "data/apo/analysis/run01/compare/figures"
    assert compare_fig_dir.is_dir()
    for name in ["rmsd_avg.png", "rg_avg.png", "rmsf_avg.png"]:
        assert (compare_fig_dir / name).is_file()
