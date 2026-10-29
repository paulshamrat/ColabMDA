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
    assert not (tmp_path / "data/atp_2mg_eb2/analysis").exists()
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


def test_sim_positional_stage_argument(monkeypatch):
    import argparse
    from unittest.mock import MagicMock

    from colabmda.cli import _dispatch_sim

    def fake_openmm_stages(*args, **kwargs):
        return None

    monkeypatch.setenv("COLABMDA_RUN_ID", "run01")
    monkeypatch.setattr("colabmda.core.workflow.verify_prepared", lambda *a: True)
    monkeypatch.setattr("colabmda.openmm.commands.openmm_em", fake_openmm_stages)
    monkeypatch.setattr("colabmda.openmm.commands.openmm_nvt", fake_openmm_stages)
    monkeypatch.setattr("colabmda.openmm.commands.openmm_npt", fake_openmm_stages)
    monkeypatch.setattr("colabmda.openmm.commands.openmm_qc", fake_openmm_stages)
    monkeypatch.setattr("colabmda.openmm.commands.openmm_md", fake_openmm_stages)

    # 1. Positional 'eq' argument
    args_eq = argparse.Namespace(
        root=None,
        state="apo",
        variant="WT",
        stage_pos="eq",
        stage="em,nvt,npt,qc",
        engine="openmm",
        replica="r1",
        total_ns=100.0,
        traj_interval=10.0,
        checkpoint_ps=1000.0,
        steps=None,
        equil_time=100.0,
        protocol="staged",
        software_test=False,
    )
    # Ensure dispatch accepts positional eq without error
    # (Mock ProjectPaths to avoid filesystem dependency)
    mock_paths = MagicMock()
    mock_paths.prepared_prmtop.return_value.exists.return_value = True
    mock_paths.prepared_inpcrd.return_value.exists.return_value = True
    mock_paths.eq_dir.return_value = MagicMock()
    mock_paths.model_pdb.return_value.exists.return_value = False
    monkeypatch.setattr("colabmda.core.project.ProjectPaths", lambda *a, **k: mock_paths)

    code = _dispatch_sim(args_eq)
    assert code == 0


def test_structure_target_dir_resolution(tmp_path):
    from colabmda.core.project import ProjectPaths

    # 1. Fresh state directory with no target yet
    paths = ProjectPaths(tmp_path, state="apo")
    assert paths.structure_dir == tmp_path / "data/apo/structure"

    # 2. When target subfolder like O76039_4BGQ is populated
    target_folder = tmp_path / "data/apo/structure/O76039_4BGQ"
    (target_folder / "models/WT").mkdir(parents=True)
    (target_folder / "STRUCTURE_DONE").touch()

    paths_auto = ProjectPaths(tmp_path, state="apo")
    assert paths_auto.structure_dir == target_folder

    # 3. Explicit target override
    paths_explicit = ProjectPaths(tmp_path, state="apo", target="O76039_4BGQ")
    assert paths_explicit.structure_dir == target_folder


def test_push_pull_cli_parsers():
    # 1. Top-level cmda push --help
    res_push = subprocess.run(
        [sys.executable, "-m", "colabmda.cli", "push", "--help"],
        capture_output=True,
        text=True,
    )
    assert res_push.returncode == 0
    assert "--session" in res_push.stdout
    assert "--project" in res_push.stdout

    # 2. Top-level cmda pull --help
    res_pull = subprocess.run(
        [sys.executable, "-m", "colabmda.cli", "pull", "--help"],
        capture_output=True,
        text=True,
    )
    assert res_pull.returncode == 0
    assert "--session" in res_pull.stdout
    assert "--stage" in res_pull.stdout
    assert "--variant" in res_pull.stdout
    assert "--traj" in res_pull.stdout
    assert "--all" in res_pull.stdout

    # 3. cmda colab push and pull
    res_colab_push = subprocess.run(
        [sys.executable, "-m", "colabmda.cli", "colab", "push", "--help"],
        capture_output=True,
        text=True,
    )
    assert res_colab_push.returncode == 0
    assert "--session" in res_colab_push.stdout

    res_colab_pull = subprocess.run(
        [sys.executable, "-m", "colabmda.cli", "colab", "pull", "--help"],
        capture_output=True,
        text=True,
    )
    assert res_colab_pull.returncode == 0
    assert "--stage" in res_colab_pull.stdout
    assert "--variant" in res_colab_pull.stdout
    assert "--traj" in res_colab_pull.stdout
    assert "--all" in res_colab_pull.stdout


def test_colab_pull_dispatch_and_fast_path(monkeypatch, tmp_path):
    from colabmda.remote.colab import colab_pull

    calls = []

    def mock_subprocess_call(cmd):
        calls.append(cmd)
        if len(cmd) > 1 and cmd[0] == "colab" and cmd[1] == "download":
            dest = cmd[-1]
            if dest.endswith(".tar.gz"):
                import tarfile

                with tarfile.open(dest, "w:gz"):
                    pass
        return 0

    monkeypatch.setattr("subprocess.call", mock_subprocess_call)

    # Fast single-file download path (testing default session='cmda')
    dummy_file = "data/apo/analysis/WT/02_eq/qc.png"
    status = colab_pull(
        root=tmp_path,
        project="kras-example",
        path=dummy_file,
    )
    assert status == 0
    assert len(calls) == 1
    assert calls[0][0] == "colab"
    assert calls[0][1] == "download"
    assert "-s" in calls[0]
    assert "cmda" in calls[0]
    assert (
        "/content/drive/MyDrive/research/projects/kras-example/data/apo/analysis/WT/02_eq/qc.png"
        in calls[0]
    )

    # Filtered packaging path
    exec_calls = []

    def mock_execute_remote(session, script, timeout=300):
        exec_calls.append((session, script))
        return 0

    monkeypatch.setattr("colabmda.remote.colab._execute_remote", mock_execute_remote)

    status_filter = colab_pull(
        session="cmda",
        root=tmp_path,
        project="kras-example",
        stage="eq",
        variant="WT",
    )
    assert status_filter == 0
    assert len(exec_calls) >= 1
    packaging_script = exec_calls[0][1]
    assert "stage = 'eq'" in packaging_script
    assert "variant = 'WT'" in packaging_script
    assert "target_stages" in packaging_script
