"""Flat source layout, pure calculations, pilot paths and command forwarding."""

import ast
import importlib
from pathlib import Path

import pytest

from colabmda.core.project import ProjectPaths
from colabmda.core.taxonomy import SystemSpec


def test_flat_imports_and_no_plotting_in_calculations():
    for name in [
        "prepare",
        "em",
        "nvt",
        "npt",
        "md",
        "qc",
        "merge",
        "analysis",
        "common",
        "run",
        "figures",
    ]:
        module = importlib.import_module("colabmda.openmm." + name)
        if name != "figures":
            tree = ast.parse(Path(module.__file__).read_text())
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    assert all(not a.name.startswith(("matplotlib", "seaborn")) for a in node.names)
                if isinstance(node, ast.ImportFrom):
                    assert not (node.module or "").startswith(("matplotlib", "seaborn"))


def test_pilot_paths_do_not_create_parallel_layouts(tmp_path, monkeypatch):
    monkeypatch.setenv("COLABMDA_RUN_ID", "pilot_01")
    p = ProjectPaths(tmp_path)
    p.init_layout()
    spec = SystemSpec.from_tag("apo", variant="WT")
    assert p.prep_dir(spec) == tmp_path / "data/02_pilot/pilot_01/apo/WT/01_prep"
    assert p.eq_dir(spec).name == "02_eq"
    assert p.replica_dir(spec).parent.name == "03_md"
    assert p.analysis_dir(spec).name == "04_check"
    assert not (tmp_path / "simulation").exists()
    assert not (tmp_path / "01_structure").exists()
    monkeypatch.setenv("COLABMDA_RUN_ID", "../bad")
    with pytest.raises(ValueError):
        ProjectPaths(tmp_path)


def test_colab_flags_reach_manager(monkeypatch):
    import sys

    from colabmda.cli import main
    from colabmda.remote import colab

    seen = []
    monkeypatch.setattr(colab, "colab_manager_cli", lambda args: seen.extend(args) or 0)
    monkeypatch.setattr(sys, "argv", ["cmda", "colab", "setup", "-s", "pilot", "--local-source"])
    with pytest.raises(SystemExit) as result:
        main()
    assert result.value.code == 0
    assert seen == ["setup", "-s", "pilot", "--local-source"]


def test_remote_kernel_error_is_not_reported_as_success(monkeypatch):
    from unittest.mock import MagicMock

    from colabmda.remote import colab

    fake = MagicMock()
    fake.stdout = iter(["RuntimeError: installation failed\n"])
    fake.wait.return_value = 0
    monkeypatch.setattr(colab.subprocess, "Popen", lambda *a, **k: fake)
    assert colab._execute_remote("test", "raise RuntimeError('failed')", 10) == 1


def test_setup_does_not_wait_for_optional_interactive_installs(tmp_path, monkeypatch):
    from colabmda.remote import colab

    archive = tmp_path / "source.zip"
    archive.write_bytes(b"test")
    seen = []
    monkeypatch.setattr(colab, "_make_source_zip", lambda root: archive)
    monkeypatch.setattr(colab.subprocess, "check_call", lambda args: 0)
    monkeypatch.setattr(
        colab, "_execute_remote", lambda session, code, timeout: seen.append(code) or 0
    )
    assert colab.colab_manager_cli(["setup", "-s", "test"]) == 0
    assert "setup.sh </dev/null" in seen[0]
    assert not archive.exists()


def test_installed_package_exposes_its_declared_version():
    from importlib.metadata import version

    import colabmda

    assert colabmda.__version__ == version("colabmda")
