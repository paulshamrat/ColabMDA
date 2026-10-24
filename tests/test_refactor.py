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


def test_rmsf_uses_native_author_numbers_and_identity_gate(tmp_path):
    from types import SimpleNamespace

    from colabmda.openmm.analysis import _residue_labels

    reference = tmp_path / "native.pdb"
    names = ["ASN", "LEU", "SER"]
    labels = [(10, " "), (12, "A"), (20, " ")]
    lines = []
    for i, (name, (number, insertion)) in enumerate(zip(names, labels), 1):
        lines.append(
            f"ATOM  {i:5d}  CA  {name} A{number:4d}{insertion}   {float(i):8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           C  "
        )
    reference.write_text("\n".join(lines) + "\nEND\n")
    atoms = [
        SimpleNamespace(residue=SimpleNamespace(name=name, resSeq=i))
        for i, name in enumerate(names)
    ]
    topology = SimpleNamespace(atom=lambda i: atoms[i])
    assert _residue_labels(topology, [0, 1, 2], reference) == [10, "12A", 20]
    atoms[1].residue.name = "ARG"
    with pytest.raises(ValueError, match="identity mismatch"):
        _residue_labels(topology, [0, 1, 2], reference)


def test_modeller_prompt_transfers_key_without_embedding_it_in_remote_code(tmp_path, monkeypatch):
    from colabmda.remote import colab

    archive = tmp_path / "source.zip"
    archive.write_bytes(b"test")
    monkeypatch.delenv("KEY_MODELLER", raising=False)
    monkeypatch.delenv("MODELLER_LICENSE", raising=False)
    monkeypatch.setattr(colab.sys.stdin, "isatty", lambda: True)
    monkeypatch.setattr(colab.getpass, "getpass", lambda prompt: "FAKE_TEST_LICENSE")
    monkeypatch.setattr(colab, "_make_source_zip", lambda root: archive)
    uploads = []

    def upload(args):
        if args[-1] == "/content/.colabmda-license":
            path = Path(args[-2])
            assert path.read_text() == "FAKE_TEST_LICENSE"
            assert path.stat().st_mode & 0o077 == 0
            uploads.append(path)

    monkeypatch.setattr(colab.subprocess, "check_call", upload)
    codes = []
    monkeypatch.setattr(
        colab, "_execute_remote", lambda session, code, timeout: codes.append(code) or 0
    )
    assert colab.colab_manager_cli(["setup", "--modeller", "-s", "test"]) == 0
    assert uploads and not uploads[0].exists()
    assert "FAKE_TEST_LICENSE" not in codes[0]
    assert "trap 'rm -f /content/.colabmda-license' EXIT" in codes[0]
