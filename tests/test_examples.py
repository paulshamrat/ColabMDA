"""Safe project creation and generated command forwarding without a local GPU."""

import hashlib
import json
import os
import subprocess
import sys

import pytest

from colabmda.core.examples import install_example
from colabmda.core.project import ProjectPaths
from colabmda.core.taxonomy import SystemSpec


def test_example_install_is_repeatable_and_preserves_user_files(tmp_path, monkeypatch):
    root = tmp_path / "project with spaces"
    install_example("trpcage", root)
    receipt = json.loads((root / "data/00_setup/example.json").read_text())
    assert len(receipt["files"]) == 7
    for name, expected in receipt["files"].items():
        assert hashlib.sha256((root / name).read_bytes()).hexdigest() == expected
    install_example("trpcage", root)
    assert not (root / "config.yaml").exists()
    monkeypatch.setenv("COLABMDA_RUN_ID", "run01")
    paths = ProjectPaths(root)
    spec = SystemSpec.from_tag("apo", variant="WT")
    assert paths.simulation_base == root / "data/02_simulation/run01"
    assert paths.analysis_dir(spec).name == "04_analysis"
    changed = root / "scripts/02_simulation/03_production.sh"
    changed.write_text("User settings\n")
    with pytest.raises(ValueError, match="preserved"):
        install_example("trpcage", root)
    assert changed.read_text() == "User settings\n"


def test_example_conflict_preflight_and_symlink_protection(tmp_path):
    root = tmp_path / "project"
    conflict = root / "scripts/02_simulation/05_check.py"
    conflict.parent.mkdir(parents=True)
    conflict.write_text("User code")
    with pytest.raises(ValueError, match="preserved"):
        install_example("trpcage", root)
    assert not (root / "scripts/01_structure").exists()
    linked = tmp_path / "linked"
    linked.mkdir()
    (linked / "scripts").symlink_to(root / "scripts", target_is_directory=True)
    with pytest.raises(ValueError, match="Symlink"):
        install_example("trpcage", linked)
    metadata_link = tmp_path / "metadata_link"
    metadata_link.mkdir()
    outside = tmp_path / "outside"
    outside.mkdir()
    (metadata_link / "data").symlink_to(outside, target_is_directory=True)
    with pytest.raises(ValueError, match="Symlink"):
        install_example("trpcage", metadata_link)
    assert not (metadata_link / "scripts").exists()
    assert not list(outside.iterdir())


def test_generated_production_command_and_duration_record(tmp_path):
    root = tmp_path / "example project"
    install_example("trpcage", root)
    inputs = root / "data/01_structure/trpcage"
    inputs.mkdir(parents=True)
    (inputs / "WT.pdb").write_bytes(b"fixture")
    (inputs / "input.json").write_text(
        json.dumps(
            {
                "input": "data/01_structure/trpcage/WT.pdb",
                "input_sha256": hashlib.sha256(b"fixture").hexdigest(),
            }
        )
    )
    system = root / "data/02_simulation/test/apo/WT"
    (system / "02_eq").mkdir(parents=True)
    (system / "02_eq/EQ_DONE").write_text("test")
    (system / "03_md/r1").mkdir(parents=True)
    binary = tmp_path / "cmda"
    binary.write_text(
        f"#!{sys.executable}\nimport os,json,sys\nfrom pathlib import Path\nPath(os.environ['CALLS']).write_text(json.dumps(sys.argv[1:]))\n(Path(os.environ['SYSTEM_DIR'])/'03_md/r1/prod.chk').write_bytes(b'checkpoint')\n"
    )
    binary.chmod(0o755)
    calls = tmp_path / "calls.json"
    env = {
        **os.environ,
        "PATH": str(tmp_path) + os.pathsep + os.environ["PATH"],
        "CALLS": str(calls),
        "COLABMDA_RUN_ID": "test",
        "COLABMDA_EXAMPLE_NS": "0.1",
    }
    subprocess.run(
        ["bash", str(root / "scripts/02_simulation/03_production.sh")], env=env, check=True
    )
    args = json.loads(calls.read_text())
    assert args[0] == "sim"
    for key, value in {
        "--root": str(root),
        "--state": "apo",
        "--stage": "md",
        "--total-ns": "0.1",
        "--checkpoint-ps": "100",
        "--traj-interval": "10",
    }.items():
        assert args[args.index(key) + 1] == value
    assert (
        json.loads((system.parents[1] / "production.json").read_text())[-1]["requested_ns"] == 0.1
    )
