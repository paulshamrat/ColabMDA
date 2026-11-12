"""Missing system tools must stop setup before downloads or environment changes."""

import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
TOOLS = ("git", "curl", "wget")


@pytest.mark.parametrize("missing", [("git",), ("curl",), ("wget",), TOOLS])
def test_setup_reports_missing_tools_before_changing_environment(tmp_path, missing):
    tools = tmp_path / "bin"
    tools.mkdir()
    for name in set(TOOLS) - set(missing):
        stub = tools / name
        stub.write_text("#!/bin/sh\nexit 97\n")
        stub.chmod(0o755)
    prefix = tmp_path / "miniforge3"
    result = subprocess.run(
        [shutil.which("bash"), str(ROOT / "scripts/setup.sh")],
        env={"HOME": str(tmp_path), "PATH": str(tools), "MINIFORGE_DIR": str(prefix)},
        capture_output=True,
        text=True,
        timeout=5,
    )
    assert result.returncode == 1
    assert "Missing required command(s): " + " ".join(missing) in result.stderr
    assert "sudo apt-get install -y git curl wget ca-certificates" in result.stderr
    assert "Then rerun this installer." in result.stderr
    assert "[STEP" not in result.stdout
    assert not prefix.exists()
    assert not (tmp_path / ".bashrc").exists()
