"""Example argument forwarding without executing a local GPU simulation."""

import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "examples/basic/run.sh"


def test_example_simulation_forwards_explicit_project_and_resume_settings(tmp_path):
    binary = tmp_path / "cmda"
    binary.write_text(
        f"#!{sys.executable}\n"
        "import json,os,sys\n"
        "with open(os.environ['EXAMPLE_CALLS'],'a') as out:\n"
        "    out.write(json.dumps({'args':sys.argv[1:],'platform':os.environ['COLABMDA_PLATFORM']})+'\\n')\n"
    )
    binary.chmod(0o755)
    calls = tmp_path / "calls.jsonl"
    env = {
        **os.environ,
        "PATH": str(tmp_path) + os.pathsep + os.environ["PATH"],
        "EXAMPLE_CALLS": str(calls),
        "COLABMDA_RUN_ID": "example_test",
        "COLABMDA_STATE": "atp_mg",
        "COLABMDA_VARIANT": "L119R",
        "COLABMDA_EXAMPLE_NS": "0.1",
        "COLABMDA_PLATFORM": "CUDA",
    }
    project = str(tmp_path / "project with spaces")
    subprocess.run(["bash", str(RUNNER), project, "unused.pdb", "simulate"], env=env, check=True)
    records = [json.loads(line) for line in calls.read_text().splitlines()]
    assert len(records) == 1
    args = records[0]["args"]
    assert args[0] == "sim"
    for option, value in {
        "--root": project,
        "--state": "atp_mg",
        "--variant": "L119R",
        "--stage": "md",
        "--total-ns": "0.1",
        "--checkpoint-ps": "100",
        "--traj-interval": "10",
    }.items():
        assert args[args.index(option) + 1] == value
    assert records[0]["platform"] == "CUDA"


def test_example_requires_run_id_and_rejects_unknown_stage(tmp_path):
    env = dict(os.environ)
    env.pop("COLABMDA_RUN_ID", None)
    result = subprocess.run(
        ["bash", str(RUNNER), str(tmp_path), "unused.pdb", "simulate"],
        env=env,
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
    assert "COLABMDA_RUN_ID" in result.stderr
    env["COLABMDA_RUN_ID"] = "example_test"
    result = subprocess.run(
        ["bash", str(RUNNER), str(tmp_path), "unused.pdb", "unknown"],
        env=env,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 2
    assert "Unknown stage" in result.stderr
