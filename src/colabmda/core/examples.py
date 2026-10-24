"""Install runnable example scripts without overwriting user files."""

import hashlib
import json
import shutil
import subprocess
import sys
from importlib.metadata import PackageNotFoundError, distribution, version
from pathlib import Path


def templates(name):
    if name != "trpcage":
        raise ValueError(f"Unknown example: {name}")
    repo = Path(__file__).resolve().parents[3]
    source = repo / "examples" / name
    if source.is_dir():
        return {
            str(p.relative_to(source)): p
            for p in source.rglob("*")
            if p.is_file() and p.suffix in {".py", ".sh"}
        }
    installed = distribution("colabmda")
    prefix = f"share/colabmda/examples/{name}/"
    result = {}
    for item in installed.files or []:
        text = str(item)
        if prefix in text and text.endswith((".py", ".sh")):
            result[text.split(prefix, 1)[1]] = Path(installed.locate_file(item))
    if not result:
        raise RuntimeError("Example resources missing; reinstall ColabMDA")
    return result


def install_example(name, root):
    root = Path(root).expanduser().resolve()
    files = {n: p.read_bytes() for n, p in templates(name).items()}
    receipt = root / "data/00_setup/example.json"
    # Check script, metadata and log destinations before writing anything.
    destinations = [root / name for name in files] + [receipt, root / "logs/02_simulation"]
    for target in destinations:
        for parent in [target, *target.parents]:
            if parent.is_symlink():
                raise ValueError(f"Symlink destination preserved: {parent}")
            if parent != target and parent.exists() and not parent.is_dir():
                raise ValueError(f"Existing non-directory preserved: {parent}")
    if (root / "logs/02_simulation").exists() and not (root / "logs/02_simulation").is_dir():
        raise ValueError("Existing log destination preserved")
    for name_, payload in files.items():
        target = root / name_
        if any(p.is_symlink() for p in [target, *target.parents] if p != root.parent):
            raise ValueError(f"Symlink destination preserved: {target}")
        if target.exists() and (not target.is_file() or target.read_bytes() != payload):
            raise ValueError(f"Existing different file preserved: {target}")
    if receipt.exists() and json.loads(receipt.read_text()).get("example") != name:
        raise ValueError("Existing different example preserved")
    record = {
        "example": name,
        "layout": "numbered",
        "tool_version": version("colabmda"),
        "files": {n: hashlib.sha256(p).hexdigest() for n, p in files.items()},
    }
    repo = Path(__file__).resolve().parents[3]
    snapshot = repo / "SNAPSHOT.json"
    if snapshot.exists():
        record["tool"] = json.loads(snapshot.read_text())
    elif (repo / ".git").is_dir():
        record["source_commit"] = subprocess.check_output(
            ["git", "-C", str(repo), "rev-parse", "HEAD"], text=True
        ).strip()
    for name_, payload in files.items():
        target = root / name_
        target.parent.mkdir(parents=True, exist_ok=True)
        if not target.exists():
            target.write_bytes(payload)
    receipt.parent.mkdir(parents=True, exist_ok=True)
    receipt.write_text(json.dumps(record, indent=2) + "\n")
    (root / "logs/02_simulation").mkdir(parents=True, exist_ok=True)
    print(f"Installed {name} example into {root}; follow the ColabMDA README")
    return 0


def record_run(root, run_id):
    """Record immutable input, runtime and executed-script fingerprints."""
    import re

    root = Path(root).resolve()
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", run_id):
        raise ValueError("Invalid run ID")
    setup = json.loads((root / "data/00_setup/example.json").read_text())
    input_record = json.loads((root / "data/01_structure/trpcage/input.json").read_text())
    coordinates = root / input_record["input"]
    if hashlib.sha256(coordinates.read_bytes()).hexdigest() != input_record["input_sha256"]:
        raise ValueError("Structure input changed")
    runtime = {"version": version("colabmda"), "python": sys.version, "packages": {}}
    for package in ["openmm", "numpy", "pandas", "MDAnalysis", "mdtraj", "biopython"]:
        try:
            runtime["packages"][package] = version(package)
        except PackageNotFoundError:
            runtime["packages"][package] = None
    if shutil.which("nvidia-smi"):
        gpu = subprocess.run(
            ["nvidia-smi", "--query-gpu=name,uuid,driver_version", "--format=csv,noheader"],
            capture_output=True,
            text=True,
        )
        runtime["gpu"] = gpu.stdout.strip() if gpu.returncode == 0 else None
    else:
        runtime["gpu"] = None
    repo = Path(__file__).resolve().parents[3]
    if (repo / "SNAPSHOT.json").exists():
        runtime["snapshot"] = json.loads((repo / "SNAPSHOT.json").read_text())
    elif (repo / ".git").is_dir():
        runtime["commit"] = subprocess.check_output(
            ["git", "-C", str(repo), "rev-parse", "HEAD"], text=True
        ).strip()
    record = {
        "input": input_record,
        "setup": setup,
        "runtime": runtime,
        "scripts": {n: hashlib.sha256((root / n).read_bytes()).hexdigest() for n in setup["files"]},
        "protocol": "quick",
        "nvt_ps": 100,
        "npt_ps": 100,
        "step_fs": 2,
        "frame_interval_ps": 10,
        "checkpoint_ps": 100,
    }
    run = root / "data/02_simulation" / run_id
    run.mkdir(parents=True, exist_ok=True)
    path = run / "provenance.json"
    if path.exists() and json.loads(path.read_text()) != record:
        raise ValueError("Run provenance changed; use a new run ID")
    path.write_text(json.dumps(record, indent=2) + "\n")
