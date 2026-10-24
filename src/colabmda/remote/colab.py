from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import tempfile
import zipfile
from pathlib import Path


def _bool_env(value: bool) -> str:
    return "1" if value else "0"


def _should_zip(path: Path, repo_root: Path) -> bool:
    rel = path.relative_to(repo_root)
    parts = set(rel.parts)
    if parts & {
        ".git",
        "data",
        ".ruff_cache",
        ".agents",
        ".codex",
        "__pycache__",
        ".pytest_cache",
    }:
        return False
    if rel.parts[:2] == ("docs", "_build"):
        return False
    if path.name == "ColabMDA-upload.zip":
        return False
    if path.suffix in {
        ".pyc",
        ".pyo",
        ".dcd",
        ".xtc",
        ".trr",
        ".nc",
        ".chk",
        ".log",
        ".pdb",
    }:
        return False
    return True


def _make_source_zip(repo_root: Path) -> Path:
    """Export committed runtime files, excluding private instructions and build copies."""

    def git(*args):
        return subprocess.check_output(["git", "-C", str(repo_root), *args])

    if git("branch", "--show-current").decode().strip() != "dev1":
        raise RuntimeError("Local development deployment requires dev1")
    allowed = [
        "src",
        "scripts/setup.sh",
        "pyproject.toml",
        "README.md",
        "LICENSE",
        "environment.yml",
    ]
    if git("status", "--porcelain", "--", *allowed).strip():
        raise RuntimeError("Commit runtime changes before deploying to Colab")
    commit = git("rev-parse", "HEAD").decode().strip()
    files = git("ls-tree", "-r", "--name-only", "HEAD", "--", *allowed).decode().splitlines()
    tmp = tempfile.NamedTemporaryFile(suffix=".zip", delete=False)
    tmp.close()
    record = {"commit": commit, "files": {}}
    with zipfile.ZipFile(tmp.name, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name in files:
            payload = git("show", f"HEAD:{name}")
            archive.writestr(name, payload)
            record["files"][name] = hashlib.sha256(payload).hexdigest()
        archive.writestr("SNAPSHOT.json", json.dumps(record, indent=2) + "\n")
    return Path(tmp.name)


def _session_args(session: str | None) -> list[str]:
    return ["-s", session] if session else []


def colab_manager_cli(argv: list[str] | None = None) -> int:
    """CLI handler for `cmda colab`."""
    repo_root = Path(__file__).resolve().parents[3]

    parser = argparse.ArgumentParser(
        description="Remote Google Colab VM provisioning and execution manager."
    )
    sub = parser.add_subparsers(dest="colab_action")

    # Status / Whoami
    sub.add_parser("whoami", help="Check Google Colab authentication status")

    # Start new VM session
    p_new = sub.add_parser("new", help="Start a new Colab GPU instance")
    p_new.add_argument(
        "-s",
        "--session",
        default=None,
        help="Session name (default: auto-detected if omitted)",
    )
    p_new.add_argument(
        "--gpu",
        default="T4",
        help="Request GPU instance (e.g. T4, L4, A100; default: T4)",
    )

    # Stop session
    p_stop = sub.add_parser("stop", help="Stop and terminate a Colab VM session")
    p_stop.add_argument(
        "-s",
        "--session",
        default=None,
        help="Session name to terminate (default: auto-detected if omitted)",
    )

    # Exec arbitrary command or script
    p_exec = sub.add_parser("exec", help="Execute command or script inside remote Colab VM")
    p_exec.add_argument(
        "-s",
        "--session",
        default=None,
        help="Target session name (default: auto-detected if omitted)",
    )
    p_exec.add_argument(
        "--timeout",
        type=float,
        default=7200,
        help="Timeout in seconds (default: 7200)",
    )
    p_exec.add_argument("cmd", nargs="*", help="Command string to run remotely")

    # Bootstrap setup.sh
    p_boot = sub.add_parser("setup", help="Bootstrap and install ColabMDA on remote VM")
    p_boot.add_argument(
        "-s",
        "--session",
        default=None,
        help="Target session name (default: auto-detected if omitted)",
    )
    p_boot.add_argument(
        "--local-source",
        action="store_true",
        default=True,
        help="Upload local ColabMDA source zip and install editable",
    )
    p_boot.add_argument(
        "--timeout",
        type=float,
        default=7200,
        help="Timeout in seconds (default: 7200)",
    )

    args = parser.parse_args(argv)

    if not args.colab_action:
        parser.print_help()
        return 0

    if args.colab_action == "whoami":
        return subprocess.call(["colab", "whoami"])

    if args.colab_action == "new":
        cmd = ["colab", "new"] + _session_args(args.session)
        if args.gpu:
            cmd.extend(["--gpu", args.gpu])
        return subprocess.call(cmd)

    if args.colab_action == "stop":
        return subprocess.call(["colab", "stop"] + _session_args(args.session))

    if args.colab_action == "exec":
        full_cmd = " ".join(args.cmd) if args.cmd else "!nvidia-smi"
        if (
            not full_cmd.startswith("!")
            and not full_cmd.startswith("%")
            and not full_cmd.startswith("import")
            and not full_cmd.startswith("print")
        ):
            full_cmd = f"!{full_cmd}"
        exec_cmd = ["colab", "exec"] + _session_args(args.session)
        if args.timeout:
            exec_cmd.extend(["--timeout", str(args.timeout)])
        return subprocess.run(
            exec_cmd,
            input=full_cmd,
            text=True,
        ).returncode

    if args.colab_action == "setup":
        setup_script = repo_root / "scripts" / "setup.sh"
        if not setup_script.exists():
            raise FileNotFoundError(f"Setup script not found: {setup_script}")

        source_zip: Path | None = None
        try:
            if args.local_source:
                source_zip = _make_source_zip(repo_root)
                print("[INFO] Uploading local ColabMDA source zip to /content/ColabMDA-upload.zip")
                subprocess.check_call(
                    ["colab", "upload"]
                    + _session_args(args.session)
                    + [
                        str(source_zip),
                        "/content/ColabMDA-upload.zip",
                    ]
                )

            # Upload setup.sh
            print("[INFO] Uploading scripts/setup.sh to /content/setup.sh")
            subprocess.check_call(
                ["colab", "upload"]
                + _session_args(args.session)
                + [
                    str(setup_script),
                    "/content/setup.sh",
                ]
            )

            remote_cmd = (
                "set -e; "
                "if [ -f /content/ColabMDA-upload.zip ]; then "
                "mkdir -p /content/ColabMDA; "
                "unzip -oq /content/ColabMDA-upload.zip -d /content/ColabMDA; "
                "fi; "
                "chmod +x /content/setup.sh; "
                "COLABMDA_SOURCE=/content/ColabMDA /content/setup.sh"
            )
            session_desc = args.session if args.session else "active"
            print(f"[INFO] Executing setup on {session_desc} session...")
            exec_cmd = ["colab", "exec"] + _session_args(args.session)
            if args.timeout:
                exec_cmd.extend(["--timeout", str(args.timeout)])
            return subprocess.run(
                exec_cmd,
                input=f"!{remote_cmd}",
                text=True,
            ).returncode
        finally:
            if source_zip is not None:
                source_zip.unlink(missing_ok=True)

    return 0
