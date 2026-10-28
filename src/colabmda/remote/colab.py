from __future__ import annotations

import argparse
import getpass
import hashlib
import json
import os
import shlex
import subprocess
import sys
import tempfile
import uuid
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

    git("rev-parse", "--verify", "HEAD")  # Any clean committed branch or release tag.
    allowed = [
        "src",
        "examples",
        "scripts/setup.sh",
        "pyproject.toml",
        "README.md",
        "docs",
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


def _execute_remote(session, code, timeout):
    """Stream output and require a success sentinel; kernel errors can exit CLI zero."""
    marker = "COLABMDA_REMOTE_OK_" + uuid.uuid4().hex
    command = ["colab", "exec", *_session_args(session), "--timeout", str(timeout)]
    process = subprocess.Popen(
        command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
    )
    process.stdin.write(code + f"\nprint({marker!r})\n")
    process.stdin.close()
    succeeded = False
    for line in process.stdout:
        print(line, end="", flush=True)
        if line.strip() == marker:
            succeeded = True
    status = process.wait()
    return status if status else (0 if succeeded else 1)


def colab_push(
    session: str | None = None,
    root: str | Path | None = None,
    project: str | None = None,
    timeout: float = 300,
) -> int:
    """Package local data/ and push directly into remote Colab Google Drive."""
    import tarfile

    local_root = Path(root or os.getcwd()).resolve()
    data_dir = local_root / "data"
    if not data_dir.is_dir():
        print(f"[ERROR] No 'data/' directory found in {local_root}")
        print("Run local modeling first (cmda fetch, cmda model) before pushing to Colab.")
        return 1

    project_name = project or local_root.name
    if project_name in {"research", "projects", "ColabMDA"}:
        project_name = "cmda-example"
    remote_base = f"/content/drive/MyDrive/research/projects/{project_name}"

    print(f"[INFO] Packaging local data/ from {local_root}...")
    tmp_tar = tempfile.NamedTemporaryFile(suffix=".tar.gz", delete=False)
    tmp_tar.close()
    try:
        with tarfile.open(tmp_tar.name, "w:gz") as tar:
            tar.add(str(data_dir), arcname="data")

        session_name = session or "active"
        print(f"[INFO] Uploading data archive to Colab session '{session_name}'...")
        remote_tar = "/content/_cmda_push.tar.gz"
        upload_cmd = ["colab", "upload", *_session_args(session), tmp_tar.name, remote_tar]
        res = subprocess.call(upload_cmd)
        if res != 0:
            print("[ERROR] colab upload failed. Ensure the session is running and authenticated.")
            return res

        print(f"[INFO] Extracting directly into Google Drive: {remote_base}/data/ ...")
        remote_script = (
            "import os, tarfile\n"
            "drive_base = '/content/drive/MyDrive'\n"
            "if not os.path.exists(drive_base):\n"
            "    raise RuntimeError('Google Drive is not mounted! Run: colab drivemount --session <name>')\n"
            f"dest = {remote_base!r}\n"
            "os.makedirs(dest, exist_ok=True)\n"
            f"with tarfile.open({remote_tar!r}, 'r:gz') as tar:\n"
            "    tar.extractall(dest)\n"
            f"if os.path.exists({remote_tar!r}):\n"
            f"    os.remove({remote_tar!r})\n"
            f"print('[OK] Synchronized data/ to {{dest}}/data')\n"
        )
        status = _execute_remote(session, remote_script, timeout=timeout)
        if status == 0:
            print("✔ Successfully pushed local models and structures to Google Drive:")
            print(f"  • Colab Linux Path:   {remote_base}/data/")
            print(
                f"  • Google Drive Web:   My Drive > research > projects > {project_name} > data/"
            )
        return status
    finally:
        Path(tmp_tar.name).unlink(missing_ok=True)


def colab_pull(
    session: str | None = None,
    root: str | Path | None = None,
    project: str | None = None,
    all_files: bool = False,
    timeout: float = 300,
) -> int:
    """Download analysis figures and reports (or full data) from Google Drive back locally."""
    import tarfile

    local_root = Path(root or os.getcwd()).resolve()
    project_name = project or local_root.name
    if project_name in {"research", "projects", "ColabMDA"}:
        project_name = "cmda-example"
    remote_base = f"/content/drive/MyDrive/research/projects/{project_name}"
    remote_tar = "/content/_cmda_pull.tar.gz"

    session_name = session or "active"
    print(f"[INFO] Packaging results from Colab session '{session_name}' ({remote_base}/data)...")

    exclude_logic = (
        ""
        if all_files
        else "def filter_heavy(tarinfo):\n"
        "    if tarinfo.name.endswith(('.dcd', '.chk', '.inpcrd')):\n"
        "        return None\n"
        "    return tarinfo\n"
    )
    filter_arg = "" if all_files else ", filter=filter_heavy"

    remote_script = (
        "import os, tarfile\n"
        f"src = os.path.join({remote_base!r}, 'data')\n"
        "if not os.path.exists(src):\n"
        "    raise FileNotFoundError(f'Remote directory not found: {src}')\n"
        f"{exclude_logic}"
        f"with tarfile.open({remote_tar!r}, 'w:gz') as tar:\n"
        f"    tar.add(src, arcname='data'{filter_arg})\n"
        f"print('[OK] Results packaged to {remote_tar}')\n"
    )
    status = _execute_remote(session, remote_script, timeout=timeout)
    if status != 0:
        print("[ERROR] Failed to package remote results.")
        return status

    tmp_tar = tempfile.NamedTemporaryFile(suffix=".tar.gz", delete=False)
    tmp_tar.close()
    try:
        print(f"[INFO] Downloading results from Colab session '{session_name}'...")
        download_cmd = ["colab", "download", *_session_args(session), remote_tar, tmp_tar.name]
        res = subprocess.call(download_cmd)
        if res != 0:
            print("[ERROR] colab download failed.")
            return res

        print(f"[INFO] Extracting results into {local_root}/data/ ...")
        with tarfile.open(tmp_tar.name, "r:gz") as tar:
            tar.extractall(local_root)

        cleanup_cmd = f"import os\nif os.path.exists({remote_tar!r}): os.remove({remote_tar!r})\n"
        _execute_remote(session, cleanup_cmd, timeout=30)
        print(f"✔ Successfully pulled results to {local_root}/data/")
        print(f"  • Local Path:         {local_root}/data/")
        print(f"  • Retrieved From:     My Drive > research > projects > {project_name} > data/")
        return 0
    finally:
        Path(tmp_tar.name).unlink(missing_ok=True)


def colab_manager_cli(argv: list[str] | None = None) -> int:
    """CLI handler for `cmda colab`."""
    repo_root = Path(__file__).resolve().parents[3]
    if not (repo_root / "scripts" / "setup.sh").is_file():
        for candidate in [Path.cwd(), Path.home() / "ColabMDA", Path(__file__).resolve().parents[4]]:
            if (candidate / "scripts" / "setup.sh").is_file():
                repo_root = candidate
                break

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

    p_boot.add_argument(
        "--modeller",
        action="store_true",
        help="Install MODELLER; securely prompt for a key if needed",
    )
    p_boot.add_argument(
        "--no-modeller", action="store_true", help="Skip the optional interactive MODELLER question"
    )

    # Push data/ to Google Drive
    p_push = sub.add_parser("push", help="Push local data/ to remote Google Drive via Colab")
    p_push.add_argument(
        "-s",
        "--session",
        default=None,
        help="Session name (default: auto-detected if omitted)",
    )
    p_push.add_argument(
        "--root",
        default=None,
        help="Local project root (default: current working directory)",
    )
    p_push.add_argument(
        "--project",
        default=None,
        help="Remote project folder name in Google Drive (default: directory name or cmda-example)",
    )
    p_push.add_argument(
        "--timeout",
        type=float,
        default=300,
        help="Push timeout in seconds (default: 300)",
    )

    # Pull results from Google Drive
    p_pull = sub.add_parser("pull", help="Pull results from remote Google Drive via Colab")
    p_pull.add_argument(
        "-s",
        "--session",
        default=None,
        help="Session name (default: auto-detected if omitted)",
    )
    p_pull.add_argument(
        "--root",
        default=None,
        help="Local project root (default: current working directory)",
    )
    p_pull.add_argument(
        "--project",
        default=None,
        help="Remote project folder name in Google Drive (default: directory name or cmda-example)",
    )
    p_pull.add_argument(
        "--all",
        dest="all_files",
        action="store_true",
        help="Download all files including heavy trajectories (.dcd)",
    )
    p_pull.add_argument(
        "--timeout",
        type=float,
        default=300,
        help="Pull timeout in seconds (default: 300)",
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
        if "\n" in full_cmd and not full_cmd.startswith(("%", "import", "print")):
            # A notebook ! escape handles only one line. IPython.system streams
            # output; a Bash cell magic buffers it until the long job finishes.
            shell = full_cmd.removeprefix("!")
            command = "bash -lc " + shlex.quote(shell)
            code = (
                f"get_ipython().system({command!r})\n"
                "if get_ipython().user_ns.get('_exit_code', 0) != 0: "
                "raise RuntimeError('Remote command failed')"
            )
            return _execute_remote(args.session, code, args.timeout)
        if (
            not full_cmd.startswith("!")
            and not full_cmd.startswith("%")
            and not full_cmd.startswith("import")
            and not full_cmd.startswith("print")
        ):
            full_cmd = f"!{full_cmd}"
        if full_cmd.startswith(("!", "%")):
            full_cmd += "\nif get_ipython().user_ns.get('_exit_code', 0) != 0: raise RuntimeError('Remote command failed')"
        return _execute_remote(args.session, full_cmd, args.timeout)

    if args.colab_action == "setup":
        setup_script = repo_root / "scripts" / "setup.sh"
        if not setup_script.exists():
            raise FileNotFoundError(f"Setup script not found: {setup_script}")

        source_zip: Path | None = None
        license_path = None
        use_modeller = args.modeller or os.environ.get("WITH_MODELLER") == "1"
        if args.modeller and args.no_modeller:
            raise ValueError("Choose --modeller or --no-modeller")
        if not use_modeller and not args.no_modeller and sys.stdin.isatty():
            use_modeller = input("Install MODELLER as well? [y/N]: ").strip().lower() in {
                "y",
                "yes",
            }
        key = None
        if use_modeller:
            key = os.environ.get("KEY_MODELLER") or os.environ.get("MODELLER_LICENSE")
            if not key and sys.stdin.isatty():
                key = getpass.getpass("Enter MODELLER license key: ")
            if not key:
                raise ValueError("MODELLER requires KEY_MODELLER or an interactive key prompt")
        try:
            if key:
                with tempfile.NamedTemporaryFile(
                    mode="w", prefix="colabmda-license-", delete=False
                ) as secret:
                    secret.write(key)
                    license_path = Path(secret.name)
                subprocess.check_call(
                    [
                        "colab",
                        "upload",
                        *_session_args(args.session),
                        str(license_path),
                        "/content/.colabmda-license",
                    ]
                )
                license_path.unlink()
                license_path = None
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

            license_setup = ""
            if use_modeller:
                license_setup = "chmod 600 /content/.colabmda-license; export WITH_MODELLER=1; export KEY_MODELLER=\"$(cat /content/.colabmda-license)\"; trap 'rm -f /content/.colabmda-license' EXIT; "
            remote_cmd = (
                "set -e; " + license_setup + "if [ -f /content/ColabMDA-upload.zip ]; then "
                "mkdir -p /content/ColabMDA; "
                "unzip -oq /content/ColabMDA-upload.zip -d /content/ColabMDA; "
                "fi; "
                "chmod +x /content/setup.sh; "
                "COLABMDA_SOURCE=/content/ColabMDA /content/setup.sh </dev/null"
            )
            session_desc = args.session if args.session else "active"
            print(f"[INFO] Executing setup on {session_desc} session...")
            code = f"!{remote_cmd}\nif get_ipython().user_ns.get('_exit_code', 0) != 0: raise RuntimeError('Remote setup failed')"
            return _execute_remote(args.session, code, args.timeout)
        finally:
            if license_path is not None:
                license_path.unlink(missing_ok=True)
            if source_zip is not None:
                source_zip.unlink(missing_ok=True)

    if args.colab_action == "push":
        return colab_push(
            session=args.session,
            root=args.root,
            project=args.project,
            timeout=args.timeout,
        )

    if args.colab_action == "pull":
        return colab_pull(
            session=args.session,
            root=args.root,
            project=args.project,
            all_files=args.all_files,
            timeout=args.timeout,
        )

    return 0
