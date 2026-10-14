#!/usr/bin/env python3
"""Run the ColabMDA bootstrap inside an existing Colab session via colab exec.

This is a local helper. It generates a temporary Python file that `colab exec`
sends to the remote Colab kernel. The remote Python code then launches the
shell bootstrap and streams stdout/stderr back through `colab exec`, while also
writing `/content/bootstrap_colabmda.log` on the VM.
"""

from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import tempfile
import zipfile
from pathlib import Path


def _bool_env(value: bool) -> str:
    return "1" if value else "0"


def _should_zip(path: Path, repo_root: Path) -> bool:
    rel = path.relative_to(repo_root)
    parts = set(rel.parts)
    if parts & {".git", "data", ".ruff_cache", ".agents", ".codex", "__pycache__"}:
        return False
    if rel.parts[:2] == ("docs", "_build"):
        return False
    if path.name == "ColabMDA-upload.zip":
        return False
    if path.suffix in {".pyc", ".pyo", ".dcd", ".xtc", ".trr", ".nc", ".chk", ".log"}:
        return False
    return True


def _make_source_zip(repo_root: Path) -> Path:
    tmp = tempfile.NamedTemporaryFile(suffix=".zip", delete=False)
    tmp.close()
    zip_path = Path(tmp.name)
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in repo_root.rglob("*"):
            if path.is_file() and _should_zip(path, repo_root):
                archive.write(path, path.relative_to(repo_root))
    return zip_path


def main() -> int:
    repo_root = Path(__file__).resolve().parents[1]

    parser = argparse.ArgumentParser(
        description="Bootstrap a running Colab session using colab exec."
    )
    parser.add_argument("-s", "--session", default="kras-sim", help="colab-cli session name")
    parser.add_argument(
        "--branch",
        default="pwpl",
        help="Git branch/ref to download bootstrap_colab_openmm_gpu.sh from",
    )
    parser.add_argument(
        "--repo",
        default="paulshamrat/ColabMDA",
        help="GitHub repository in owner/name form",
    )
    parser.add_argument(
        "--with-modeller",
        action="store_true",
        help="Install MODELLER environment. Requires KEY_MODELLER or MODELLER_LICENSE.",
    )
    parser.add_argument(
        "--with-ligand",
        action="store_true",
        help="Install OpenFF/openmmforcefields ligand parameterization tools.",
    )
    parser.add_argument(
        "--with-docking",
        action="store_true",
        help="Install RDKit, Meeko, and AutoDock Vina for docking/virtual screening.",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=7200,
        help="colab exec timeout in seconds",
    )
    parser.add_argument(
        "--key-modeller",
        default=None,
        help="MODELLER license key. Prefer setting KEY_MODELLER in your shell.",
    )
    parser.add_argument(
        "--local-bootstrap",
        action="store_true",
        help="Use scripts/bootstrap_colab_openmm_gpu.sh from this local checkout.",
    )
    parser.add_argument(
        "--local-source",
        action="store_true",
        help="Upload this local checkout and install it editable after bootstrap.",
    )
    parser.add_argument(
        "--from-local",
        action="store_true",
        help="Shortcut for --local-bootstrap --local-source.",
    )
    args = parser.parse_args()

    if args.from_local:
        args.local_bootstrap = True
        args.local_source = True

    key_modeller = ""
    if args.with_modeller:
        key_modeller = (
            args.key_modeller
            or os.environ.get("KEY_MODELLER")
            or os.environ.get("MODELLER_LICENSE", "")
        )

    if args.with_modeller and not key_modeller:
        raise SystemExit(
            "ERROR: --with-modeller requires a MODELLER license key.\n"
            "Set it locally first, e.g.:\n"
            "  export KEY_MODELLER='YOUR_MODELLER_LICENSE_KEY'\n"
            "Then rerun this command."
        )

    bootstrap_lines = ["cd /content"]
    bootstrap_text = ""
    if args.local_bootstrap:
        bootstrap_path = repo_root / "scripts" / "bootstrap_colab_openmm_gpu.sh"
        bootstrap_text = bootstrap_path.read_text()
        bootstrap_lines.append("chmod +x /content/bootstrap_colab_openmm_gpu.sh")
    else:
        bootstrap_url = (
            f"https://raw.githubusercontent.com/{args.repo}/{args.branch}/"
            "scripts/bootstrap_colab_openmm_gpu.sh"
        )
        bootstrap_lines.extend(
            [
                f"curl -fsSL {shlex.quote(bootstrap_url)} -o bootstrap_colab_openmm_gpu.sh",
                "chmod +x bootstrap_colab_openmm_gpu.sh",
            ]
        )

    remote_shell_lines = [
        "set -euo pipefail",
        *bootstrap_lines,
        f"export WITH_MODELLER={_bool_env(args.with_modeller)}",
        f"export WITH_LIGAND={_bool_env(args.with_ligand)}",
        f"export WITH_DOCKING={_bool_env(args.with_docking)}",
        f"export INSTALL_REF={shlex.quote(args.branch)}",
        "bash ./bootstrap_colab_openmm_gpu.sh 2>&1 | tee /content/bootstrap_colabmda.log",
    ]

    if args.local_source:
        remote_shell_lines.extend(
            [
                "rm -rf /content/ColabMDA",
                "mkdir -p /content/ColabMDA",
                "unzip -q /content/ColabMDA-upload.zip -d /content/ColabMDA",
                'source "$HOME/miniforge3/etc/profile.d/conda.sh"',
                "conda activate openmm_env",
                "python -m pip install -e /content/ColabMDA",
                'if [[ "${WITH_MODELLER}" == "1" ]]; then',
                "  conda activate modeller_env",
                "  python -m pip install -e /content/ColabMDA",
                "  conda activate openmm_env",
                "fi",
                'if [[ "${WITH_DOCKING}" == "1" ]]; then',
                "  conda activate docking_env",
                "  python -m pip install -e /content/ColabMDA",
                "  conda activate openmm_env",
                "fi",
                "colabmda --help >/dev/null",
                'echo "[OK] Local checkout installed editable from /content/ColabMDA"',
            ]
        )

    remote_shell = "\n".join(remote_shell_lines)

    remote_py = (
        "import pathlib\n"
        "import os\n"
        "import subprocess\n"
        f"bootstrap_text = {bootstrap_text!r}\n"
        "if bootstrap_text:\n"
        "    pathlib.Path('/content/bootstrap_colab_openmm_gpu.sh').write_text(bootstrap_text)\n"
        "env = os.environ.copy()\n"
        f"key_modeller = {key_modeller!r}\n"
        "if key_modeller:\n"
        "    env['KEY_MODELLER'] = key_modeller\n"
        f"cmd = {remote_shell!r}\n"
        '_result = subprocess.run(["bash", "-lc", cmd], check=True, env=env)\n'
        "print('[OK] Remote bootstrap command finished with returncode', _result.returncode)\n"
    )

    source_zip: Path | None = None
    with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False) as handle:
        handle.write(remote_py)
        tmp_path = Path(handle.name)

    try:
        if args.local_source:
            source_zip = _make_source_zip(repo_root)
            upload_cmd = [
                "colab",
                "upload",
                "-s",
                args.session,
                str(source_zip),
                "/content/ColabMDA-upload.zip",
            ]
            print("[INFO] Uploading local source zip to /content/ColabMDA-upload.zip")
            subprocess.check_call(upload_cmd)

        cmd = [
            "colab",
            "exec",
            "-s",
            args.session,
            "-f",
            str(tmp_path),
            "--timeout",
            str(args.timeout),
        ]
        print("[INFO] Running:", " ".join(shlex.quote(part) for part in cmd))
        print("[INFO] Remote log: /content/bootstrap_colabmda.log")
        return subprocess.call(cmd)
    finally:
        tmp_path.unlink(missing_ok=True)
        if source_zip is not None:
            source_zip.unlink(missing_ok=True)


if __name__ == "__main__":
    raise SystemExit(main())
