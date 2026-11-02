"""Public source exports and installed wheel resources work without research files."""

import os
import subprocess
import sys
import zipfile
from pathlib import Path

import pytest

from colabmda.remote.colab import _make_source_zip


@pytest.mark.parametrize("branch", ["main", "dev1", None])
def test_committed_exports_allow_public_branches_and_release_commits(tmp_path, branch):
    repo = tmp_path / "repo"
    repo.mkdir()
    subprocess.run(
        ["git", "init", "-b", branch or "main", str(repo)], check=True, capture_output=True
    )
    for name, value in [("user.email", "test@example.invalid"), ("user.name", "Test")]:
        subprocess.run(["git", "-C", str(repo), "config", name, value], check=True)
    (repo / "README.md").write_text("Public guide\n")
    (repo / "AGENTS.md").write_text("Private instructions\n")
    (repo / "data").mkdir()
    (repo / "data/private.txt").write_text("Private runtime file\n")
    subprocess.run(["git", "-C", str(repo), "add", "README.md", "AGENTS.md"], check=True)
    subprocess.run(
        ["git", "-C", str(repo), "commit", "-m", "Synthetic fixture"],
        check=True,
        capture_output=True,
    )
    if branch is None:
        subprocess.run(
            ["git", "-C", str(repo), "checkout", "--detach"], check=True, capture_output=True
        )
    archive = _make_source_zip(repo)
    try:
        with zipfile.ZipFile(archive) as bundle:
            assert set(bundle.namelist()) == {"README.md", "SNAPSHOT.json"}
    finally:
        archive.unlink()
    (repo / "README.md").write_text("Uncommitted runtime change\n")
    with pytest.raises(RuntimeError, match="Commit runtime changes"):
        _make_source_zip(repo)


def test_wheel_retains_nucleotide_resources_without_editable_checkout(tmp_path):
    pytest.importorskip("build")
    repo = Path(__file__).resolve().parents[1]
    source = tmp_path / "source"
    source.mkdir()
    # Export only tracked build inputs; avoid old egg-info masking missing resources.
    tracked = subprocess.check_output(["git", "-C", str(repo), "ls-files"], text=True).splitlines()
    for name in tracked:
        original = repo / name
        if original.is_file():
            target = source / name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(original.read_bytes())
    subprocess.run(
        [
            sys.executable,
            "-m",
            "build",
            "--wheel",
            "--no-isolation",
            "--outdir",
            str(tmp_path / "dist"),
        ],
        cwd=source,
        check=True,
        capture_output=True,
    )
    wheel = next((tmp_path / "dist").glob("*.whl"))
    with zipfile.ZipFile(wheel) as archive:
        expected = {
            str(p.relative_to(repo / "src"))
            for p in (repo / "src/colabmda/amber/params").glob("*")
            if p.suffix in {".frcmod", ".prepi"}
        }
        assert len(expected) == 7
        assert expected <= set(archive.namelist())
        installed = tmp_path / "installed"
        archive.extractall(installed)
    # -I/-S prevents the checkout or editable finder from satisfying imports.
    code = "import sys;sys.path.insert(0,sys.argv[1]);from colabmda.amber.leap import get_param_path;assert get_param_path('ATP-B3.prepi').is_file();assert get_param_path('ATP-B3.frcmod').is_file()"
    subprocess.run(
        [sys.executable, "-I", "-S", "-c", code, str(installed)],
        cwd=tmp_path,
        check=True,
        env={k: v for k, v in os.environ.items() if k != "PYTHONPATH"},
    )
