"""Credential hooks must reject real-looking values without printing them."""

import os
import shutil
import subprocess
import uuid
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def git(repo, *args):
    return subprocess.run(["git", "-C", str(repo), *args], text=True, capture_output=True)


@pytest.fixture
def guarded_repo(tmp_path):
    scanner = shutil.which("gitleaks")
    if not scanner:
        pytest.skip("Gitleaks required for credential hook acceptance")
    repo = tmp_path / "repo"
    assert git(tmp_path, "init", "-b", "main", str(repo)).returncode == 0
    for key, value in [("user.name", "Test"), ("user.email", "test@example.invalid")]:
        assert git(repo, "config", key, value).returncode == 0
    (repo / "README.md").write_text("Synthetic repository\n")
    assert git(repo, "add", ".").returncode == 0
    assert git(repo, "commit", "-m", "Initial fixture").returncode == 0
    shutil.copy2(ROOT / ".gitleaks.toml", repo / ".gitleaks.toml")
    shutil.copytree(ROOT / ".githooks", repo / ".githooks")
    assert git(repo, "config", "core.hooksPath", ".githooks").returncode == 0
    assert git(repo, "config", "colabmda.gitleaksPath", scanner).returncode == 0
    return repo


@pytest.mark.parametrize(
    "prefix",
    ['KEY_MODELLER="', 'license = r"', 'MODELLER_LICENSE="', 'SAPISID="', '__Secure-3PAPISID="'],
)
def test_commit_hook_rejects_credential_without_disclosure(guarded_repo, prefix):
    secret = "fixture_" + uuid.uuid4().hex
    (guarded_repo / "fixture.txt").write_text(prefix + secret + '"\n')
    assert git(guarded_repo, "add", "fixture.txt").returncode == 0
    result = git(guarded_repo, "commit", "-m", "Attempted synthetic credential")
    assert result.returncode != 0
    assert secret not in result.stdout + result.stderr
    (guarded_repo / "fixture.txt").write_text(prefix + "YOUR_LICENSE_KEY" + '"\n')
    assert git(guarded_repo, "add", "fixture.txt").returncode == 0
    result = git(guarded_repo, "commit", "-m", "Safe placeholder")
    assert result.returncode == 0, result.stderr


def test_private_push_guard_and_deleted_historical_secret(guarded_repo):
    repo = guarded_repo
    hook = repo / ".githooks/pre-push"
    assert git(repo, "config", "colabmda.privateCheckout", "true").returncode == 0
    result = subprocess.run([str(hook)], cwd=repo, capture_output=True, text=True)
    assert result.returncode != 0
    assert "private research checkout" in result.stderr
    assert git(repo, "config", "colabmda.privateCheckout", "false").returncode == 0
    secret = "fixture_" + uuid.uuid4().hex
    p = repo / "fixture.txt"
    p.write_text('KEY_MODELLER="' + secret + '"\n')
    assert git(repo, "add", "fixture.txt").returncode == 0
    assert (
        git(
            repo, "-c", "core.hooksPath=/dev/null", "commit", "-m", "Synthetic old exposure"
        ).returncode
        == 0
    )
    p.unlink()
    assert git(repo, "add", "-u").returncode == 0
    assert (
        git(
            repo, "-c", "core.hooksPath=/dev/null", "commit", "-m", "Removed from current tree"
        ).returncode
        == 0
    )
    result = subprocess.run(
        [str(hook)], cwd=repo, capture_output=True, text=True, env=os.environ.copy()
    )
    assert result.returncode != 0
    assert secret not in result.stdout + result.stderr


def test_initial_commit_scans_the_staged_index(tmp_path):
    scanner = shutil.which("gitleaks")
    if not scanner:
        pytest.skip("Gitleaks required for credential hook acceptance")
    repo = tmp_path / "initial"
    assert git(tmp_path, "init", "-b", "main", str(repo)).returncode == 0
    for key, value in [
        ("user.name", "Test"),
        ("user.email", "test@example.invalid"),
        ("core.hooksPath", ".githooks"),
        ("colabmda.gitleaksPath", scanner),
    ]:
        assert git(repo, "config", key, value).returncode == 0
    shutil.copytree(ROOT / ".githooks", repo / ".githooks")
    shutil.copy2(ROOT / ".gitleaks.toml", repo / ".gitleaks.toml")
    secret = "fixture_" + uuid.uuid4().hex
    p = repo / "fixture.txt"
    p.write_text('KEY_MODELLER="' + secret + '"\n')
    assert git(repo, "add", ".").returncode == 0
    result = git(repo, "commit", "-m", "Reject first-commit credential")
    assert result.returncode != 0
    assert secret not in result.stdout + result.stderr
    # A working-tree edit must not replace the safe staged snapshot being checked.
    p.write_text('KEY_MODELLER="YOUR_LICENSE_KEY"\n')
    assert git(repo, "add", "fixture.txt").returncode == 0
    p.write_text('KEY_MODELLER="' + secret + '"\n')
    result = git(repo, "commit", "-m", "Safe staged initial commit")
    assert result.returncode == 0, result.stderr
    assert secret not in git(repo, "show", "HEAD:fixture.txt").stdout


def test_installer_keeps_license_out_of_arguments_and_private_files(tmp_path):
    script = (ROOT / "scripts/setup.sh").read_text()
    code = script.split("python - <<'PY_LICENSE'\n", 1)[1].split("\nPY_LICENSE", 1)[0]
    prefix = tmp_path / "environment"
    (prefix / "conda-meta").mkdir(parents=True)
    cfg = prefix / "lib/modeller-fixture/modlib/modeller/config.py"
    cfg.parent.mkdir(parents=True)
    cfg.write_text("license = 'YOUR_LICENSE_KEY'\n")
    secret = "fixture_" + uuid.uuid4().hex
    result = subprocess.run(
        [os.sys.executable, "-c", code],
        env={**os.environ, "CONDA_PREFIX": str(prefix), "KEY_MODELLER": secret},
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, result.stderr
    assert secret not in result.stdout + result.stderr
    assert prefix.stat().st_mode & 0o777 == 0o700
    for p in [cfg, prefix / "conda-meta/state"]:
        assert p.stat().st_mode & 0o777 == 0o600
        assert secret in p.read_text()
    assert "conda env config vars set KEY_MODELLER=" not in script
