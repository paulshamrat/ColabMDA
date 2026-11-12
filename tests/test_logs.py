"""Completion must be verified through stage logs, without standalone flags."""

import pytest

from colabmda.core.logs import StageLog, migrate_project_logs


def preparation(root):
    folder = root / "data/apo/simulation/WT/prep"
    folder.mkdir(parents=True)
    for name in ["system.prmtop", "system.inpcrd", "solvated.pdb"]:
        (folder / name).write_text("synthetic output: " + name)
    return folder


def test_completion_log_mirrors_data_and_detects_changed_or_missing_outputs(tmp_path):
    folder = preparation(tmp_path)
    record = StageLog(folder, "PREP_DONE")
    record.write_text("Preparation checks passed\n")
    assert record.path == tmp_path / "logs/apo/simulation/WT/prep/prep.log"
    assert record.exists()
    assert not list((tmp_path / "data").rglob("*_DONE"))
    (folder / "system.inpcrd").write_text("different coordinates")
    assert not record.exists()
    record.write_text("Reverified after intentional rebuild\n")
    assert record.exists()
    (folder / "solvated.pdb").unlink()
    assert not record.exists()
    with pytest.raises(ValueError, match="required outputs"):
        record.write_text("Must not claim success\n")
    assert not record.exists()


def test_invalidation_preserves_legacy_metadata_and_blocks_legacy_fallback(tmp_path):
    folder = preparation(tmp_path)
    legacy = folder / "PREP_DONE"
    legacy.write_text("Important original metadata\n")
    record = StageLog(folder, "PREP_DONE")
    assert record.exists()
    record.unlink(missing_ok=True)
    assert not legacy.exists()
    assert "Important original metadata" in record.path.read_text()
    assert not record.exists()
    # A stale flag cannot override a failed/restarted stage's current log.
    legacy.write_text("stale")
    assert not record.exists()
    with pytest.raises(ValueError, match="Current stage log is incomplete"):
        migrate_project_logs(tmp_path)
    assert legacy.exists()


def test_project_migration_preserves_scientific_bytes_and_metadata(tmp_path):
    folder = preparation(tmp_path)
    legacy = folder / "PREP_DONE"
    details = "Legacy preparation verified\nseed: 42\n"
    legacy.write_text(details)
    before = {p.name: p.read_bytes() for p in folder.iterdir() if p.name != "PREP_DONE"}
    paths = migrate_project_logs(tmp_path)
    assert paths == ["logs/apo/simulation/WT/prep/prep.log"]
    assert not list((tmp_path / "data").rglob("*_DONE"))
    assert {p.name: p.read_bytes() for p in folder.iterdir()} == before
    record = StageLog(folder, "PREP_DONE")
    assert record.exists() and details.strip() in record.read_text()
    assert migrate_project_logs(tmp_path) == []


def test_migration_refuses_incomplete_runs_before_removing_any_flags(tmp_path):
    folder = preparation(tmp_path)
    legacy = folder / "PREP_DONE"
    legacy.write_text("verified")
    (folder / "solvated.pdb").unlink()
    with pytest.raises(ValueError, match="Incomplete outputs"):
        migrate_project_logs(tmp_path)
    assert legacy.exists()
    chunks = folder.parent / "md/r1/chunks"
    chunks.mkdir(parents=True)
    (chunks / "0002.dcd.tmp").write_text("unfinished")
    with pytest.raises(ValueError, match="Unfinished production"):
        migrate_project_logs(tmp_path)
    assert legacy.exists()


def test_completion_log_refuses_an_artifact_outside_its_data_folder(tmp_path):
    folder = preparation(tmp_path)
    record = StageLog(folder, "PREP_DONE")
    record.write_text("verified")
    text = record.path.read_text()
    record.path.write_text(text.replace("\tsystem.prmtop\n", "\t../../secret\n"))
    assert not record.exists()


def test_colab_push_keeps_completion_logs_with_scientific_inputs(tmp_path, monkeypatch):
    import tarfile

    from colabmda.remote import colab

    folder = preparation(tmp_path)
    StageLog(folder, "PREP_DONE").write_text("verified inputs")
    members = []

    def upload(command):
        assert command[:2] == ["colab", "upload"]
        with tarfile.open(command[-2]) as archive:
            members.extend(archive.getnames())
        return 0

    monkeypatch.setattr(colab.subprocess, "call", upload)
    monkeypatch.setattr(colab, "_execute_remote", lambda *a, **k: 0)
    assert colab.colab_push(root=tmp_path) == 0
    assert "data/apo/simulation/WT/prep/system.prmtop" in members
    assert "logs/apo/simulation/WT/prep/prep.log" in members
    assert not any(name.endswith("_DONE") for name in members)
