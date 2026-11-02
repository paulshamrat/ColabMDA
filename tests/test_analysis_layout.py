"""Analysis separation, replica isolation and verified migration regressions."""

import hashlib
import json
from argparse import Namespace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from colabmda.core.project import ProjectPaths
from colabmda.core.taxonomy import SystemSpec
from colabmda.openmm.analysis import aggregate_reports


def completed_legacy(tmp_path, monkeypatch):
    monkeypatch.setenv("COLABMDA_RUN_ID", "old")
    paths = ProjectPaths(tmp_path, state="apo")
    spec = SystemSpec()
    source = paths.variant_dir(spec) / "04_analysis"
    for name in [
        "ANALYSIS_DONE",
        "analysis.json",
        "rmsd.csv",
        "rg.csv",
        "rmsf.csv",
        "figures/rmsd_vs_time.png",
        "figures/rg_vs_time.png",
        "figures/rmsf_per_residue.png",
    ]:
        p = source / name
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(("original:" + name).encode())
    return paths, spec, source


def test_completed_analysis_migration_preserves_every_byte_and_receipt(tmp_path, monkeypatch):
    paths, spec, source = completed_legacy(tmp_path, monkeypatch)
    before = {
        str(p.relative_to(source)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in source.rglob("*")
        if p.is_file()
    }
    result = paths.migrate_analysis(spec)
    assert result["files"] == before
    assert result["destination"] == "data/apo/analysis/old/WT/r1"
    assert not source.exists()
    target = paths.replica_analysis_dir(spec)
    assert {
        str(p.relative_to(target)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in target.rglob("*")
        if p.is_file()
    } == before
    receipt = tmp_path / "data/00_setup/analysis_migrations/apo_old_WT.json"
    assert json.loads(receipt.read_text()) == result
    assert paths.migrate_analysis(spec) is None


def test_migration_refuses_conflicts_incomplete_reports_and_symlinks(tmp_path, monkeypatch):
    paths, spec, source = completed_legacy(tmp_path, monkeypatch)
    target = paths.replica_analysis_dir(spec)
    target.mkdir(parents=True)
    (target / "rmsd.csv").write_text("different")
    with pytest.raises(ValueError, match="Conflicting"):
        paths.migrate_analysis(spec)
    assert (source / "ANALYSIS_DONE").is_file()
    (target / "rmsd.csv").unlink()
    (source / "extra").symlink_to(source / "rmsd.csv")
    with pytest.raises(ValueError, match="symlink"):
        paths.migrate_analysis(spec)
    (source / "extra").unlink()
    (source / "ANALYSIS_DONE").unlink()
    with pytest.raises(ValueError, match="completed"):
        paths.migrate_analysis(spec)


def test_migration_changed_source_keeps_original_without_new_completion(tmp_path, monkeypatch):
    import shutil

    paths, spec, source = completed_legacy(tmp_path, monkeypatch)
    original = shutil.copy2

    def changing_copy(src, dst):
        result = original(src, dst)
        if Path(src).name == "rmsd.csv":
            Path(src).write_text("changed while copying")
        return result

    monkeypatch.setattr(shutil, "copy2", changing_copy)
    with pytest.raises(ValueError, match="changed during migration"):
        paths.migrate_analysis(spec)
    assert source.is_dir()
    assert not (paths.replica_analysis_dir(spec) / "ANALYSIS_DONE").exists()


def test_aggregation_rejects_unmatched_times_and_native_labels(tmp_path):
    folders = [tmp_path / "r1", tmp_path / "r2"]
    for i, folder in enumerate(folders):
        folder.mkdir()
        for metric in ["rmsd", "rg"]:
            pd.DataFrame({"Time (ps)": [10, 20], "Value": [i + 1, i + 3]}).to_csv(
                folder / f"{metric}.csv", index=False
            )
        pd.DataFrame({"Residue": ["118", "119A"], "Value": [i + 1, i + 3]}).to_csv(
            folder / "rmsf.csv", index=False
        )
    out = tmp_path / "aggregate"
    aggregate_reports(folders, out)
    values = pd.read_csv(out / "rmsd_avg.csv")
    np.testing.assert_allclose(values["Mean"], [1.5, 3.5])
    np.testing.assert_allclose(values["Std"], [0.5, 0.5])
    frame = pd.read_csv(folders[1] / "rmsd.csv")
    frame.iloc[1, 0] = 30
    frame.to_csv(folders[1] / "rmsd.csv", index=False)
    with pytest.raises(ValueError, match="axes differ"):
        aggregate_reports(folders, out)
    frame.iloc[1, 0] = 20
    frame.to_csv(folders[1] / "rmsd.csv", index=False)
    frame = pd.read_csv(folders[1] / "rmsf.csv")
    frame.iloc[1, 0] = "119B"
    frame.to_csv(folders[1] / "rmsf.csv", index=False)
    with pytest.raises(ValueError, match="axes differ"):
        aggregate_reports(folders, out)
    assert not (out / "AGGREGATE_DONE").exists()


def test_real_amber_replica_analysis_routes_reports_outside_simulation(tmp_path, monkeypatch):
    import shutil
    import subprocess

    import MDAnalysis as mda
    from Bio.PDB import PDBIO, PDBParser

    from colabmda.cli import _dispatch_analyze

    tleap = shutil.which("tleap")
    if not tleap:
        pytest.skip("AmberTools required for real replica routing acceptance")
    monkeypatch.setenv("COLABMDA_RUN_ID", "routing")
    paths = ProjectPaths(tmp_path, state="apo")
    spec = SystemSpec()
    prep = paths.prep_dir(spec)
    prep.mkdir(parents=True)
    deck = prep / "build.in"
    deck.write_text(
        "source leaprc.protein.ff19SB\nm = sequence { ALA LEU GLY }\nsaveamberparm m tiny.prmtop tiny.rst7\nsavepdb m raw.pdb\nquit\n"
    )
    subprocess.run([tleap, "-f", str(deck)], cwd=prep, check=True, capture_output=True)
    model = PDBParser(QUIET=True).get_structure("model", prep / "raw.pdb")
    for i, residue in enumerate(model.get_residues(), 118):
        residue.id = (" ", i, " ")
    writer = PDBIO()
    writer.set_structure(model)
    writer.save(str(prep / "model.pdb"))
    u = mda.Universe(str(prep / "tiny.prmtop"), str(prep / "tiny.rst7"), format="INPCRD")
    xyz = u.atoms.positions.copy()
    for rep, seed in [("r1", 1), ("r2", 2)]:
        folder = paths.replica_dir(spec, rep)
        folder.mkdir(parents=True)
        shutil.copy2(prep / "tiny.prmtop", folder / "stripped.prmtop")
        rng = np.random.default_rng(seed)
        with mda.Writer(str(folder / "prod_keep_ligand.dcd"), n_atoms=len(u.atoms)) as out:
            for _ in range(4):
                u.atoms.positions = xyz + rng.normal(0, 0.05, xyz.shape)
                out.write(u.atoms)
    args = Namespace(root=str(tmp_path), state="apo", variant="WT", replicas="r1,r2")
    assert _dispatch_analyze(args) == 0
    for rep in ["r1", "r2"]:
        folder = paths.replica_analysis_dir(spec, rep)
        assert (folder / "ANALYSIS_DONE").is_file()
        assert pd.read_csv(folder / "rmsf.csv")["Residue"].tolist() == [118, 119, 120]
        assert pd.read_csv(folder / "rmsd.csv")["Time (ps)"].tolist() == [10, 20, 30, 40]
    assert (paths.analysis_dir(spec) / "aggregate/AGGREGATE_DONE").is_file()
    assert not list(paths.variant_dir(spec).rglob("rmsd.csv"))
    # Failure must clear completion and invalidate a previously computed aggregate.
    monkeypatch.setattr(
        "colabmda.openmm.commands.openmm_analysis",
        lambda **kw: (_ for _ in ()).throw(RuntimeError("interrupted")),
    )
    args.replicas = "r1"
    with pytest.raises(RuntimeError, match="interrupted"):
        _dispatch_analyze(args)
    assert not (paths.replica_analysis_dir(spec) / "ANALYSIS_DONE").exists()
    assert not (paths.analysis_dir(spec) / "aggregate/AGGREGATE_DONE").exists()
