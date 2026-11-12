"""Check lossless filename migration and trajectory selection across states."""

import hashlib
import shutil
import subprocess

import numpy as np
import pytest

from colabmda.core.logs import StageLog
from colabmda.openmm.merge import find_trajectory_pair, rename_trajectory_files


@pytest.fixture
def products(tmp_path):
    import MDAnalysis as mda
    import parmed as pmd

    tleap = shutil.which("tleap")
    if not tleap:
        pytest.skip("AmberTools is required for real topology/trajectory fixtures")
    deck = tmp_path / "build.in"
    deck.write_text(
        "source leaprc.protein.ff19SB\nm = sequence { ALA LEU GLY }\n"
        "saveamberparm m system.prmtop system.inpcrd\nquit\n"
    )
    subprocess.run([tleap, "-f", str(deck)], cwd=tmp_path, check=True, capture_output=True)
    structure = pmd.load_file(str(tmp_path / "system.prmtop"))
    structure.strip("!@CA")
    structure.save(str(tmp_path / "ca.prmtop"))
    universe = mda.Universe(str(tmp_path / "system.prmtop"), str(tmp_path / "system.inpcrd"))
    for top, trajectory, atoms in [
        ("stripped.prmtop", "prod_keep_ligand.dcd", universe.atoms),
        ("ca.prmtop", "prod_ca.dcd", universe.select_atoms("name CA")),
    ]:
        if top == "stripped.prmtop":
            shutil.copy2(tmp_path / "system.prmtop", tmp_path / top)
        with mda.Writer(str(tmp_path / trajectory), n_atoms=len(atoms)) as writer:
            writer.write(atoms)
            atoms.positions += np.array([0.1, 0.2, 0.3])
            writer.write(atoms)
    (tmp_path / "prod_full.log").write_text("Time (ps),Temperature (K)\n10,300\n20,300\n")
    return tmp_path


def snapshot(folder):
    return {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in folder.iterdir() if p.is_file()
    }


@pytest.mark.parametrize("state", ["apo", "atpmg"])
def test_rename_preserves_bytes_and_pairing_and_is_repeatable(products, tmp_path, state):
    replica = tmp_path / "data" / state / "simulation/WT/md/r1"
    replica.mkdir(parents=True)
    before = snapshot(products)
    mapping = {
        "stripped.prmtop": "top.prmtop",
        "prod_keep_ligand.dcd": "traj.dcd",
        "ca.prmtop": "topCa.prmtop",
        "prod_ca.dcd": "trajCa.dcd",
        "prod_full.log": "md.log",
    }
    for old in mapping:
        shutil.copy2(products / old, replica / old)
    assert rename_trajectory_files(replica) == 2
    assert snapshot(replica) == {new: before[old] for old, new in mapping.items()}
    assert find_trajectory_pair(replica) == (replica / "top.prmtop", replica / "traj.dcd")
    assert rename_trajectory_files(replica) == 2


def test_conflicting_or_uncommitted_outputs_are_never_removed(products):
    before = snapshot(products)
    (products / "top.prmtop").write_bytes(b"different")
    with pytest.raises(ValueError, match="Conflicting"):
        rename_trajectory_files(products)
    assert all(
        hashlib.sha256((products / k).read_bytes()).hexdigest() == v for k, v in before.items()
    )
    (products / "top.prmtop").unlink()
    (products / "chunks").mkdir()
    (products / "chunks/0002.dcd.tmp").write_bytes(b"unfinished")
    with pytest.raises(ValueError, match="Uncommitted"):
        rename_trajectory_files(products)
    assert snapshot(products) == before


def test_missing_or_mismatched_pairs_fail_before_any_rename(products):
    before = snapshot(products)
    (products / "prod_ca.dcd").unlink()
    with pytest.raises(ValueError, match="Both matching"):
        rename_trajectory_files(products)
    assert (products / "prod_keep_ligand.dcd").is_file()
    assert not (products / "traj.dcd").exists()
    assert (
        hashlib.sha256((products / "stripped.prmtop").read_bytes()).hexdigest()
        == before["stripped.prmtop"]
    )


@pytest.mark.parametrize("state", ["apo", "atpmg"])
def test_view_without_run_id_loads_the_movie_instead_of_a_starting_model(
    products, monkeypatch, state
):
    import MDAnalysis as mda

    from colabmda.openmm.figures import openmm_view

    monkeypatch.delenv("COLABMDA_RUN_ID", raising=False)
    root = products / "project"
    replica = root / "data" / state / "simulation/WT/md/r1"
    prep = replica.parent.parent / "prep"
    replica.mkdir(parents=True)
    prep.mkdir()
    for old in [
        "stripped.prmtop",
        "prod_keep_ligand.dcd",
        "ca.prmtop",
        "prod_ca.dcd",
        "prod_full.log",
    ]:
        shutil.copy2(products / old, replica / old)
    rename_trajectory_files(replica)
    model = mda.Universe(str(products / "system.prmtop"), str(products / "system.inpcrd"))
    model.residues.resids = [118, 119, 120]
    model.atoms.write(str(prep / "model.pdb"))
    openmm_view(state=state, variant="WT", root=root, frame="all", no_launch=True)
    output = root / "data" / state / "analysis/WT/r1"
    assert (output / "view.pdb").is_file()
    script = (output / "visualize.pml").read_text()
    assert str(replica / "traj.dcd") in script
    assert "mset 1 -2" in script and "mplay" in script
    assert "source_kind: trajectory" in StageLog(output, "VISUALIZE_DONE").path.read_text()
    assert not list(output.glob("*.json"))


def test_status_reads_amber_atom_and_residue_counts(products):
    import MDAnalysis as mda

    from colabmda.openmm.commands import _read_topology_info

    universe = mda.Universe(str(products / "stripped.prmtop"))
    info = _read_topology_info(products / "stripped.prmtop")
    assert info["atoms"] == universe.atoms.n_atoms > 0
    assert info["protein_residues"] == 3
    assert info["residues"] == 3
    assert info["waters"] == info["ions"] == 0
