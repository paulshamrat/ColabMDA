"""Native-reference display checks; fixtures are synthetic, not archived proteins."""

import numpy as np
import pytest

from colabmda.openmm.figures import CDKL5_CAMERAS, _view_dataset


def protein(path, shift=0, mutation=False):
    rows = []
    for i in range(100, 136):
        x, y, z = np.sin(i) * 4, np.cos(i) * 4, (i - 100) * 0.4
        name = "ARG" if mutation and i == 119 else ("LEU" if i == 119 else "ALA")
        rows.append(
            f"ATOM  {i:5d}  CA  {name} A{i:4d}    {x + shift:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n"
        )
    path.write_text("".join(rows) + "END\n")
    return path


def test_native_labels_fit_and_no_invented_zero_time(tmp_path):
    ref = protein(tmp_path / "ref.pdb")
    moving = protein(tmp_path / "moving.pdb", shift=100)
    u, panels, record = _view_dataset(moving, reference=ref)
    assert len(panels) == 1
    assert panels[0]["time_ps"] is None
    assert panels[0]["core_rmsd_angstrom"] < 1e-4
    assert record["source_kind"] == "structure"
    assert record["native_residues"][19] == ("A", 119, "")
    assert u.select_atoms("protein and name CA").resids.tolist() == list(range(100, 136))


def test_declared_variant_identity_gate(tmp_path):
    ref = protein(tmp_path / "ref.pdb")
    mutant = protein(tmp_path / "mut.pdb", mutation=True)
    with pytest.raises(ValueError, match="Declared variant is absent"):
        _view_dataset(ref, reference=ref, variant="L119R")
    with pytest.raises(ValueError, match="sequence mismatch"):
        _view_dataset(mutant, reference=ref)
    _, panels, _ = _view_dataset(mutant, reference=ref, variant="L119R")
    assert panels[0]["core_rmsd_angstrom"] < 1e-4
    with pytest.raises(ValueError, match="Declared variant"):
        _view_dataset(mutant, reference=ref, variant="L118R")


def test_invalid_frame_and_absent_native_highlight(tmp_path):
    ref = protein(tmp_path / "ref.pdb")
    with pytest.raises(ValueError, match="outside trajectory"):
        _view_dataset(ref, frames="2")
    with pytest.raises(ValueError, match="absent"):
        _view_dataset(ref, resi=12)


def test_activation_camera_rotation_is_proper():
    # Stored PyMOL views carry about 2e-5 accumulated rotation drift.
    # Preserve the activation values rather than silently replacing its view.
    for view in CDKL5_CAMERAS.values():
        rotation = np.array(view[:9]).reshape(3, 3)
        np.testing.assert_allclose(rotation @ rotation.T, np.eye(3), atol=5e-5)
        assert np.linalg.det(rotation) == pytest.approx(1, abs=5e-5)


def test_log_times_are_used_and_mismatched_counts_fail(tmp_path):
    ref = protein(tmp_path / "ref.pdb")
    import MDAnalysis as mda

    u = mda.Universe(str(ref))
    dcd = tmp_path / "sample.dcd"
    with mda.Writer(str(dcd), n_atoms=len(u.atoms)) as writer:
        writer.write(u.atoms)
        u.atoms.positions += 50
        writer.write(u.atoms)
    log = tmp_path / "prod_full.log"
    log.write_text("Time (ps),Temperature (K)\n10,300\n100,300\n")
    _, panels, _ = _view_dataset(ref, dcd, ref, log=log)
    assert [p["time_ps"] for p in panels] == [10, 100]
    assert all(p["core_rmsd_angstrom"] < 1e-4 for p in panels)
    log.write_text("Time (ps)\n10\n")
    with pytest.raises(ValueError, match="times do not match"):
        _view_dataset(ref, dcd, ref, log=log)


def test_movie_script_keeps_native_labels_and_common_reference(tmp_path):
    import json

    from colabmda.openmm.figures import openmm_view

    ref = protein(tmp_path / "ref.pdb")
    import MDAnalysis as mda

    u = mda.Universe(str(ref))
    dcd = tmp_path / "sample.dcd"
    with mda.Writer(str(dcd), n_atoms=len(u.atoms)) as writer:
        writer.write(u.atoms)
    openmm_view(
        pdb_dir=tmp_path,
        topology=ref,
        trajectory=dcd,
        reference=ref,
        frame="all",
        no_launch=True,
    )
    script = (tmp_path / "visualize.pml").read_text()
    assert "load_traj" in script
    assert "zoom reference and polymer.protein" in script
    assert "resi 119" in script
    assert "Switch" not in script
    assert json.loads((tmp_path / "visualize.json").read_text())["display_frames"] == "all"


def test_custom_snapshot_times_require_a_production_log(tmp_path):
    import json

    from colabmda.openmm.figures import openmm_snapshots

    config = tmp_path / "config.json"
    config.write_text(json.dumps({"systems": {"old": {"times": [0, 1]}}}))
    with pytest.raises(ValueError, match="matching production log"):
        openmm_snapshots(config_path=config)
