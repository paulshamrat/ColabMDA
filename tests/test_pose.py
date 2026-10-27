"""Public workflow gates; no archived coordinates or MODELLER key in fixtures."""

import json
from argparse import Namespace

import numpy as np
import pytest
from Bio.PDB import Atom, Chain, Residue

from colabmda.amber.leap import generate_tleap_script
from colabmda.core.pose import fit_donor
from colabmda.core.taxonomy import SystemSpec
from colabmda.core.workflow import model, sha


def ca_chain(points, offset=0):
    chain = Chain.Chain("A")
    for i, point in enumerate(points, 1):
        residue = Residue.Residue((" ", i + offset, " "), "ALA", " ")
        residue.add(Atom.Atom("CA", np.array(point), 0, 1, " ", " CA ", i, element="C"))
        chain.add(residue)
    return chain


def test_homologous_fit_uses_sequence_not_equal_residue_numbers():
    points = np.random.default_rng(17).normal(size=(60, 3)) * 10
    target = ca_chain(points)
    donor = ca_chain(points + [4, -2, 3], offset=100)
    donor[(" ", 120, " ")]["CA"].coord += 10
    rotation, translation, report = fit_donor(target, donor, set(range(1, 61)), [1, 60])
    assert report["retained_pairs"] == 59
    assert report["core_rmsd_angstrom"] < 1e-10
    assert report["fitted_residues"][0]["donor"] == 101
    np.testing.assert_allclose((points + [4, -2, 3]) @ rotation + translation, points, atol=1e-10)


def test_fit_rejects_insufficient_experimental_coverage():
    chain = ca_chain(np.zeros((60, 3)))
    with pytest.raises(ValueError, match="Insufficient"):
        fit_donor(chain, chain, set(range(1, 20)), [1, 60])


def test_dry_pocket_topology_adds_no_solvent_or_salt(tmp_path):
    script = generate_tleap_script(
        SystemSpec.from_tag("atpmg"),
        tmp_path / "model.pdb",
        tmp_path / "p.prmtop",
        tmp_path / "p.inpcrd",
        solvate=False,
    )
    assert "solvateBox" not in script
    assert "addIonsRand" not in script
    assert "ionslm_126_opc" in script
    assert "ATP-B3.prepi" in script


def test_model_wrapper_checks_inputs_and_reuses_only_matching_request(tmp_path, monkeypatch):
    from colabmda.core import workflow

    folder = tmp_path / "data/apo/structure"
    (folder / "mapping").mkdir(parents=True)
    (folder / "raw").mkdir()
    (folder / "STRUCTURE_DONE").write_text("Verified public inputs\n")
    for name in ["template.pdb", "uniprot.fasta"]:
        (folder / "raw" / name).write_text("public-input fixture\n")
    inputs = dict(
        pdb="4BGQ",
        uniprot="O76039",
        chain="A",
        range=[1, 1],
        sequence="A",
        files={"raw/template.pdb": sha(folder / "raw/template.pdb")},
    )
    (folder / "mapping/input.json").write_text(json.dumps(inputs))
    calls = []

    def build(argv, **kwargs):
        from pathlib import Path

        calls.append(argv)
        out_dir = Path(argv[argv.index("--outdir") + 1])
        out_dir.mkdir(parents=True, exist_ok=True)
        out = out_dir / "WT.pdb"
        out.write_text(
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n"
        )

    monkeypatch.setattr("subprocess.run", build)
    args = Namespace(root=tmp_path, state="apo", variant="WT", seed=13031)
    assert workflow.model(args) == 0
    assert model(args) == 0
    assert len(calls) == 1
    args.seed = 13032
    with pytest.raises(ValueError, match="different inputs/settings"):
        model(args)
    (folder / "raw/template.pdb").write_text("changed input\n")
    with pytest.raises(ValueError, match="Changed structure input"):
        model(args)


def test_missing_state_mutant_does_not_fall_back_to_wt(tmp_path, monkeypatch):
    from colabmda.cli import _dispatch_prep

    monkeypatch.setenv("COLABMDA_RUN_ID", "public01")
    folder = tmp_path / "data/apo/structure/models/WT"
    folder.mkdir(parents=True)
    (folder / "WT.pdb").write_text("WT must not be used for L119R\n")
    args = Namespace(
        root=tmp_path,
        state="apo",
        protein="cdkl5",
        variant="L119R",
        receptor=None,
        ligand=None,
        ions=None,
        phos=None,
        peptide_seq=None,
        peptide_name=None,
    )
    with pytest.raises(FileNotFoundError, match="L119R"):
        _dispatch_prep(args)


def test_invalid_public_variant_rejected_before_filesystem(tmp_path):
    with pytest.raises(ValueError, match="Invalid variant"):
        model(Namespace(root=tmp_path, state="apo", variant="../WT", seed=1))


def test_model_refuses_mutation_without_experimental_coordinate_density(tmp_path):
    folder = tmp_path / "data/apo/structure/mapping"
    folder.mkdir(parents=True)
    (folder.parent / "STRUCTURE_DONE").write_text("Verified public inputs\n")
    (folder / "input.json").write_text(json.dumps(dict(files={})))
    (folder / "residues.json").write_text(
        json.dumps([dict(uniprot=119, resolved=False, amino_acid="L")])
    )
    with pytest.raises(ValueError, match="density gate"):
        model(Namespace(root=tmp_path, state="apo", variant="L119R", seed=13031))


def test_model_requires_completed_fetch(tmp_path):
    with pytest.raises(ValueError, match="STRUCTURE_DONE"):
        model(Namespace(root=tmp_path, state="apo", variant="WT", seed=13031))
