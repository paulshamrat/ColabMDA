"""Numbering, archive fallback and mutation gates for natural-protein tutorials."""

import io
import json
from urllib.error import HTTPError

import pytest

from colabmda.core.examples import install_example
from colabmda.core.taxonomy import SystemSpec
from colabmda.modeller.structure import fetch_archive, verify_sequence
from colabmda.modeller.substrate import assemble_system_model


def coordinate(residue="GLY", position=2, chain="A"):
    return f"ATOM      1  CA  {residue} {chain}{position:4d}       1.000   2.000   3.000  1.00  0.00           C\nEND\n"


def test_examples_coexist_without_replacing_trpcage(tmp_path):
    for name in ["trpcage", "fkbp", "hras"]:
        install_example(name, tmp_path)
    original = (tmp_path / "data/00_setup/example.json").read_bytes()
    for name in ["fkbp", "hras"]:
        receipt = json.loads((tmp_path / f"data/00_setup/{name}.json").read_text())
        assert len(receipt["files"]) == 7
        install_example(name, tmp_path)
        script = tmp_path / f"scripts/01_structure/{name}/02_model.py"
        script.write_text("user edit")
        with pytest.raises(ValueError, match="preserved"):
            install_example(name, tmp_path)
    assert (tmp_path / "data/00_setup/example.json").read_bytes() == original


def test_mmcif_failure_falls_back_and_retains_archive(tmp_path, monkeypatch):
    urls = []

    def get(url, timeout):
        urls.append(url)
        if url.endswith(".cif"):
            raise HTTPError(url, 404, "missing", {}, None)
        return io.BytesIO(coordinate().encode())

    monkeypatch.setattr("colabmda.modeller.structure.urlopen", get)
    structure, path, url = fetch_archive("test", tmp_path)
    assert urls[0].endswith(".cif") and url.endswith(".pdb")
    assert path.read_text() == coordinate()
    assert structure[0]["A"][2].resname == "GLY"
    path.write_text("user data")
    with pytest.raises(ValueError, match="preserved"):
        fetch_archive("test", tmp_path)


def test_sequence_checks_fail_on_identity_numbering_chain_and_density(tmp_path):
    path = tmp_path / "model.pdb"
    path.write_text(coordinate())
    verify_sequence(path, "A", "G", 2)
    for chain, sequence, start in [("B", "G", 2), ("A", "V", 2), ("A", "G", 1), ("A", "GG", 2)]:
        with pytest.raises(ValueError):
            verify_sequence(path, chain, sequence, start)
    path.write_text(coordinate().replace(" CA ", " N  "))
    with pytest.raises(ValueError, match="C-alpha"):
        verify_sequence(path, "A", "G", 2)


def test_prebuilt_mutant_is_preserved_and_wrong_residue_fails(tmp_path):
    path, out = tmp_path / "mutant.pdb", tmp_path / "assembled.pdb"
    path.write_text(coordinate("VAL", 12))
    spec = SystemSpec.from_tag("apo", variant="G12V")
    assemble_system_model(spec, path, out)
    assert "VAL" in out.read_text()
    path.write_text(coordinate("ALA", 12))
    with pytest.raises(ValueError, match="identity mismatch"):
        assemble_system_model(spec, path, out)


def test_modeling_cli_reports_failure_for_rejected_mutation(tmp_path, monkeypatch):
    from colabmda.modeller import engine

    path = tmp_path / "wt.pdb"
    path.write_text(coordinate("GLY", 12))
    monkeypatch.setattr(engine, "_MODELLER_OK", True)
    monkeypatch.setattr(engine, "mutate_once", lambda *a, **k: (None, None, None, "G", "FAIL"))
    monkeypatch.setattr("sys.argv", ["model", "--pdb-in", str(path), "--mut", "G12V"])
    with pytest.raises(SystemExit, match="mutations failed"):
        engine.main()
