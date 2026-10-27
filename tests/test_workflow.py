"""Guards for direct scientific workflows; fixtures contain no archived structures."""

import json
from argparse import Namespace

import numpy as np
import pytest
from Bio.PDB import PDBParser
from openmm import System, unit
from openmm.app import Topology, element

from colabmda.core.project import ProjectPaths
from colabmda.core.taxonomy import SystemSpec
from colabmda.core.workflow import check, sha
from colabmda.openmm.common import add_solute_heavy_atom_restraints, find_highest_completed_stage
from colabmda.openmm.qc import analyze_logs


def test_state_layout_and_existing_run_compatibility(tmp_path, monkeypatch):
    monkeypatch.setenv("COLABMDA_RUN_ID", "trial01")
    paths = ProjectPaths(tmp_path, state="atp_1mg")
    spec = SystemSpec.from_tag("atp_1mg", variant="L119R")
    assert paths.structure_dir == tmp_path / "data/atpmg/structure"
    assert paths.prep_dir(spec) == tmp_path / "data/atpmg/simulation/trial01/L119R/01_prep"
    assert paths.analysis_dir(spec) == tmp_path / "data/atpmg/analysis/trial01/L119R"
    assert paths.replica_analysis_dir(spec, "2").name == "r2"
    assert paths.compare_dir(spec.state_tag) == tmp_path / "data/atpmg/analysis/trial01/compare"
    assert paths.log_dir(spec) == tmp_path / "logs/atpmg/trial01/L119R"
    (tmp_path / "data/02_simulation/trial01").mkdir(parents=True)
    legacy = ProjectPaths(tmp_path, state="apo")
    assert legacy.state_dir is None
    assert legacy.prep_dir(SystemSpec()) == tmp_path / "data/02_simulation/trial01/apo/WT/01_prep"


def test_fresh_state_uses_flat_simulation_tree_without_run_id(tmp_path, monkeypatch):
    monkeypatch.delenv("COLABMDA_RUN_ID", raising=False)
    paths = ProjectPaths(tmp_path, state="apo")
    assert paths.prep_dir(SystemSpec()) == tmp_path / "data/apo/simulation/WT/01_prep"


def test_restraints_exclude_amber_salt_but_include_cofactor():
    topology = Topology()
    chain = topology.addChain("A")
    system = System()
    for name, symbol in [
        ("ALA", "C"),
        ("WAT", "O"),
        ("Na+", "Na"),
        ("Cl-", "Cl"),
        ("K+", "K"),
        ("MG", "Mg"),
    ]:
        residue = topology.addResidue(name, chain)
        topology.addAtom(symbol, element.get_by_symbol(symbol), residue)
        system.addParticle(1)
    _, count = add_solute_heavy_atom_restraints(
        system, topology, np.zeros((6, 3)) * unit.nanometer, "k"
    )
    force = system.getForce(0)
    assert count == 2
    assert [force.getParticleParameters(i)[0] for i in range(2)] == [0, 5]


def test_recovery_rejects_a_gap(tmp_path):
    folder = tmp_path / "equilibration"
    folder.mkdir()
    for stage in [19, 21]:
        (folder / f"stage_{stage:02d}.state.xml").write_text("<State/>")
        (folder / f"stage_{stage:02d}.json").write_text(json.dumps({"stage": stage}))
    with pytest.raises(RuntimeError, match="Non-contiguous"):
        find_highest_completed_stage(str(tmp_path), 19, 24)


def test_model_check_removes_marker_on_failed_identity(tmp_path):
    structure = tmp_path / "data/apo/structure"
    folder = structure / "models/WT"
    folder.mkdir(parents=True)
    mapping = structure / "mapping"
    mapping.mkdir()
    (mapping / "input.json").write_text(
        json.dumps({"files": {}, "sequence": "A", "range": [1, 1], "chain": "A"})
    )
    (folder / "WT.pdb").write_text(
        "ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n"
    )
    (folder / "MODEL_DONE").touch()
    with pytest.raises(ValueError, match="amino acids"):
        check(Namespace(root=str(tmp_path), state="apo", variant="WT", model=None))
    assert not (folder / "MODEL_DONE").exists()


def test_research_rejects_short_equilibration(tmp_path, monkeypatch):
    from colabmda.cli import _dispatch_sim

    monkeypatch.setenv("COLABMDA_RUN_ID", "trial")
    args = Namespace(
        root=str(tmp_path),
        state="apo",
        variant="WT",
        stage="md",
        protocol="quick",
        software_test=False,
    )
    with pytest.raises(SystemExit, match="Quick equilibration requires"):
        _dispatch_sim(args)
    args.protocol = "staged"
    with pytest.raises(SystemExit, match="PREP_DONE"):
        _dispatch_sim(args)
    monkeypatch.setattr("colabmda.core.workflow.verify_prepared", lambda *a: None)
    with pytest.raises(SystemExit, match="completed staged equilibration"):
        _dispatch_sim(args)
    monkeypatch.setenv("COLABMDA_TEST_1_24", "1")
    with pytest.raises(SystemExit, match="Reduced 1-24"):
        _dispatch_sim(args)


def test_strict_qc_requires_real_durations(tmp_path):
    (tmp_path / "nvt.log").write_text("Time (ps),Temperature (K)\n360,300\n460,300\n")
    (tmp_path / "npt.log").write_text("Time (ps),Density (g/mL)\n1940,1\n1960,1\n")
    folder = tmp_path / "equilibration"
    folder.mkdir()
    times = [0] * 3 + list(range(10, 61, 10)) + [60] * 5 + list(range(160, 961, 100)) + [1960]
    assert len(times) == 24
    for stage, time in enumerate(times, 1):
        (folder / f"stage_{stage:02d}.json").write_text(
            json.dumps(
                {
                    "stage": stage,
                    "time_ps": time,
                    "potential_energy_kj_mol": -10,
                    "test_mode": False,
                }
            )
        )
        (folder / f"stage_{stage:02d}.state.xml").write_text(f'<State time="{time}"/>')
    assert analyze_logs(str(tmp_path), strict=True, protocol="varmdyn")
    (folder / "stage_24.state.xml").write_text('<State time="0.2"/>')
    assert not analyze_logs(str(tmp_path), strict=True, protocol="varmdyn")
    record = json.loads((tmp_path / "equilibration_qc.json").read_text())
    assert record["varmdyn_invalid_timing_stages"] == [24]


def test_noncovalent_components_are_terminated(tmp_path):
    from colabmda.modeller.substrate import terminate_components

    path = tmp_path / "model.pdb"
    path.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\nHETATM    2  P   ATP A   2       2.000   0.000   0.000  1.00  0.00           P\nHETATM    3 MG    MG A   3       4.000   0.000   0.000  1.00  0.00          MG\nEND\n"
    )
    terminate_components(path)
    assert path.read_text().count("TER\n") == 2
    digest = sha(path)
    terminate_components(path)
    assert sha(path) == digest
    assert len(list(PDBParser(QUIET=True).get_structure("a", path).get_atoms())) == 3


def test_fetch_maps_label_offset_and_preserves_insertion(tmp_path, monkeypatch):
    import io

    from colabmda.core import workflow

    raw = tmp_path / "data/apo/structure/raw"
    raw.mkdir(parents=True)
    archive = raw / "test.cif"
    archive.write_text("synthetic archive")
    pdb = "ATOM      1  CA  ALA A  42A      0.000   0.000   0.000  1.00  0.00           C\nEND\n"
    structure = PDBParser(QUIET=True).get_structure("test", io.StringIO(pdb))
    monkeypatch.setattr(
        workflow,
        "fetch_archive",
        lambda *a: (structure, archive, "https://example.invalid/test.cif"),
    )
    monkeypatch.setattr(workflow, "urlopen", lambda *a, **k: io.BytesIO(b">synthetic\nGA\n"))
    monkeypatch.setattr(
        workflow,
        "MMCIF2Dict",
        lambda *a: {
            "_struct_ref_seq.pdbx_strand_id": ["A"],
            "_struct_ref_seq.pdbx_db_accession": ["TEST"],
            "_struct_ref_seq.seq_align_beg": ["2"],
            "_struct_ref_seq.seq_align_end": ["3"],
            "_struct_ref_seq.db_align_beg": ["1"],
            "_atom_site.auth_asym_id": ["A"],
            "_atom_site.label_atom_id": ["CA"],
            "_atom_site.auth_seq_id": ["42"],
            "_atom_site.pdbx_PDB_ins_code": ["A"],
            "_atom_site.label_seq_id": ["3"],
        },
    )
    args = Namespace(
        root=str(tmp_path),
        state="apo",
        pdb_id="TEST",
        uniprot_id="TEST",
        chain="A",
        range=[1, 2],
        mutations="A2V",
    )
    assert workflow.fetch(args) == 0
    res_path = (raw.parent / "residues.json") if (raw.parent / "residues.json").exists() else (raw.parent / "mapping/residues.json")
    rows = json.loads(res_path.read_text())
    assert rows[0]["resolved"] is False
    assert rows[1] == {
        "uniprot": 2,
        "amino_acid": "A",
        "resolved": True,
        "author": 42,
        "insertion_code": "A",
    }
    tmpl_path = (raw.parent / "template.pdb") if (raw.parent / "template.pdb").exists() else (raw / "template.pdb")
    assert "  42A" in tmpl_path.read_text()
    args.mutations = "G1A"
    with pytest.raises(ValueError, match="coordinate gate"):
        workflow.fetch(args)
    assert not (raw.parent / "STRUCTURE_DONE").exists()


def test_preparation_cannot_reuse_changed_inputs(tmp_path, monkeypatch):
    from colabmda.core.workflow import preparation_request

    monkeypatch.setenv("COLABMDA_RUN_ID", "trial")
    paths = ProjectPaths(tmp_path, state="apo")
    receptor = tmp_path / "model.pdb"
    receptor.write_text("synthetic input")
    assert preparation_request(paths, SystemSpec(), receptor, None) is False
    receptor.write_text("changed input")
    with pytest.raises(ValueError, match="fresh COLABMDA_RUN_ID"):
        preparation_request(paths, SystemSpec(), receptor, None)


def test_fetch_creates_flat_directory_structure(tmp_path, monkeypatch):
    import io
    from colabmda.core import workflow

    archive_content = "synthetic archive\n"
    pdb = "ATOM      1  CA  ALA A  42A      0.000   0.000   0.000  1.00  0.00           C\nEND\n"
    structure = PDBParser(QUIET=True).get_structure("test", io.StringIO(pdb))

    def mock_fetch(pdb_id, dest):
        f = dest / "test.cif"
        f.write_text(archive_content)
        return structure, f, "https://example.invalid/test.cif"

    monkeypatch.setattr(workflow, "fetch_archive", mock_fetch)
    monkeypatch.setattr(workflow, "urlopen", lambda *a, **k: io.BytesIO(b">synthetic\nGA\n"))
    monkeypatch.setattr(
        workflow,
        "MMCIF2Dict",
        lambda *a: {
            "_struct_ref_seq.pdbx_strand_id": ["A"],
            "_struct_ref_seq.pdbx_db_accession": ["TEST"],
            "_struct_ref_seq.seq_align_beg": ["2"],
            "_struct_ref_seq.seq_align_end": ["3"],
            "_struct_ref_seq.db_align_beg": ["1"],
            "_atom_site.auth_asym_id": ["A"],
            "_atom_site.label_atom_id": ["CA"],
            "_atom_site.auth_seq_id": ["42"],
            "_atom_site.pdbx_PDB_ins_code": ["A"],
            "_atom_site.label_seq_id": ["3"],
        },
    )
    args = Namespace(
        root=str(tmp_path),
        state="apo",
        pdb_id="TEST",
        uniprot_id="TEST",
        chain="A",
        range=[1, 2],
        mutations="A2V",
    )
    assert workflow.fetch(args) == 0
    folder = tmp_path / "data/apo/structure/TEST_TEST"
    assert (folder / "STRUCTURE_DONE").exists()
    assert (folder / "input.json").exists()
    assert (folder / "residues.json").exists()
    assert (folder / "template.pdb").exists()
    assert (folder / "uniprot.fasta").exists()
    assert not (folder / "raw").exists()
    assert not (folder / "mapping").exists()
