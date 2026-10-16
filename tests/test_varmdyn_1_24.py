"""Integration tests for the complete 1-24 VarMDyn equilibration pipeline."""

import json
from pathlib import Path

import pytest

from colabmda.openmm.modular.em import run_em
from colabmda.openmm.modular.md import run_md
from colabmda.openmm.modular.npt import run_npt
from colabmda.openmm.modular.nvt import run_nvt
from colabmda.openmm.modular.qc import analyze_logs


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_varmdyn_1_24_full_pipeline(tmp_path: Path, monkeypatch):
    monkeypatch.setenv("COLABMDA_TEST_1_24", "1")

    # Setup minimal dipeptide structure
    pdb_lines = [
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n",
        "ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C\n",
        "ATOM      3  C   ALA A   1       2.009   1.428   0.000  1.00  0.00           C\n",
        "ATOM      4  O   ALA A   1       1.272   2.408   0.000  1.00  0.00           O\n",
        "ATOM      5  CB  ALA A   1       2.009  -0.784  -1.204  1.00  0.00           C\n",
        "ATOM      6  N   ALA A   2       3.328   1.536   0.000  1.00  0.00           N\n",
        "ATOM      7  CA  ALA A   2       3.985   2.842   0.000  1.00  0.00           C\n",
        "ATOM      8  C   ALA A   2       5.498   2.704   0.000  1.00  0.00           C\n",
        "ATOM      9  O   ALA A   2       6.235   3.684   0.000  1.00  0.00           O\n",
        "ATOM     10  CB  ALA A   2       3.534   3.626   1.204  1.00  0.00           C\n",
        "ATOM     11  OXT ALA A   2       5.985   1.576   0.000  1.00  0.00           O\n",
        "TER\n",
        "END\n",
    ]
    eq_dir = tmp_path / "eq"
    eq_dir.mkdir(parents=True)
    clean_pdb = eq_dir / "ala2_cleaned.pdb"
    clean_pdb.write_text("".join(pdb_lines))

    # 1. EM
    em_ok = run_em(str(eq_dir), "ala2", padding_nm="1.4", max_iterations=20)
    assert em_ok is True
    assert (eq_dir / "system.xml").exists()
    assert (eq_dir / "solvated.pdb").exists()
    assert (eq_dir / "em.state.xml").exists()

    # 2. NVT: stages 01-18
    nvt_ok = run_nvt(str(eq_dir), "ala2", protocol="varmdyn")
    assert nvt_ok is True
    assert (eq_dir / "nvt.state.xml").exists()

    # 3. NPT: stages 19-24
    npt_ok = run_npt(str(eq_dir), "ala2", protocol="varmdyn")
    assert npt_ok is True
    assert (eq_dir / "npt.state.xml").exists()
    assert (eq_dir / "equilibrated.pdb").exists()
    assert (eq_dir / "eq.pdb").exists()

    # Verify all 24 stages were generated with metadata
    eq_artifacts = eq_dir / "equilibration"
    assert eq_artifacts.is_dir()
    for stage_idx in range(1, 25):
        s_json = eq_artifacts / f"stage_{stage_idx:02d}.json"
        s_xml = eq_artifacts / f"stage_{stage_idx:02d}.state.xml"
        assert s_json.exists(), f"Missing stage_{stage_idx:02d}.json"
        assert s_xml.exists(), f"Missing stage_{stage_idx:02d}.state.xml"
        with open(s_json) as f:
            meta = json.load(f)
        assert meta["stage"] == stage_idx
        assert "description" in meta
        assert "potential_energy_kj_mol" in meta

    # 4. QC: analyzes logs and checks 24/24 stages
    analyze_logs(str(eq_dir), strict=False, protocol="varmdyn")
    with open(eq_dir / "equilibration_qc.json") as f:
        qc_data = json.load(f)

    assert qc_data["protocol"] == "varmdyn"
    assert qc_data["varmdyn_stages_total"] == 24
    assert qc_data["varmdyn_stages_verified"] == 24
    assert qc_data["varmdyn_stages_complete"] is True

    # 5. Production replica branching
    r1_dir = tmp_path / "r1"
    r1_dir.mkdir(parents=True)
    md_ok = run_md(
        str(r1_dir),
        "ala2",
        total_ns=0.00004,
        traj_interval_ps=0.01,
        checkpoint_ps=0.02,
    )
    assert md_ok is not False
    assert (r1_dir / "prod.chk").exists()
    assert (r1_dir / "npt.state.xml").exists()


def test_varmdyn_stage_resumption(tmp_path: Path, monkeypatch):
    """Test that NPT recovers from intermediate stages (e.g. stage 23) without restarting."""
    from colabmda.openmm.modular.utils import find_highest_completed_stage

    monkeypatch.setenv("COLABMDA_TEST_1_24", "1")
    eq_dir = tmp_path / "eq"
    eq_dir.mkdir(parents=True)
    eq_artifacts = eq_dir / "equilibration"
    eq_artifacts.mkdir(parents=True)

    # Empty dir returns 0
    highest, latest = find_highest_completed_stage(str(eq_dir), start_stage=19, end_stage=24)
    assert highest == 0
    assert latest is None

    # Simulate completed stages 19-23
    for s in range(19, 24):
        (eq_artifacts / f"stage_{s:02d}.state.xml").write_text("<State/>")
        (eq_artifacts / f"stage_{s:02d}.json").write_text(json.dumps({"stage": s}))

    highest, latest = find_highest_completed_stage(str(eq_dir), start_stage=19, end_stage=24)
    assert highest == 23
    assert latest == str(eq_artifacts / "stage_23.state.xml")
