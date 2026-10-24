"""Quality Control and stability verification unit tests."""

import json
from pathlib import Path

import pandas as pd

from colabmda.openmm.figures import plot_equilibration
from colabmda.openmm.qc import analyze_logs, parse_em_log


def test_parse_em_log(tmp_path: Path):
    # Missing em.chk
    assert parse_em_log(str(tmp_path)) == "Failed"

    # Present em.chk
    (tmp_path / "em.chk").touch()
    assert parse_em_log(str(tmp_path)) == "Minimized"


def test_analyze_logs_passing(tmp_path: Path):
    # Setup mock NVT log
    df_nvt = pd.DataFrame(
        {
            "Time (ps)": [0.0, 10.0, 20.0, 30.0, 40.0, 50.0],
            "Temperature (K)": [298.0, 299.5, 300.2, 300.0, 299.8, 300.1],
            "Potential Energy (kJ/mole)": [-1000.0, -1500.0, -1800.0, -1850.0, -1860.0, -1855.0],
        }
    )
    df_nvt.to_csv(tmp_path / "nvt.log", index=False)

    # Setup mock NPT log
    df_npt = pd.DataFrame(
        {
            "Time (ps)": [0.0, 10.0, 20.0, 30.0, 40.0, 50.0],
            "Density (g/mL)": [0.995, 0.998, 1.001, 1.002, 1.002, 1.002],
            "Potential Energy (kJ/mole)": [-2000.0, -2100.0, -2150.0, -2160.0, -2162.0, -2160.0],
            "Box Volume (nm^3)": [120.0, 118.0, 117.5, 117.2, 117.1, 117.1],
        }
    )
    df_npt.to_csv(tmp_path / "npt.log", index=False)

    passed = analyze_logs(str(tmp_path))
    assert passed is True

    plot_equilibration(str(tmp_path))
    # Check generated artifacts
    assert (tmp_path / "figures/equilibration_qc.png").exists()
    assert (tmp_path / "equilibration_qc.json").exists()

    with open(tmp_path / "equilibration_qc.json") as f:
        data = json.load(f)
    assert data["passed"] is True
    assert 295.0 <= data["temperature_mean_k"] <= 305.0
    assert 0.95 <= data["density_mean_g_ml"] <= 1.05


def test_analyze_logs_varmdyn_passing(tmp_path: Path):
    df_nvt = pd.DataFrame(
        {
            "Time (ps)": [0.0, 10.0, 20.0, 30.0],
            "Temperature (K)": [299.5, 300.0, 300.1, 299.9],
            "Potential Energy (kJ/mole)": [-1800.0, -1850.0, -1860.0, -1855.0],
        }
    )
    df_nvt.to_csv(tmp_path / "nvt.log", index=False)

    df_npt = pd.DataFrame(
        {
            "Time (ps)": [0.0, 10.0, 20.0, 30.0],
            "Density (g/mL)": [0.999, 1.000, 1.001, 1.001],
            "Potential Energy (kJ/mole)": [-2150.0, -2160.0, -2162.0, -2160.0],
            "Box Volume (nm^3)": [117.5, 117.2, 117.1, 117.1],
        }
    )
    df_npt.to_csv(tmp_path / "npt.log", index=False)

    # Set protocol to varmdyn
    with open(tmp_path / "npt_protocol.json", "w") as f:
        json.dump({"protocol": "varmdyn"}, f)

    # Mock all 24 stages
    eq_dir = tmp_path / "equilibration"
    eq_dir.mkdir(parents=True, exist_ok=True)
    for st in range(1, 25):
        (eq_dir / f"stage_{st:02d}.json").write_text(json.dumps({"stage": st}))
        (eq_dir / f"stage_{st:02d}.state.xml").write_text("<State/>")

    passed = analyze_logs(str(tmp_path))
    assert passed is True

    with open(tmp_path / "equilibration_qc.json") as f:
        data = json.load(f)
    assert data["passed"] is True
    assert data["protocol"] == "varmdyn"
    assert data["varmdyn_stages_verified"] == 24
    assert data["varmdyn_stages_complete"] is True


def test_analyze_logs_varmdyn_missing_stages(tmp_path: Path):
    df_nvt = pd.DataFrame(
        {
            "Time (ps)": [0.0, 10.0],
            "Temperature (K)": [300.0, 300.0],
            "Potential Energy (kJ/mole)": [-1850.0, -1855.0],
        }
    )
    df_nvt.to_csv(tmp_path / "nvt.log", index=False)

    df_npt = pd.DataFrame(
        {
            "Time (ps)": [0.0, 10.0],
            "Density (g/mL)": [1.000, 1.000],
            "Potential Energy (kJ/mole)": [-2160.0, -2160.0],
            "Box Volume (nm^3)": [117.1, 117.1],
        }
    )
    df_npt.to_csv(tmp_path / "npt.log", index=False)

    with open(tmp_path / "npt_protocol.json", "w") as f:
        json.dump({"protocol": "varmdyn"}, f)

    # Mock only 18 stages (stages 19-24 missing)
    eq_dir = tmp_path / "equilibration"
    eq_dir.mkdir(parents=True, exist_ok=True)
    for st in range(1, 19):
        (eq_dir / f"stage_{st:02d}.json").write_text(json.dumps({"stage": st}))
        (eq_dir / f"stage_{st:02d}.state.xml").write_text("<State/>")

    passed = analyze_logs(str(tmp_path))
    assert passed is False

    with open(tmp_path / "equilibration_qc.json") as f:
        data = json.load(f)
    assert data["passed"] is False
    assert data["protocol"] == "varmdyn"
    assert data["varmdyn_stages_verified"] == 18
    assert data["varmdyn_stages_complete"] is False
    assert data["varmdyn_missing_stages"] == [19, 20, 21, 22, 23, 24]
