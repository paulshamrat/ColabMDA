"""Quality Control and stability verification unit tests."""

import json
from pathlib import Path

import pandas as pd

from colabmda.openmm.modular.qc import analyze_logs, parse_em_log


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

    # Check generated artifacts
    assert (tmp_path / "equilibration_qc.png").exists()
    assert (tmp_path / "equilibration_qc.pdf").exists()
    assert (tmp_path / "equilibration_qc.json").exists()

    with open(tmp_path / "equilibration_qc.json") as f:
        data = json.load(f)
    assert data["passed"] is True
    assert 295.0 <= data["temperature_mean_k"] <= 305.0
    assert 0.95 <= data["density_mean_g_ml"] <= 1.05
