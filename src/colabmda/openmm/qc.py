import argparse
import json
import math
import os

import pandas as pd
from scipy.stats import linregress


def parse_em_log(workdir):
    # EM doesn't have a CSV log, but we can capture the final minimized energy
    # from the standard output if we piped it, or check the system.xml state.
    # For now, we look for 'em.chk' existence and potentially a log if added.
    return "Minimized" if os.path.exists(os.path.join(workdir, "em.chk")) else "Failed"


def analyze_logs(workdir, strict=False, protocol=None):
    workdir = os.path.abspath(workdir)
    old_cwd = os.getcwd()
    try:
        os.chdir(workdir)
        nvt_log = "nvt.log"
        npt_log = "npt.log"

        print("▶ Running Stability Check …")

        # Process NVT
        df_nvt = pd.read_csv(nvt_log)
        last_20_pct_nvt = max(2 if len(df_nvt) >= 2 else 1, int(len(df_nvt) * 0.2))
        avg_temp = df_nvt["Temperature (K)"].tail(last_20_pct_nvt).mean()
        std_temp = df_nvt["Temperature (K)"].tail(last_20_pct_nvt).std()
        if pd.isna(std_temp):
            std_temp = 0.0

        # Process NPT
        df_npt = pd.read_csv(npt_log)
        last_20_pct_npt = max(2 if len(df_npt) >= 2 else 1, int(len(df_npt) * 0.2))
        avg_dens = df_npt["Density (g/mL)"].tail(last_20_pct_npt).mean()

        subset = df_npt.tail(last_20_pct_npt)
        if len(subset) >= 2:
            slope, _, _, _, _ = linregress(subset["Time (ps)"], subset["Density (g/mL)"])
        else:
            slope = 0.0 if len(subset) == 1 else math.nan

        # Detect protocol
        effective_protocol = protocol
        if effective_protocol is None:
            for proto_candidate in ["npt_protocol.json", "nvt_protocol.json"]:
                if os.path.exists(proto_candidate):
                    try:
                        with open(proto_candidate) as pf:
                            effective_protocol = json.load(pf).get("protocol")
                        if effective_protocol:
                            break
                    except Exception:
                        pass

        # Check 1-24 stages if varmdyn
        varmdyn_stages_found = []
        varmdyn_stages_missing = []
        if effective_protocol == "varmdyn":
            stages_dir = os.path.join(workdir, "stages")
            legacy_dir = os.path.join(workdir, "equilibration")
            eq_dir = stages_dir if os.path.exists(stages_dir) else legacy_dir
            for st in range(1, 25):
                candidates = [
                    os.path.join(eq_dir, f"state_{st:02d}.xml"),
                    os.path.join(eq_dir, f"stage_{st:02d}.state.xml"),
                ]
                st_xml = next(
                    (p for p in candidates if os.path.isfile(p) and os.path.getsize(p) > 0), None
                )
                if st_xml is not None:
                    varmdyn_stages_found.append(st)
                else:
                    varmdyn_stages_missing.append(st)

        # Plotting
        passed = all(math.isfinite(value) for value in (avg_temp, avg_dens, slope))
        reasons = []
        if not passed:
            reasons.append("non-finite or insufficient thermodynamic data")
        if abs(avg_temp - 300) > 10:
            print("⚠ WARNING: Temperature unstable.")
            passed = False
            reasons.append("mean final temperature differs from 300 K by more than 10 K")
        if not 0.90 <= avg_dens <= 1.10:
            print("⚠ WARNING: Density outside the expected explicit-water range.")
            passed = False
            reasons.append("mean final density is outside 0.90-1.10 g/mL")
        if abs(slope) > 1e-4:
            print("⚠ WARNING: Density not converged.")
            passed = False
            reasons.append("absolute final density slope exceeds 1e-4 g/mL/ps")

        invalid_stage_timing = []
        if effective_protocol == "varmdyn":
            for stage in varmdyn_stages_found:
                expected = 0
                if 4 <= stage <= 9:
                    expected = (stage - 3) * 10
                elif 10 <= stage <= 14:
                    expected = 60
                elif 15 <= stage <= 23:
                    expected = 60 + (stage - 14) * 100
                elif stage == 24:
                    expected = 1960
                try:
                    import xml.etree.ElementTree as ET

                    candidates = [
                        os.path.join(eq_dir, f"state_{stage:02d}.xml"),
                        os.path.join(eq_dir, f"stage_{stage:02d}.state.xml"),
                    ]
                    st_xml = next(
                        (p for p in candidates if os.path.isfile(p) and os.path.getsize(p) > 0),
                        None,
                    )
                    xml_time = float(ET.parse(st_xml).getroot().attrib["time"]) if st_xml else -1.0

                    stem = os.path.join(eq_dir, f"stage_{stage:02d}")
                    if os.path.exists(stem + ".json"):
                        with open(stem + ".json") as handle:
                            metadata = json.load(handle)
                        valid = (
                            not metadata.get("test_mode", False)
                            and math.isfinite(metadata["potential_energy_kj_mol"])
                            and abs(metadata["time_ps"] - expected) < 0.001
                            and abs(xml_time - expected) < 0.001
                        )
                    else:
                        is_test = os.environ.get("COLABMDA_TEST_1_24", "").lower() in {
                            "1",
                            "true",
                        }
                        valid = not is_test and abs(xml_time - expected) < 0.001
                except (KeyError, ValueError, OSError, ET.ParseError):
                    valid = False
                if not valid:
                    invalid_stage_timing.append(stage)
            if strict and invalid_stage_timing:
                passed = False
                reasons.append(f"staged timing/energy evidence invalid: {invalid_stage_timing}")
            if varmdyn_stages_missing:
                print(
                    f"⚠ WARNING: Missing {len(varmdyn_stages_missing)}/24 VarMDyn stages: {varmdyn_stages_missing}"
                )
                passed = False
                reasons.append(
                    f"missing {len(varmdyn_stages_missing)}/24 VarMDyn stages: {varmdyn_stages_missing}"
                )
            else:
                print("  • Verified all 24/24 VarMDyn equilibration stages")

        qc_data = {
            "passed": passed,
            "protocol": effective_protocol,
            "temperature_mean_k": float(avg_temp),
            "temperature_std_k": float(std_temp),
            "density_mean_g_ml": float(avg_dens),
            "density_slope_g_ml_ps": float(slope),
            "reasons": reasons,
        }
        if effective_protocol == "varmdyn":
            qc_data["varmdyn_timing_verified"] = (
                not invalid_stage_timing and not varmdyn_stages_missing
            )
            qc_data["varmdyn_invalid_timing_stages"] = invalid_stage_timing
            qc_data["varmdyn_stages_total"] = 24
            qc_data["varmdyn_stages_verified"] = len(varmdyn_stages_found)
            qc_data["varmdyn_stages_complete"] = len(varmdyn_stages_missing) == 0
            if varmdyn_stages_missing:
                qc_data["varmdyn_missing_stages"] = varmdyn_stages_missing

        with open("equilibration_qc.json", "w") as handle:
            json.dump(qc_data, handle, indent=2)
            handle.write("\n")

        with open("qc.txt", "w") as f:
            f.write(f"Status: {'PASSED' if passed else 'FAILED'}\n")
            f.write(f"Temperature: {avg_temp:.2f} +/- {std_temp:.2f} K\n")
            f.write(f"Density: {avg_dens:.4f} g/mL (slope: {slope:.6f} g/mL/ps)\n")

        return passed
    finally:
        os.chdir(old_cwd)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("workdir")
    p.add_argument("--protocol", default=None)
    p.add_argument("--strict", action="store_true")
    args = p.parse_args()
    if not analyze_logs(args.workdir, strict=args.strict, protocol=args.protocol):
        raise SystemExit(1)
