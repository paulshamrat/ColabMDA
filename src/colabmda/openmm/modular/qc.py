import argparse
import json
import math
import os

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
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

        print("▶ Running Comprehensive Stability Check & Plotting …")

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
            eq_dir = os.path.join(workdir, "equilibration")
            for st in range(1, 25):
                st_json = os.path.join(eq_dir, f"stage_{st:02d}.json")
                st_xml = os.path.join(eq_dir, f"stage_{st:02d}.state.xml")
                if os.path.exists(st_json) and os.path.exists(st_xml):
                    varmdyn_stages_found.append(st)
                else:
                    varmdyn_stages_missing.append(st)

        # Plotting
        fig, axes = plt.subplots(3, 2, figsize=(12, 15))

        # 1. NVT Temp
        axes[0, 0].plot(df_nvt["Time (ps)"], df_nvt["Temperature (K)"], color="tab:red")
        axes[0, 0].axhline(300, ls="--", color="black", alpha=0.5)
        axes[0, 0].set_title(f"NVT Temperature (Avg: {avg_temp:.1f}K)")
        axes[0, 0].set_ylabel("K")

        # 2. NVT Energy
        axes[0, 1].plot(df_nvt["Time (ps)"], df_nvt["Potential Energy (kJ/mole)"], color="tab:blue")
        axes[0, 1].set_title("NVT Potential Energy")
        axes[0, 1].set_ylabel("kJ/mole")

        # 3. NPT Density
        axes[1, 0].plot(df_npt["Time (ps)"], df_npt["Density (g/mL)"], color="tab:green")
        axes[1, 0].set_title(f"NPT Density (Avg: {avg_dens:.4f} g/mL)")
        axes[1, 0].set_ylabel("g/mL")

        # 4. NPT Energy
        axes[1, 1].plot(
            df_npt["Time (ps)"], df_npt["Potential Energy (kJ/mole)"], color="tab:purple"
        )
        axes[1, 1].set_title("NPT Potential Energy")
        axes[1, 1].set_ylabel("kJ/mole")

        # 5. NPT Volume
        axes[2, 0].plot(df_npt["Time (ps)"], df_npt["Box Volume (nm^3)"], color="tab:orange")
        axes[2, 0].set_title("NPT Box Volume")
        axes[2, 0].set_ylabel("nm^3")

        # 6. Status Info Text
        axes[2, 1].axis("off")
        proto_line = f"Protocol    : {effective_protocol}\n" if effective_protocol else ""
        varmdyn_line = (
            f"VarMDyn 1-24: {len(varmdyn_stages_found)}/24 stages\n"
            if effective_protocol == "varmdyn"
            else ""
        )
        status_text = (
            f"EQUILIBRATION SUMMARY\n"
            f"---------------------\n"
            f"{proto_line}"
            f"{varmdyn_line}"
            f"EM Status   : {parse_em_log(workdir)}\n"
            f"NVT Temp    : {avg_temp:.2f} K (±{std_temp:.2f})\n"
            f"NPT Density : {avg_dens:.4f} g/mL\n"
            f"Dens. Slope : {slope:.2e}\n"
        )
        axes[2, 1].text(0.1, 0.5, status_text, fontsize=12, family="monospace", va="center")

        plt.tight_layout()
        plt.savefig("equilibration_qc.png", dpi=300)
        plt.savefig("equilibration_qc.pdf")
        print(f"  • NVT Temp: {avg_temp:.1f}K | NPT Density: {avg_dens:.4f} | Slope: {slope:.2e}")
        print("  • Saved comprehensive equilibration_qc.png")

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

        if effective_protocol == "varmdyn":
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
            qc_data["varmdyn_stages_total"] = 24
            qc_data["varmdyn_stages_verified"] = len(varmdyn_stages_found)
            qc_data["varmdyn_stages_complete"] = len(varmdyn_stages_missing) == 0
            if varmdyn_stages_missing:
                qc_data["varmdyn_missing_stages"] = varmdyn_stages_missing

        with open("equilibration_qc.json", "w") as handle:
            json.dump(qc_data, handle, indent=2)
            handle.write("\n")

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
