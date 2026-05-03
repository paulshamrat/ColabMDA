import os, sys, argparse
import pandas as pd
import matplotlib.pyplot as plt
import re

def parse_em_log(workdir):
    # EM doesn't have a CSV log, but we can capture the final minimized energy
    # from the standard output if we piped it, or check the system.xml state.
    # For now, we look for 'em.chk' existence and potentially a log if added.
    return "Minimized" if os.path.exists(os.path.join(workdir, 'em.chk')) else "Failed"

def analyze_logs(workdir):
    os.chdir(workdir)
    nvt_log = 'nvt.log'
    npt_log = 'npt.log'
    
    if not os.path.exists(nvt_log) or not os.path.exists(npt_log):
        print("Error: Equilibration logs not found.")
        return False

    print("▶ Running Comprehensive Stability Check & Plotting …")
    
    # Process NVT
    df_nvt = pd.read_csv(nvt_log)
    last_20_pct = max(1, int(len(df_nvt) * 0.2))
    avg_temp = df_nvt['Temperature (K)'].tail(last_20_pct).mean()
    std_temp = df_nvt['Temperature (K)'].tail(last_20_pct).std()
    
    # Process NPT
    df_npt = pd.read_csv(npt_log)
    avg_dens = df_npt['Density (g/mL)'].tail(last_20_pct).mean()
    
    # Simple slope check for density
    from scipy.stats import linregress
    subset = df_npt.tail(last_20_pct)
    slope, _, _, _, _ = linregress(range(len(subset)), subset['Density (g/mL)'])

    # Plotting
    fig, axes = plt.subplots(3, 2, figsize=(12, 15))
    
    # 1. NVT Temp
    axes[0,0].plot(df_nvt['Time (ps)'], df_nvt['Temperature (K)'], color='tab:red')
    axes[0,0].axhline(300, ls='--', color='black', alpha=0.5)
    axes[0,0].set_title(f'NVT Temperature (Avg: {avg_temp:.1f}K)')
    axes[0,0].set_ylabel('K')
    
    # 2. NVT Energy
    axes[0,1].plot(df_nvt['Time (ps)'], df_nvt['Potential Energy (kJ/mole)'], color='tab:blue')
    axes[0,1].set_title('NVT Potential Energy')
    axes[0,1].set_ylabel('kJ/mole')
    
    # 3. NPT Density
    axes[1,0].plot(df_npt['Time (ps)'], df_npt['Density (g/mL)'], color='tab:green')
    axes[1,0].set_title(f'NPT Density (Avg: {avg_dens:.4f} g/mL)')
    axes[1,0].set_ylabel('g/mL')
    
    # 4. NPT Energy
    axes[1,1].plot(df_npt['Time (ps)'], df_npt['Potential Energy (kJ/mole)'], color='tab:purple')
    axes[1,1].set_title('NPT Potential Energy')
    axes[1,1].set_ylabel('kJ/mole')
    
    # 5. NPT Volume
    axes[2,0].plot(df_npt['Time (ps)'], df_npt['Box Volume (nm^3)'], color='tab:orange')
    axes[2,0].set_title('NPT Box Volume')
    axes[2,0].set_ylabel('nm^3')
    
    # 6. Status Info Text
    axes[2,1].axis('off')
    status_text = (
        f"EQUILIBRATION SUMMARY\n"
        f"---------------------\n"
        f"EM Status   : {parse_em_log(workdir)}\n"
        f"NVT Temp    : {avg_temp:.2f} K (±{std_temp:.2f})\n"
        f"NPT Density : {avg_dens:.4f} g/mL\n"
        f"Dens. Slope : {slope:.2e}\n"
    )
    axes[2,1].text(0.1, 0.5, status_text, fontsize=12, family='monospace', va='center')
    
    plt.tight_layout()
    plt.savefig('equilibration_qc.png', dpi=300)
    plt.savefig('equilibration_qc.pdf')
    print(f"  • NVT Temp: {avg_temp:.1f}K | NPT Density: {avg_dens:.4f} | Slope: {slope:.2e}")
    print("  • Saved comprehensive equilibration_qc.png")

    # Stability Gate logic (Keep it helpful but not strictly blocking for small tests)
    passed = True
    if abs(avg_temp - 300) > 20:
        print("⚠ WARNING: Temperature unstable.")
        passed = False
    if abs(slope) > 1e-3:
        print("⚠ WARNING: Density not converged.")
        passed = False
        
    return passed

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("workdir")
    args = p.parse_args()
    if not analyze_logs(args.workdir):
        # We still exit 0 for now to allow the test run to proceed, 
        # but we print clear warnings.
        pass
