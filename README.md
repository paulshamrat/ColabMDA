# ColabMDA

**ColabMDA** is a specialized tool that lets you run high-quality Molecular Dynamics simulations on Google Colab without the fear of losing your work. The biggest problem with Colab is that it disconnects, often destroying hours of simulation data. **ColabMDA fixes this** with a "resume-safe" system that automatically saves your progress to Google Drive. If your session expires, you can resume exactly where you left off with one simple command. From modeling protein mutations to generating publication-ready analysis, ColabMDA handles the complex setup for you, making it the easiest way to get high-quality MD results using free cloud GPUs.

📖 **Full Documentation:** Visit our official manual at [colabmda.readthedocs.io](https://colabmda.readthedocs.io/)

---

## 🛠 Project Information
| Category | Details |
| :--- | :--- |
| **Release** | [![GitHub tag](https://img.shields.io/github/v/tag/paulshamrat/ColabMDA)](https://github.com/paulshamrat/ColabMDA/tags) |
| **Availability** | [![GitHub](https://img.shields.io/badge/GitHub-Repository-blue?logo=github)](https://github.com/paulshamrat/ColabMDA) |
| **Documentation** | [![Documentation Status](https://readthedocs.org/projects/colabmda/badge/?version=latest)](https://colabmda.readthedocs.io/en/latest/?badge=latest) |
| **Workflows** | [![Python CI](https://github.com/paulshamrat/ColabMDA/actions/workflows/python-test.yml/badge.svg?branch=main)](https://github.com/paulshamrat/ColabMDA/actions/workflows/python-test.yml) |
| **Issues** | [![GitHub issues](https://img.shields.io/github/issues/paulshamrat/ColabMDA)](https://github.com/paulshamrat/ColabMDA/issues) |
| **License** | [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) |
| **Style / Lint** | [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black) [![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff) |
| **Dependencies** | `OpenMM`, `AmberTools (tleap, cpptraj)`, `Modeller`, `MDAnalysis`, `MDTraj`, `PDB2PQR / PROPKA`, `OpenFF (Optional)` |
| **Platform** | `Linux` · `Google Colab (CLI / Terminal)` · `HPC (SLURM)` |
| **Structure** | [Source Layout](#-project-structure) |

### 📂 Project Structure
ColabMDA is organized into clear, functional modules:
*   `src/colabmda/cli.py`: Central unified CLI (`colabmda` and shortcut `cmda`).
*   `src/colabmda/core/`: System taxonomy (`SystemSpec`), project layout (`ProjectPaths`), and divalent ion parameters.
*   `src/colabmda/amber/`: AmberTools parameterization (`tleap`, Hu 2024 quantum nucleotide parameters, `cpptraj` autoimaging & solvent stripping).
*   `src/colabmda/modeller/`: Homology modeling engine (UniProt biological numbering, mutagenesis, geometric substrate docking).
*   `src/colabmda/openmm/`: OpenMM simulation engines (Modular prep/EM/NVT/NPT/QC/MD/analysis).
*   `src/colabmda/remote/`: Google Colab GPU remote session management via `colab` CLI.
*   `scripts/`: Single standalone bootstrap script (`scripts/setup.sh`).
*   `tests/`: Unit and regression test suite (`pytest`).
*   `docs/`: Sphinx and ReadTheDocs scientific documentation source.

---

## 1. Quick Start (Google Colab Installation)

> 💡 **Terminal Access:** All bash commands should be run in the **Colab Terminal** (Open via the **⋮** menu -> **Terminal**).

### 1.1. Setup Colab Runtime & Drive
Before starting, ensure your environment is ready:
1. **Enable GPU:** Go to `Runtime` -> `Change runtime type` and select **T4 GPU**.
2. **Verify GPU:** Run `!nvidia-smi` in a cell to confirm GPU access.
3. **Mount Drive:** Click the **Folder icon** 📂 in the left sidebar, then click the **Drive icon** (Mount Drive), or run the code block below:

```python
from google.colab import drive
drive.mount('/content/drive')
```

### 1.2. Environment & Package Installation (Required)
Run the following in the **Colab Terminal** (⋮ → Terminal).
*Estimated time: ~3–5 minutes.*

```bash
# 1. Download and run the unified bootstrap script (provisions OpenMM, AmberTools, MDAnalysis, and ColabMDA)
cd /content
curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/main/scripts/setup.sh -o setup.sh
bash setup.sh latest

# Optional: To include MODELLER for homology modeling & mutagenesis:
# WITH_MODELLER=1 bash setup.sh latest
```

> [!IMPORTANT]
> **Modeller License (Optional):** If you enable `WITH_MODELLER=1`, the script will prompt for your Modeller License Key.
> 🔑 **Get a Free License:** Register at [salilab.org/modeller/registration.html](https://salilab.org/modeller/registration.html) (Academic licenses are free and sent instantly via email).

---

### 💡 Tip: How to Resume After a Timeout
If your Google Colab session expires:
1. Re-run **Required Step 1.2** (`bash setup.sh latest`) to restore the environment.
2. Run your simulation command (`cmda sim ...` or `colabmda openmm run ...`).
3. The engine automatically discovers existing `.chk` files and state XMLs and resumes seamless execution.

---

## 2. Workstation & HPC Installation

<details>
<summary><b>💻 A. Local Workstation Setup (Conda / Mamba)</b></summary>

ColabMDA runs natively on Linux workstations with NVIDIA GPUs. Use the provided `environment.yml` to provision the production environment:
```bash
git clone https://github.com/paulshamrat/ColabMDA.git
cd ColabMDA
mamba env create -f environment.yml
conda activate colabmda
```
</details>

<details>
<summary><b>🧪 B. Running the Test Suite</b></summary>

To verify force fields, peptide docking, Amber leap parameterization, and CLI execution:
```bash
pytest tests/ -v
```
</details>

<details>
<summary><b>🏢 C. HPC Usage (SLURM)</b></summary>

ColabMDA generates decoupled SLURM batch scripts automatically under `<project>/scripts/submit.slurm`. Since trajectories are chunked with automatic checkpointing, jobs run reliably within cluster walltime limits.
</details>

---

## 3. Recommended Directory Strategy & Architectural Principles

ColabMDA enforces strict architectural separation to ensure reproducibility across local workstations, cloud GPUs, and HPC clusters:

1. **Zero Data in `src/`:** Structural inputs (`*.pdb`), topologies (`*.prmtop`), and trajectories belong strictly in the user's project workspace. Non-standard residues like phosphotyrosine (`pY171` -> `PTR171`) are constructed mathematically with zero hardcoded PDBs in source code.
2. **Decoupled Execution Scripts:** Execution shell scripts and SLURM batch scripts reside exclusively in `<project>/scripts/`, never inside simulation data folders.
3. **Single Bootstrap Script:** Standalone environment provisioning is handled by a single script: `scripts/setup.sh`.

```text
research/projects/<project_name>/         # Your project root (e.g. cdkl5, kras, 4ldj)
├── structure/                            # Input structures (Zero raw PDBs in src/)
│   ├── protein/                          # Target receptor structures (e.g. cdkl5_wt.pdb, cdkl5_r59w.pdb)
│   └── template/                         # Structural templates (e.g. 3X2V.pdb)
├── simulation/                           # Staged & executed simulation data
│   └── <state>/<variant>/                # E.g. atp_2mg_eb2/WT/
│       ├── eq/                           # Shared equilibration hub (system.prmtop, em, nvt, npt, qc)
│       ├── r1/                           # Production replica 1 (chunks + trajectory)
│       ├── r2/                           # Production replica 2 (independent seed from eq/)
│       └── r3/                           # Production replica 3
├── scripts/                              # Decoupled execution scripts
│   ├── run.sh                            # Standalone runner for local/remote execution
│   └── submit.slurm                      # HPC SLURM batch submission script
└── analysis/                             # Trajectory metrics, CSVs, and publication plots
    └── <state>/<variant>/                # E.g. atp_2mg_eb2/WT/ (rmsd.csv, rmsf.csv, rg.csv)
```

> 💡 **Why Shared `eq/`?**
> Equilibration (EM, NVT, NPT, QC) only needs to be run once per system. Replicas (`r1`, `r2`, `r3`) fork directly from `eq/npt.state.xml` with distinct random velocity seeds, saving substantial GPU time while maintaining statistical independence.

---

## 4. Unified Kinase-Substrate Simulation Engine (`cmda`)

ColabMDA provides a streamlined, single-keyword CLI accessible via `cmda` (or `colabmda`) that unifies structure preparation, parameterization, OpenMM simulation dispatch, trajectory post-processing, and multi-replica analysis.

### 4.1. Core Capabilities: Kinase-Substrate Dynamics
* **Phosphotyrosine Modeling (`pY171`):** Mathematically converts target tyrosine residues to Amber `PTR` with sp3-hybridized oxygens and $-2.0e$ net charge, dynamically routing Amber `phosaa19SB` (or `phosaa14SB`) force fields.
* **Quantum Nucleotide Parameters:** Integrates Hu et al. 2024 quantum-calculated parameters for ATP, ADP, and GDP.
* **Dual Divalent Metal Ions:** Supports catalytic dual $Mg^{2+}$ or $Mn^{2+}$ coordination using Amber 12-6-4 / standard LJ parameters.
* **Substrate Peptide Docking:** Docks substrate peptides (e.g. EB2 `KRPSSAPSG` at Ser222) via continuous geometric backbone vector propagation ($d_{\text{C-N}} = 1.33\text{ Å}$), eliminating atom overlap and NaN minimization errors.
* **Ultra-Fast Solvation & Neutralization:** Uses Amber `tleap` with randomized ion placement (`addIonsRand`), solvating systems in under 2 seconds.
* **Dual-Stream Trajectory Stripping:** Automates `cpptraj` PBC autoimaging and strips bulk water, producing lightweight dry (`keep_ligand`) and hydration-shell trajectories (>90% storage savings).

---

### 4.2. Single-Keyword Command Reference (`cmda`)

| Command | Syntax | Description |
| :--- | :--- | :--- |
| **`stage`** | `cmda stage --state <tag> --variant <var>` | Generates decoupled project trees and runner scripts (`run.sh`, `submit.slurm`). |
| **`prep`** | `cmda prep --state <tag> --variant <var>` | Assembles coordinates and parameterizes Amber topology (`system.prmtop`, `system.inpcrd`). |
| **`sim`** | `cmda sim --state <tag> --variant <var> --stage <stages>` | Dispatches OpenMM simulations (`em`, `nvt`, `npt`, `qc`, `md`). |
| **`merge`** | `cmda merge --state <tag> --variant <var>` | Wraps PBC and strips bulk solvent via `cpptraj` (`stripped.system.prmtop`). |
| **`analyze`** | `cmda analyze --state <tag> --variant <var>` | Computes RMSD, $R_g$, RMSF, PCA, and contact maps. |
| **`compare`** | `cmda compare --state <tag> --variants <v1,v2>` | Renders multi-replica or cross-variant comparative overlay plots. |
| **`colab`** | `cmda colab <status\|new\|run\|stop>` | Controls remote Google Colab GPU sessions directly from the terminal. |

#### Complete End-to-End Walkthrough:
```bash
# 1. Stage the decoupled directory structure and scripts
cmda stage --state atp_2mg_eb2 --variant WT --root research/projects/cdkl5

# 2. Build coordinates and Amber topology (pY171 + ATP + 2 Mg2+ + EB2 peptide)
cmda prep --state atp_2mg_eb2 --variant WT --root research/projects/cdkl5

# 3. Run equilibration and QC gate in shared eq/ hub
cmda sim --state atp_2mg_eb2 --variant WT --stage em,nvt,npt,qc --protocol varmdyn --root research/projects/cdkl5

# 4. Run resume-safe production MD across triplicates (r1, r2, r3)
cmda sim --state atp_2mg_eb2 --variant WT --stage md --replica r1 --total-ns 100.0 --root research/projects/cdkl5
cmda sim --state atp_2mg_eb2 --variant WT --stage md --replica r2 --total-ns 100.0 --root research/projects/cdkl5
cmda sim --state atp_2mg_eb2 --variant WT --stage md --replica r3 --total-ns 100.0 --root research/projects/cdkl5

# 5. Autoimage PBC and strip bulk solvent with cpptraj
cmda merge --state atp_2mg_eb2 --variant WT --replicas r1,r2,r3 --root research/projects/cdkl5

# 6. Extract trajectory metrics and plots
cmda analyze --state atp_2mg_eb2 --variant WT --root research/projects/cdkl5

# 7. Compare Wild-Type vs Mutant (e.g. R59W)
cmda compare --state atp_2mg_eb2 --variants WT,R59W --root research/projects/cdkl5
```

---

## 5. Granular & Legacy Workflows (`colabmda openmm` / `modeller`)

For users requiring fine-grained control over individual modeling steps or legacy PDB workflows:

### 5.1. Homology Modeling & Point Mutations (Modeller)
**Environment:** `modeller_env`

The Modeller workflow includes **Biological Residue Numbering** and **Automated Quality Control**:
*   **`--uniprot-numbering`**: Physically renumbers residues to match UniProt biological indexing.
*   **Automatic Alignment Summary**: Confirms sequence alignment before building.
*   **Post-Build QC**: Verifies modeled coordinates and stereochemical energies.

```bash
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate modeller_env

# 1. Build Wild-Type KRAS (Starting at Residue 1 matching UniProt)
colabmda modeller build --pdb-id 4ldj --uniprot-id P01116 --chain A --range 1 169 --uniprot-numbering --outdir data/structure/4ldj/wt

# 2. Create Point Mutant (e.g. G12D) preserving biological numbering
colabmda modeller mutate --pdb-in data/structure/4ldj/wt/target.B99990001_with_cryst.pdb --chain A --mut G12D --outdir-mut data/structure/4ldj/mutants/4ldj_g12d
```

---

### 5.2. Structure Preparation & 3D pKa Titration (PROPKA at pH 7.4)
**Environment:** `openmm_env`

Accurate protonation states are critical for biophysical simulation fidelity. ColabMDA natively integrates **PDBFixer**, **PDB2PQR**, and **PROPKA 3** to calculate local 3D electrostatic microenvironments and assign explicit AMBER titratable residue forms (`HID`, `HIE`, `HIP`, `ASH`, `GLH`, `LYN`, `CYX`):

```bash
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate openmm_env

# Run 3D pKa titration at physiological pH 7.4
colabmda openmm prep --pdb-file data/structure/4ldj/wt/4ldj_wt.pdb --outdir data/structure/4ldj/wt/prep --name 4ldj_wt --ph 7.4
colabmda openmm prep --pdb-file data/structure/4ldj/mutants/4ldj_g12d/4ldj_g12d.pdb --outdir data/structure/4ldj/mutants/4ldj_g12d/prep --name 4ldj_g12d --ph 7.4
```

*Key features:*
* **Automated Coordinate Repair:** PDBFixer resolves missing heavy atoms, non-standard residues, and chain breaks.
* **3D Electrostatic Titration:** Evaluates local hydrogen-bonding networks using PROPKA 3 at target pH (default: `7.4`).
* **Explicit AMBER Forms:** Differentiates neutral Histidine tautomers (`HID`/`HIE`), charged Histidine (`HIP`), protonated Asp/Glu (`ASH`/`GLH`), neutral Lysine (`LYN`), and disulfide bridges (`CYX`).
* **Audit Manifest:** Produces `protonation_summary.json` alongside `<name>_cleaned.pdb` recording every assigned titration state.

---

### 5.3. Simulation Execution: Modular Stages vs. All-in-One

ColabMDA provides both fine-grained single-keyword commands and an automated all-in-one runner:

#### Single-Keyword Modular Stages

| Stage | Command | Description |
| :--- | :--- | :--- |
| **Stage** | `colabmda openmm stage --pdb-file <file> --name <system> --replica r1` | Initialize folder structure in `data/simulation/` and `eq/`. |
| **EM** | `colabmda openmm em --name <system> --protein-ff ff19SB --water-model opc` | Energy minimization & box solvation in `eq/`. |
| **NVT** | `colabmda openmm nvt --name <system> --equil-time 100` | NVT thermalization to 300 K. |
| **NPT** | `colabmda openmm npt --name <system> --equil-time 100` | NPT density equilibration with Monte Carlo barostat. |
| **QC** | `colabmda openmm qc --name <system>` | Stability verification & publication QC plots (`qc.png`). |
| **MD** | `colabmda openmm md --name <system> --replica r1 --total-ns 100 --traj-interval 10` | Resume-safe chunked production MD. |
| **Merge** | `colabmda openmm merge --name <system> --replica r1 --center --wrap` | Concatenate chunks into centered, PBC-wrapped trajectory. |
| **Status** | `colabmda openmm status --name <system> --replica r1` | Verify chunks, topology, frame counts, and resume readiness. |
| **Analysis** | `colabmda openmm analysis --name <system> --replica r1` | Single-system trajectory metrics (RMSD, $R_g$, RMSF). |

#### Step-by-Step Modular Execution:

```bash
# 1. Initialize simulation folder and eq hub
colabmda openmm stage --pdb-file data/structure/4ldj/wt/prep/4ldj_wt_cleaned.pdb --name 4ldj_wt --replica r1

# 2. Energy Minimization in shared eq/ hub (ff19SB / OPC 4-point water)
colabmda openmm em --name 4ldj_wt --protein-ff ff19SB --water-model opc

# 3. NVT Thermalization (100 ps heating to 300 K)
colabmda openmm nvt --name 4ldj_wt --equil-time 100

# 4. NPT Density Equilibration (100 ps with Monte Carlo barostat)
colabmda openmm npt --name 4ldj_wt --equil-time 100

# 5. Equilibration Quality Control Gating
colabmda openmm qc --name 4ldj_wt

# 6. Production MD across Triplicates (r1, r2, r3 fork automatically from eq/)
colabmda openmm md --name 4ldj_wt --replica r1 --total-ns 100.0 --traj-interval 10
colabmda openmm md --name 4ldj_wt --replica r2 --total-ns 100.0 --traj-interval 10
colabmda openmm md --name 4ldj_wt --replica r3 --total-ns 100.0 --traj-interval 10

# 7. Merge Trajectory Chunks for all replicas (Center + Wrap)
colabmda openmm merge --name 4ldj_wt --replica r1 --center --wrap
colabmda openmm merge --name 4ldj_wt --replica r2 --center --wrap
colabmda openmm merge --name 4ldj_wt --replica r3 --center --wrap

# 8. Check Status & Frame Counts
colabmda openmm status --name 4ldj_wt --replica r1

# 9. Biophysical Trajectory Analysis (RMSD, Rg, RMSF)
colabmda openmm analysis --name 4ldj_wt --replica r1
colabmda openmm analysis --name 4ldj_wt --replica r2
colabmda openmm analysis --name 4ldj_wt --replica r3
```

#### All-in-One Alternative (`run`):
For automated, hands-off staging and execution from EM through Production MD in one call:
```bash
colabmda openmm run --name 4ldj_wt --replica r1 --total-ns 100.0 --traj-interval 10 --equil-time 1000 --checkpoint-ps 1000
```

> 💡 **Storage Tip:** For a typical system (e.g., KRAS in water, ~30,000 atoms), a 100ns run at high resolution (1ps) can produce over **36GB** of data. On a free 15GB Google Drive, we recommend using **`--traj-interval 10`** to reduce this to ~3.6GB. Always calculate your storage needs based on your specific system size before starting long runs.

> 🛡️ **Stability Gate:** The `qc` command (and built-in gating in `run`) automatically analyzes equilibration logs and aborts if the system has not stabilized, preventing wasted cloud GPU time.

---

### 5.4. Merge and Center Trajectories

Combine trajectory chunks into a single DCD, apply periodic boundary condition (PBC) correction, and center the protein using the robust MDAnalysis-based engine:

```bash
# Standard Merge (Center + Wrap, All Atoms)
colabmda openmm merge --name 4ldj_wt --replica r1 --center --wrap
colabmda openmm merge --name 4ldj_wt --replica r2 --center --wrap
colabmda openmm merge --name 4ldj_wt --replica r3 --center --wrap

# For Mutant (G12D):
colabmda openmm merge --name 4ldj_g12d --replica r1 --center --wrap
colabmda openmm merge --name 4ldj_g12d --replica r2 --center --wrap
colabmda openmm merge --name 4ldj_g12d --replica r3 --center --wrap
```

#### Output Files:
*   `prod_full.dcd`: The final, concatenated, and centered/wrapped trajectory file.
*   `prod_full.log`: The consolidated log file containing energy and temperature statistics.

> 💡 **Pro-Tip for Long Runs:**
> Merging processes trajectories frame-by-frame, so it will not crash your RAM. You can merge without striding (`--stride 1`) for full resolution, or use `--stride 10` to create a lightweight file for local viewing.
>
> When you merge with `--stride 10`, the logs are automatically strided to match, and the analysis command (`colabmda openmm analysis`) will read the logs to correctly infer the time scale.

> 💡 **Protein-Only Trajectory (Optional):**
> To save >85% storage, add the `--protein-only` flag to extract only protein atoms, discarding solvent water and ions:
> `colabmda openmm merge --name 4ldj_wt --replica r1 --center --wrap --protein-only`
>
> **Outputs with `--protein-only`:**
> * `prod_full.dcd`: Thinned trajectory containing protein atoms only.
> * `prod_full.pdb`: Matching protein-only topology PDB file (crucial for subsequent analysis/visualization to avoid atom mismatch).
> * `prod_full.log`: Consolidated log file.

---

### 5.5. Post-Merge Guidelines & Verification (FAQ)

#### 5.5.1. 🔍 Verifying Trajectory Frames
To quickly check the status, simulation time, and number of frames in your merged files, run:
```bash
colabmda openmm status --name 4ldj_wt --replica r1
```
This prints a comprehensive status report including topology stats, chunks, and exact frame counts:
```text
[STATUS]
  Workdir          : data/simulation/4ldj_wt/r1
  Chunks (DCD/log) : 100 / 100
  Topology File    : solvated.pdb
                     └─ 27273 atoms, 8388 residues
                     └─ (169 protein, 8170 water, 49 ions)
  Trajectory File  : prod_full.dcd (10000 frames)
  Log File         : prod_full.log (10000 frames)
  Frames (from logs): 10000
```

Alternatively, query the frame count directly with MDAnalysis:
```bash
python3 -c "import MDAnalysis as mda; u = mda.Universe('data/simulation/4ldj_wt/r1/prod_full.pdb', 'data/simulation/4ldj_wt/r1/prod_full.dcd'); print('Frames:', len(u.trajectory))"
```

#### 5.5.2. ⏱️ Understanding Simulation Time vs. Frame Counts
Use this formula to calculate expected frames:
$$\text{Expected Frames} = \frac{\text{Total Time (ps)}}{\text{Trajectory Saving Interval (ps)}}$$
*Example:* A $100.0\text{ ns}$ ($100,000\text{ ps}$) simulation with a $10.0\text{ ps}$ saving interval:
$$\text{Expected Frames} = \frac{100,000\text{ ps}}{10\text{ ps}} = 10,000\text{ frames}$$

*Striding on Merged Trajectories:*
$$\text{Merged Frames} = \frac{\text{Total Frames}}{\text{Stride}} = \frac{10,000}{10} = 1,000\text{ frames}$$

#### 5.5.3. 🧬 Reference Topology Guidelines
When merging trajectories, the reference topology file used and the resulting output files depend on your merge options:

| Merge Mode | Command Flags | Reference Topology | Generated Topology File | Subsequent Command Usage |
| :--- | :--- | :--- | :--- | :--- |
| **Standard Merge** | `colabmda openmm merge --center --wrap` | `solvated.pdb` | *None* (only `prod_full.dcd` is written) | Use `solvated.pdb` for analysis/visualization. |
| **MDAnalysis Merge** | `colabmda openmm merge --mda --center --wrap` | `solvated.pdb` | `prod_full.pdb` (all atoms) | Use `prod_full.pdb` (or `solvated.pdb`). |
| **Protein-Only Merge** | `colabmda openmm merge --protein-only` | `solvated.pdb` | `prod_full.pdb` (protein atoms only) | **Must** use `prod_full.pdb` (since `prod_full.dcd` contains only protein coordinates). |

> ⚠️ **Why does this matter?**
> If you merge using `--protein-only`, `prod_full.dcd` contains only protein coordinates (~2,600 atoms). Attempting to load this trajectory along with `solvated.pdb` (~27,000 atoms) in PyMOL or MDAnalysis will result in a **fatal atom mismatch error**. Always match the trajectory file with its corresponding topology file as shown in the table above.

---

### 5.6. Analysis & Comparison

#### 5.6.1. Single System Analysis (RMSD, Rg, RMSF)
```bash
# For Wild-Type Triplicates:
colabmda openmm analysis --name 4ldj_wt --replica r1
colabmda openmm analysis --name 4ldj_wt --replica r2
colabmda openmm analysis --name 4ldj_wt --replica r3

# For Mutant Triplicates (G12D):
colabmda openmm analysis --name 4ldj_g12d --replica r1
colabmda openmm analysis --name 4ldj_g12d --replica r2
colabmda openmm analysis --name 4ldj_g12d --replica r3
```

> ⚠️ **Analysis Tip:** If your plots show the wrong time scale (e.g., 10ns instead of 100ns), provide the frame interval manually with `--interval 10`.

#### 5.6.2. WT vs Mutant Comparison (Triplicate Statistics)
```bash
colabmda openmm compare \
  --series "WT=data/analysis/4ldj_wt/r1,data/analysis/4ldj_wt/r2,data/analysis/4ldj_wt/r3" \
  --series "G12D=data/analysis/4ldj_g12d/r1,data/analysis/4ldj_g12d/r2,data/analysis/4ldj_g12d/r3" \
  --outdir data/analysis/compare/wt_vs_g12d
```

---

## 6. Trajectory Visualization in PyMOL

To visualize molecular dynamics trajectories in PyMOL with full secondary structure (ribbon/cartoon) and detailed sidechain representations, use the built-in `colabmda openmm view` command.

### 6.1. Local Setup and Installation
Since molecular dynamics simulations are performed on Google Colab, trajectory visualization is run locally on your workstation/laptop.

> [!NOTE]
> **For Windows Users:** It is highly recommended to run visualization locally either using **WSL (Windows Subsystem for Linux)** or using native **PyMOL on Windows** (by running the CLI to generate `visualize.pml`, then opening it in the Windows PyMOL GUI).

Create a dedicated conda environment containing PyMOL and ColabMDA:

```bash
# 1. Create a clean environment and install PyMOL
conda create -y -n colabmda_env -c conda-forge python=3.11 pymol-open-source

# 2. Activate the environment
conda activate colabmda_env

# 3. Install ColabMDA package (Lightweight local installation)
python3 -m pip install --upgrade "git+https://github.com/paulshamrat/ColabMDA.git@main"

# 4. Verify PyMOL installation
pymol --version
```

### 6.2. Running the Visualization in PyMOL
Run the view command from your simulation folder (it looks for `prod_full.pdb` and `prod_full.dcd`, generates a PyMOL script, and launches PyMOL):

```bash
# For Wild-Type:
colabmda openmm view --pdb-dir data/simulation/4ldj_wt/r1

# For Mutant (G12D, highlighting residue 12):
colabmda openmm view --pdb-dir data/simulation/4ldj_g12d/r1 --resi 12
```

The command automatically generates a `visualize.pml` script in the folder.

<details>
<summary><b>Click to view the generated PyMOL script configuration (<code>visualize.pml</code>)</b></summary>

```text
# 1. Clear out memory
reinitialize
bg_color white

# 2. Load merged trajectory files
load prod_full.pdb, kras
load_traj prod_full.dcd, kras

# 2b. Align trajectory to crystal reference (e.g., 4ldj_wt.pdb)
load 4ldj_wt.pdb, crystal_ref
align kras and name CA, crystal_ref and name CA, mobile_state=1, target_state=1
delete crystal_ref

# 3. Hide solvent water and ions
hide everything, all
hide nonbonded, all
hide nb_spheres, all

# 4. Freeze backbone tumbling frame-by-frame
intra_fit kras and name CA

# 5. Generate crisp secondary structure ribbon
dss kras
cartoon automatic, kras
show cartoon, kras

# 6. Apply smooth cyan color scheme
color cyan, kras and name C*
util.cnc("kras")

# 7. Highlight mutated residue 12 side chain sticks
show sticks, kras and resi 12 and not name N+C+O+H
color yellow, kras and resi 12 and name SG
set stick_radius, 0.25

# 8. Focus camera onto protein
zoom kras and polymer, buffer=4
```

</details>

### 6.3. Generating Comparative Snapshot Grids
Render publication-quality structural snapshot grids comparing transitions over trajectories (such as WT vs Mutants across specific frames):

```bash
# Generate the default 3x8 transition snapshot grid (WT vs G12C vs G12D)
colabmda openmm snapshots
```

#### Customization via JSON Configuration:
```bash
colabmda openmm snapshots --config my_config.json --output figures/comparison_grid.png
```

```json
{
  "align_ref_pdb": "data/structure/4ldj/wt/target.B99990001_with_cryst.pdb",
  "stable_core_sel": "resi 1-10 or resi 40-55 or resi 80-169",
  "camera_view": [
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0,
    1.0, 0.0, 0.0,
    0.0, 0.0, -50.0
  ],
  "systems": {
    "WT": {
      "pdb": "data/simulation/4ldj_wt/r1/prod_full.pdb",
      "dcd": "data/simulation/4ldj_wt/r1/prod_full.dcd",
      "states": [1, 200],
      "times": ["0.00 ns", "2.00 ns"],
      "mut_residue": 12
    }
  }
}
```

---

## Acknowledgements
- [OpenMM](https://openmm.org) & [PDBFixer](https://github.com/openmm/pdbfixer)
- [Modeller](https://salilab.org/modeller/)
- [MDAnalysis](https://www.mdanalysis.org) & [MDTraj](https://www.mdtraj.org)
- [PDB2PQR & PROPKA](https://server.poissonboltzmann.org/)
- [NumPy](https://numpy.org), [Matplotlib](https://matplotlib.org), [Biopython](https://biopython.org)
- [Google Colab](https://colab.research.google.com) & [Miniforge/Conda](https://github.com/conda-forge/miniforge)

---

## Changelog

### v0.3.0 (Kinase-Substrate MD & Unified Architecture)
*   **Unified `cmda` CLI:** Single-keyword entrypoints (`stage`, `prep`, `sim`, `merge`, `analyze`, `compare`, `colab`) accessible via both `cmda` and `colabmda`.
*   **Kinase-Substrate Dynamics Engine:** Automated CDKL5 phosphotyrosine (`pY171` -> Amber `PTR` with `phosaa19SB`), quantum-parameterized nucleotides (Hu 2024 ATP/ADP/GDP), catalytic dual divalent ions ($Mg^{2+}$, $Mn^{2+}$), and continuous geometric substrate peptide docking (EB2 `KRPSSAPSG` at Ser222).
*   **High-Performance AmberTools Integration:** Automated Amber `tleap` parameterization with randomized ion neutralization (`addIonsRand`, solvating in <2s) and automated `cpptraj` PBC autoimaging and dual-stream solvent stripping (>90% storage savings).
*   **Strict Decoupled Architecture:** Guaranteed zero data files in `src/`, decoupled execution scripts in `<project>/scripts/` (`run.sh`, `submit.slurm`), and a single bootstrap script (`scripts/setup.sh`).
*   **Comprehensive Test Suite:** 40+ unit and regression tests passing across all modeling, parameterization, and simulation modules.

### v0.2.0 (Modular Refactor & Resilience)
*   **Single-Keyword Nomenclature:** Unified modular commands (`prep`, `stage`, `em`, `nvt`, `npt`, `qc`, `md`) and standardized `equil/` directory layout.
*   **PROPKA 3D pKa Titration:** Integrated PDBFixer and PDB2PQR/PROPKA 3 for physiological pH 7.4 protonation state assignment.
*   **Equilibration QC Gating:** Added automated `colabmda openmm qc` module generating publication-grade stability verification (`qc.png`).
*   **Decoupled Modules:** Isolated virtual screening workflows onto dedicated `docking` branch for lightweight performance.
*   **Automated Test Suite:** Added unit and regression test suite under `tests/` (`test_cli.py`, `test_forcefields.py`, `test_qc.py`).
*   **Force Field Resilience:** Integrated Amber14/ff19SB fallbacks and OPC 4-point water support.

### v0.1.0 (Initial Beta)
*   **Modular Pipeline:** New modular CLI for EM, NVT, NPT, and Production MD.
*   **Resume-Safe Engine:** Integrated checkpointing logic for fail-safe simulations on Google Colab.
*   **Modeling:** Automated Wild-Type building and mutation support via Modeller.
*   **Analysis:** Robust trajectory merging and comparative RMSD/Rg/RMSF analysis tools.
*   **Professional Standards:** Added CI/CD workflows, Black formatting, and Ruff linting.

---

## Citation

Initially, this repository was based on a GROMACS-on-Colab workflow, but we transformed it to utilize OpenMM to improve feasibility, performance, and stability for automated simulation recovery on cloud platforms.

If you use this tool, please consider citing the underlying studies:

> Paul SK, Saddam M, Rahaman KA, Choi JG, Lee SS, Hasan M. **Molecular modeling, molecular dynamics simulation, and essential dynamics analysis of grancalcin: An upregulated biomarker in experimental autoimmune encephalomyelitis mice.** *Heliyon*. 2022 Oct 23;8(10):e11232. doi: [10.1016/j.heliyon.2022.e11232](https://doi.org/10.1016/j.heliyon.2022.e11232). PMID: 36340004; PMCID: PMC9626934.
>
> [![Cited by](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fapi.openalex.org%2Fworks%2Fhttps%3A%2F%2Fdoi.org%2F10.1016%2Fj.heliyon.2022.e11232&query=%24.cited_by_count&label=cited%20by&color=blue)](https://doi.org/10.1016/j.heliyon.2022.e11232)
>
> Paul SK, Saddam M, Tabassum N, Hasan M. **Molecular dynamics simulation of wild and mutant proteasome subunit beta type 8 (PSMB8) protein: Implications for restoration of inflammation in experimental autoimmune encephalomyelitis pathogenesis.** *Heliyon*. 2024 Dec 15;11(1):e41166. doi: [10.1016/j.heliyon.2024.e41166](https://doi.org/10.1016/j.heliyon.2024.e41166). PMID: 39802026; PMCID: PMC11719297.
>
> [![Cited by](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fapi.openalex.org%2Fworks%2Fhttps%3A%2F%2Fdoi.org%2F10.1016%2Fj.heliyon.2024.e41166&query=%24.cited_by_count&label=cited%20by&color=blue)](https://doi.org/10.1016/j.heliyon.2024.e41166)

---

<p align="center">
  <img src="https://komarev.com/ghpvc/?username=paulshamrat&repo=ColabMDA&label=Visitors&color=blue&style=flat-square" alt="Visitors" />
</p>
