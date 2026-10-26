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
| **Dependencies** | `OpenMM`, `AmberTools`, `Modeller`, `MDAnalysis`, `MDTraj`, `PDB2PQR / PROPKA` |
| **Platform** | `Linux` · `Google Colab (CLI / Terminal)` · `HPC (SLURM)` |
| **Structure** | [Source Layout](#-project-structure) |

### 📂 Project Structure
ColabMDA is organized into clear, functional modules:
*   `src/colabmda/modeller/`: Homology modeling engine (Biological UniProt numbering supported).
*   `src/colabmda/openmm/`: OpenMM simulation engines (Modular prep/EM/NVT/NPT/QC/MD/analysis).
*   `src/colabmda/amber/`: AmberTools parameterization and tleap leap engines.
*   `envs/`: Automated installation scripts for scientific environments.
*   `scripts/`: Quick-start bootstrap scripts for Google Colab (`setup.sh`).
*   `tests/`: Unit and regression test suite (`pytest`).
*   `docs/`: Sphinx and ReadTheDocs scientific documentation source.

---

## 1. Quick Start (Google Colab Installation)

You can run ColabMDA directly inside the **Colab Terminal** or via the **Google Colab CLI** from your local machine.

### 1.1. Connect via Google Colab CLI (Recommended)

Run these commands in your local terminal:

```bash
# 1. Install Colab CLI & authenticate (one-time browser login)
pip install google-colab-cli
colab sessions

# 2. Launch Google Colab GPU session & mount Google Drive
colab new -s colabmda-example --gpu T4
# Visit the OAuth URL printed in terminal to grant Drive access, then press Enter:
colab drivemount --session colabmda-example

# 3. Connect to the interactive Colab console
colab console -s colabmda-example
```

*(Alternatively, in the Colab web interface: select `Runtime` → `Change runtime type` → **T4 GPU**, mount Drive, and open the terminal via the **⋮** menu → **Terminal**).*

---

### 1.2. Install ColabMDA & Activate Environment

Now **inside the Colab console / terminal**, install the environment and activate it:

```bash
# 1. Clone ColabMDA on the cloud GPU instance
git clone --branch dev1 https://github.com/paulshamrat/ColabMDA.git /content/ColabMDA
cd /content/ColabMDA

# 2. Run the automated installer with MODELLER support
WITH_MODELLER=1 bash scripts/setup.sh

# 3. Activate the scientific environment
source /root/miniforge3/etc/profile.d/conda.sh
conda activate colabmda

# 4. Verify installation
cmda --help
```

> [!IMPORTANT]
> **Modeller License Prompt:** During Step 2, the installer will ask you for your Modeller License Key. Academic licenses are completely free and sent instantly via email from [salilab.org/modeller/registration.html](https://salilab.org/modeller/registration.html).
> Paste your key and press **Enter** to finish setup.

> 💡 **Environment Activation:** All commands (`cmda` and `colabmda`) run in the single unified **`colabmda`** conda environment. If you ever open a new terminal or reconnect after a timeout, simply activate it:
> ```bash
> source /root/miniforge3/etc/profile.d/conda.sh
> conda activate colabmda
> ```

---

### 💡 Tip: How to Resume After a Timeout
If your Google Colab session expires:
1. Re-mount Drive and run **Step 1.2** to re-provision the environment.
2. Navigate to your project directory on Drive and run the **exact same `cmda sim md` command**.
3. ColabMDA automatically detects existing checkpoint files (`.chk`) and resumes without repeating completed work.

---

## 2. Manual / Detailed Installation (Advanced)

> ⚠️ **Note:** For standard Google Colab runs, use the **Quick Start (Section 1)** above. This section covers local workstations, HPC clusters, and test verification.

<details>
<summary><b>💻 A. Local Workstation Setup (Linux with NVIDIA GPU)</b></summary>

ColabMDA works on any Linux workstation with an NVIDIA GPU using the provided `environment.yml`:
```bash
mamba env create -f environment.yml
conda activate colabmda
pip install -e .
```
</details>

<details>
<summary><b>🧪 B. Running the Test Suite</b></summary>

Verify code integrity, force field mechanics, and CLI parsing:
```bash
pytest tests/ -v
```
</details>

<details>
<summary><b>🏢 C. HPC Usage (SLURM)</b></summary>

ColabMDA scripts can be integrated directly into SLURM job submission arrays. Because chunked simulations write discrete `.chk` checkpoints, jobs finish cleanly within standard cluster walltime limits.
</details>

---

## 3. Recommended Directory Strategy

ColabMDA uses a clean, reproducible directory hierarchy stored directly in your Google Drive (or local project directory). Biochemical states (**`apo`** and **`atpmg`**) are organized into parallel branches under `data/`:

```text
research/projects/colabmda-example/       # Your project root on Google Drive
└── data/
    ├── apo/                              # Track 1: Apo state (protein in water)
    │   ├── structure/
    │   │   ├── raw/                      # Downloaded template (e.g. 4BGQ) and UniProt FASTA
    │   │   ├── mapping/                  # Residue mapping and missing loop registers
    │   │   └── models/
    │   │       ├── WT/                   # Modeled WT coordinates (WT.pdb)
    │   │       └── L119R/                # Modeled L119R coordinates (L119R.pdb)
    │   │
    │   ├── simulation/
    │   │   ├── WT/
    │   │   │   ├── 01_prep/              # Solvated Amber system (prmtop / inpcrd / pdb)
    │   │   │   ├── 02_eq/                # 24-stage equilibration hub + qc.png
    │   │   │   ├── r1/                   # Production replica 1 (chunks + full trajectory)
    │   │   │   ├── r2/                   # Production replica 2 (independent seed from 02_eq)
    │   │   │   └── r3/                   # Production replica 3 (triplicate simulation)
    │   │   └── L119R/
    │   │       ├── 01_prep/
    │   │       ├── 02_eq/
    │   │       └── r1/, r2/, r3/
    │   │
    │   └── analysis/
    │       ├── WT/                       # Single-replica reports (RMSD, Rg, RMSF)
    │       ├── L119R/
    │       └── compare/figures/          # Triplicate Mean ± SD overlay plots
    │
    └── atpmg/                            # Track 2: Holo state (ATP + Mg²⁺ complex)
        ├── structure/models/             # Assembled complexes (WT and L119R with ATP + Mg²⁺)
        ├── simulation/                   # WT/ and L119R/ (01_prep, 02_eq, r1, r2, r3)
        └── analysis/                     # WT/, L119R/, and compare/figures/
```

> 💡 **Why Shared Equilibration (`02_eq/`)?**
> The 24-stage equilibration ladder (EM, thermalization, density relaxation) only needs to be calculated once per variant. Replicas (`r1`, `r2`, `r3`) fork from `02_eq/` with distinct random velocity seeds, saving substantial GPU time while guaranteeing true statistical independence.

---

## 4. Simulation Workflow

### Pipeline at a Glance:
1. **Model:** Fetch template coordinates, map UniProt numbering, and build WT and mutant models.
2. **Titrate & Prep:** Assign physiological protonation states with PROPKA and solvate with Amber ff19SB/OPC.
3. **Stage:** Initialize simulation folders under `data/simulation/`.
4. **Equilibrate:** Run 24-stage equilibration in `02_eq/` with progressive restraint release.
5. **QC Gate:** Verify equilibration density, temperature, and stability (`qc.png`).
6. **Production MD:** Execute resume-safe chunked production MD across triplicates (`r1`, `r2`, `r3`).
7. **Merge:** Combine trajectory chunks into centered, wrapped, solvent-stripped trajectories.
8. **Analyze & Compare:** Compute RMSD, $R_g$, and RMSF, and overlay WT vs Mutant triplicate statistics.

---

### Step 0: Create Project Directory & Activate Environment

Before running workflow commands, switch to your Google Drive project workspace and ensure your environment is active:

```bash
# 1. Create and navigate to your project directory
mkdir -p /content/drive/MyDrive/research/projects/colabmda-example
cd /content/drive/MyDrive/research/projects/colabmda-example

# 2. Activate the ColabMDA environment
source /root/miniforge3/etc/profile.d/conda.sh
conda activate colabmda
```

---

### 4.1. Homology Modeling & Point Mutations (Modeller)
**Environment:** `colabmda`

The modeling workflow enforces strict **Biological UniProt Numbering**, automated loop completion, and coordinate validation gates:

```bash
# 1. Fetch template coordinates (PDB 4BGQ) and map UniProt (O76039, residues 1–303)
cmda fetch --state apo --pdb-id 4BGQ --uniprot-id O76039 --range 1 303

# 2. Build full Wild-Type model (fills all 27 missing loop residues)
cmda model --state apo --variant WT

# 3. Build Point Mutant (L119R) on the verified WT backbone
cmda model --state apo --variant L119R
```

*What happens behind the scenes:*
* **`cmda fetch`**: Downloads experimental coordinates (mmCIF/PDB fallback), aligns them to canonical UniProt sequence, identifies missing loop segments, and writes validation manifests (`mapping/input.json`).
* **`cmda model WT`**: Invokes Modeller to build full coordinates for residues 1–303, closing missing loops and preserving author residue numbering.
* **`cmda model L119R`**: Introduces the point mutation onto the verified WT template and optimizes local sidechain rotamers.
* **Pre-built / AlphaFold models:** If you already have pre-built coordinates, verify them directly via `cmda check --state apo --variant WT --model /path/to/model.pdb`.

---

### 4.2. Structure Preparation & 3D pKa Titration (PROPKA at pH 7.4)
**Environment:** `colabmda`

Accurate protonation states are critical for biophysical fidelity. ColabMDA integrates **PDBFixer**, **PROPKA 3**, and **AmberTools `tleap`** to assign physiological protonation states at pH 7.4, solvate with OPC 4-point water, neutralize with 150 mM NaCl, and build Amber `ff19SB` topologies:

```bash
# 1. Parameterize and solvate Wild-Type
cmda prep --state apo --variant WT

# 2. Parameterize and solvate Mutant (L119R)
cmda prep --state apo --variant L119R
```

*Key outputs in `data/apo/simulation/<variant>/01_prep/`:*
* `system.prmtop`: Amber topology file containing all atomic parameters.
* `system.inpcrd`: Solvated, charge-neutralized input coordinates.
* `system.pdb`: Reference starting coordinate file.

---

### 4.3. Simulation Execution: 24-Stage Equilibration & Production MD
**Environment:** `colabmda`

ColabMDA uses a standardized 24-stage equilibration ladder (1.96 ns total) that progressively heats the system to 300 K and relaxes positional restraints before production MD:

#### Step 1: Initialize Stage & Shared Equilibration Hub (`02_eq/`)
```bash
# Initialize simulation directories
cmda stage --state apo --variant WT

# Run full 24-stage equilibration (EM, NVT heating, NPT density relaxation)
cmda sim eq --state apo --variant WT

# Check equilibration stability QC gate (generates qc.png)
cmda sim qc --state apo --variant WT
```

*(You can also run individual modular stages if preferred: `cmda sim em`, `cmda sim nvt`, `cmda sim npt`).*

#### Step 2: Production MD Across Triplicates (`r1`, `r2`, `r3`)
Replicas fork from `02_eq/` with distinct random velocity seeds:

```bash
# Replica 1 (e.g. 10 ns chunked production run)
cmda sim md --state apo --variant WT --replica r1 --length-ns 10

# Replica 2
cmda sim md --state apo --variant WT --replica r2 --length-ns 10

# Replica 3
cmda sim md --state apo --variant WT --replica r3 --length-ns 10
```

*(Repeat the same commands for mutant: `cmda stage --state apo --variant L119R`, `cmda sim eq ...`, and `cmda sim md ...`).*

> 🛡️ **Resume Safety:** Production MD saves checkpoints (`prod.chk`) at regular intervals. If Colab disconnects, simply re-run the same command and ColabMDA will resume from the exact last completed chunk.

---

### 4.4. Merge Trajectories & Apply PBC Autoimaging
**Environment:** `colabmda`

Concatenate trajectory chunks, apply periodic boundary condition (PBC) autoimaging, center the protein, and strip solvent:

```bash
# Merge Wild-Type replicas
cmda merge --state apo --variant WT

# Merge Mutant replicas
cmda merge --state apo --variant L119R
```

*Key outputs in `data/apo/simulation/<variant>/r1/`:*
* `prod_full.dcd`: The final, concatenated, and centered/wrapped trajectory file.
* `prod_full.pdb`: Matching reference topology PDB file.

---

### 4.5. Trajectory Analysis & Direct Overlay Comparison
**Environment:** `colabmda`

#### Step 1: Single-Variant Analysis (RMSD, Rg, RMSF)
Calculate backbone metrics for each replica:

```bash
# Analyze Wild-Type replicas
cmda analysis --state apo --variant WT

# Analyze Mutant replicas
cmda analysis --state apo --variant L119R
```

#### Step 2: Direct Comparative Overlay (WT vs Mutant)
Automatically aggregate triplicate statistics and generate publication-grade comparison plots:

```bash
# Compare WT vs L119R (computes Triplicate Mean ± SD)
cmda compare --state apo --variants WT,L119R
```

*Publication-grade figures generated in `data/apo/analysis/compare/figures/`:*
* `rmsd_avg.png`: Backbone RMSD overlay across time (Mean ± SD band).
* `rg_avg.png`: Radius of gyration ($R_g$) compaction profiles.
* `rmsf_avg.png`: Per-residue flexibility plot highlighting dynamic differences.

---

### 4.6. Track 2: Holo Dynamics (ATP + $\text{Mg}^{2+}$ Complex)

To simulate the catalytically active holo kinase complex with ATP and $\text{Mg}^{2+}$ bound in the active site:

#### 1. Assemble Holo Complex from Donor Template (PDB 8FP5)
Transfer ATP and $\text{Mg}^{2+}$ from donor structure [8FP5](https://www.rcsb.org/structure/8FP5) onto your verified WT and mutant models:

```bash
# Assemble Holo WT
cmda assemble --state atpmg --template 8FP5 --variant WT --relax-sidechains

# Assemble Holo Mutant (L119R)
cmda assemble --state atpmg --template 8FP5 --variant L119R --relax-sidechains
```

#### 2. Execute Holo Workflow
Every subsequent command is identical to the apo workflow, simply specifying `--state atpmg`:

```bash
# Prep & Solvate
cmda prep --state atpmg --variant WT
cmda prep --state atpmg --variant L119R

# Stage & Equilibrate
cmda stage --state atpmg --variant WT
cmda sim eq --state atpmg --variant WT

# Production MD
cmda sim md --state atpmg --variant WT --replica r1 --length-ns 10

# Merge & Compare
cmda merge --state atpmg --variant WT
cmda merge --state atpmg --variant L119R
cmda compare --state atpmg --variants WT,L119R
```

---

## 5. Trajectory Visualization in PyMOL

Visualize molecular dynamics trajectories in PyMOL with full secondary structure (ribbon/cartoon) and active site sidechains:

### 5.1. Local Setup
On your local computer or laptop (with PyMOL installed):
```bash
# Install ColabMDA lightweight CLI locally
pip install "git+https://github.com/paulshamrat/ColabMDA.git@dev1"
```

### 5.2. Launch Trajectory View
```bash
# View Wild-Type trajectory frame
cmda view --state apo --variant WT --replica r1

# Generate 3x8 transition snapshot grids comparing WT vs Mutant
cmda snapshots --state apo
```

---

## Acknowledgements
- [OpenMM](https://openmm.org) & [PDBFixer](https://github.com/openmm/pdbfixer)
- [Modeller](https://salilab.org/modeller/)
- [AmberTools](https://ambermd.org/AmberTools.php)
- [MDAnalysis](https://www.mdanalysis.org) & [MDTraj](https://www.mdtraj.org)
- [PDB2PQR & PROPKA](https://server.poissonboltzmann.org/)
- [NumPy](https://numpy.org), [Matplotlib](https://matplotlib.org), [Biopython](https://biopython.org)
- [Google Colab](https://colab.research.google.com) & [Miniforge/Conda](https://github.com/conda-forge/miniforge)

---

## Changelog

### v0.2.0 (Kinase Architectures, Amber ff19SB & Resilience)
*   **Unified CLI (`cmda`):** Intuitive top-level commands (`fetch`, `model`, `assemble`, `prep`, `stage`, `sim`, `merge`, `analysis`, `compare`).
*   **Dual State Architecture:** Clean separation between unliganded **`apo`** dynamics and catalytic **`atpmg`** holo dynamics.
*   **Automated Triplicate Comparison:** Single-call `cmda compare --state <state> --variants WT,MUT` generating publication-grade Mean ± SD overlay plots.
*   **24-Stage Equilibration Hub:** Robust heating and restraint release ladder guaranteeing stable density and temperature.
*   **PROPKA 3D pKa Titration & AmberTools `tleap`:** Physiological pH 7.4 protonation assignment with Amber ff19SB and OPC 4-point water.
*   **Resume-Safe Checkpointing:** Automated recovery from Google Colab session timeouts.

### v0.1.0 (Initial Beta)
*   Modular pipeline for EM, NVT, NPT, and Production MD.
*   Checkpointing logic for fail-safe simulations on Google Colab.
*   Wild-Type building and mutation support via Modeller.

---

## Citation

If you use ColabMDA in your research, please cite the underlying studies:

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
