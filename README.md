# ColabMDA

[![Release](https://img.shields.io/github/v/tag/paulshamrat/ColabMDA)](https://github.com/paulshamrat/ColabMDA/tags) [![Python CI](https://github.com/paulshamrat/ColabMDA/actions/workflows/python-test.yml/badge.svg?branch=main)](https://github.com/paulshamrat/ColabMDA/actions/workflows/python-test.yml) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![GitHub](https://img.shields.io/badge/GitHub-Repository-blue?logo=github)](https://github.com/paulshamrat/ColabMDA) [![Issues](https://img.shields.io/github/issues/paulshamrat/ColabMDA)](https://github.com/paulshamrat/ColabMDA/issues) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black) [![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

**ColabMDA** is a specialized tool that lets you run high-quality Molecular Dynamics simulations on Google Colab without the fear of losing your work. The biggest problem with Colab is that it disconnects, often destroying hours of simulation data. **ColabMDA fixes this** with a "resume-safe" system that automatically saves your progress to Google Drive. If your session expires, you can resume exactly where you left off with one simple command. From modeling protein mutations to generating publication-ready analysis, ColabMDA handles the complex setup for you, making it the easiest way to get high-quality MD results using free cloud GPUs.

---

## 1. Quick Start (Google Colab CLI)

ColabMDA is designed to be operated directly from your terminal using the official **Google Colab CLI (`google-colab-cli`)**. This eliminates browser disconnects, clunky web editors, and lost tabs by streaming an interactive GPU terminal session directly into your shell.

### 1.1. Launch GPU Session & Mount Drive (Local Terminal)

Run these three commands in your local terminal:

```bash
# 1. Install Colab CLI & authenticate (one-time browser login)
pip install google-colab-cli
colab sessions

# 2. Launch Google Colab GPU session & mount Google Drive
colab new -s colabmda-example --gpu T4
# Visit the OAuth link printed in terminal to approve Drive access, then press Enter:
colab drivemount --session colabmda-example

# 3. Connect to the interactive Colab console terminal
colab console -s colabmda-example
```

> 💡 **Why Colab CLI?**
> * **No Browser Tab Needed:** You interact with the cloud GPU directly in your native terminal.
> * **Instant Reconnect:** If your connection drops, simply run `colab console -s colabmda-example` to reattach to the running instance.
> * **Full Terminal Power:** Native keyboard shortcuts, tmux/screen, and bash history without browser conflicts.
>
> *(Optional Browser Fallback: If you prefer the Google Colab web interface, select `Runtime` → `Change runtime type` → **T4 GPU**, mount Drive, and click the **⋮** menu → **Terminal**).*

---

### 1.2. Install ColabMDA & Activate Environment (Inside Colab Console)

Now that you are attached to the Colab GPU terminal via `colab console`, install the scientific environment:

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
1. Launch a session and re-mount Drive:
   ```bash
   colab new -s colabmda-example --gpu T4
   colab drivemount --session colabmda-example
   colab console -s colabmda-example
   ```
2. Re-run **Step 1.2** to re-provision the environment.
3. Navigate to your project directory on Drive and run the **exact same `cmda sim md` command**.
4. ColabMDA automatically detects existing checkpoint files (`.chk`) and resumes without repeating completed work.

---

## 2. Local Workstation Setup & Hybrid Cloud Architecture

ColabMDA supports a high-efficiency **Hybrid Workflow**:
1. **Local Workstation / Laptop (CPU):** Run fast homology modeling with Modeller (~20 seconds) and perform interactive 3D visual inspection in local PyMOL.
2. **Google Drive Sync:** Push your lightweight modeled structures (`data/`) to Google Drive using `rclone` or desktop sync.
3. **Google Colab (GPU):** Mount Drive and execute heavy 24-stage equilibration and 100 ns OpenMM MD on free cloud GPUs (NVIDIA T4/A100).

---

### 2.1. Local Installation (Modeller + PyMOL)

On your local laptop or workstation (macOS, Linux, or Windows WSL):

```bash
# 1. Create a local environment
conda create -n colabmda-local python=3.12 -y
conda activate colabmda-local

# 2. Install MODELLER (CPU homology modeling)
conda install -c salilab modeller -y
# When prompted or in ~/.bashrc, set your free academic license:
# export KEY_MODELLER="YOUR_ACADEMIC_KEY"

# 3. Install Open-Source PyMOL (Free interactive 3D visualization)
conda install -c conda-forge pymol-open-source -y

# 4. Install ColabMDA CLI
pip install git+https://github.com/paulshamrat/ColabMDA.git@dev1
```

*(Alternative PyMOL install: Ubuntu/Debian `sudo apt install pymol`, macOS `brew install brewsci/bio/pymol`).*

---

### 2.2. Transferring Local Data to Google Drive

After generating and inspecting your models locally, transfer the lightweight `data/` folder to Google Drive:

#### Option A: Terminal CLI via `rclone` (Recommended)
```bash
# Push local data to Google Drive project workspace:
rclone copy data/ gdrive:research/projects/colabmda-example/data/ \
  --exclude ".git/**" --exclude ".~lock.*#" --fast-list -P
```

#### Option B: Google Drive for Desktop
Work directly inside or symlink into your Google Drive mount (`~/Google Drive/MyDrive/research/projects/colabmda-example/`). Any file generated locally is automatically uploaded in the background.

---

### 2.3. Advanced Environments & HPC

<details>
<summary><b>🧪 Running the Test Suite Locally</b></summary>

Verify code integrity, force field mechanics, and CLI parsing:
```bash
pytest tests/ -v
```
</details>

<details>
<summary><b>🏢 HPC Usage (SLURM)</b></summary>

ColabMDA scripts can be integrated directly into SLURM job submission arrays. Because chunked simulations write discrete `.chk` checkpoints, jobs finish cleanly within standard cluster walltime limits.
</details>

---

## 3. Recommended Directory Strategy

ColabMDA uses a clean, reproducible directory hierarchy stored directly in your Google Drive (or local project directory). Biochemical states (**`apo`** and **`atpmg`**) are organized into parallel branches under `data/`:

```text
research/projects/colabmda-example/       # Project root (or custom project folder)
└── data/
    ├── apo/                              # Track 1: Apo state (protein in water)
    │   ├── structure/
    │   │   └── O76039_4BGQ/              # Biological Target (UniProt) + Structural Template (PDB)
    │   │       ├── raw/                  # Downloaded archive (4bgq.cif), template.pdb, uniprot.fasta
    │   │       ├── mapping/              # Residue mapping (residues.json) and manifest (input.json)
    │   │       ├── models/
    │   │       │   ├── WT/               # Modeled WT coordinates (WT.pdb)
    │   │       │   └── L119R/            # Modeled L119R coordinates (L119R.pdb)
    │   │       └── STRUCTURE_DONE        # Verification completion marker
    │   │
    │   ├── simulation/
    │   │   └── run01/                    # Run session (defaults automatically to run01)
    │   │       ├── WT/
    │   │       │   ├── 01_prep/          # Solvated Amber system (model.pdb, system.prmtop, system.inpcrd)
    │   │       │   ├── 02_eq/            # 24-stage equilibration hub + qc.png
    │   │       │   └── 03_md/            # Production MD replicas
    │   │       │       ├── r1/           # Replica 1 (chunks + full trajectory)
    │   │       │       ├── r2/           # Replica 2 (independent seed from 02_eq)
    │   │       │       └── r3/           # Replica 3 (triplicate simulation)
    │   │       └── L119R/
    │   │           ├── 01_prep/
    │   │           ├── 02_eq/
    │   │           └── 03_md/{r1, r2, r3}/
    │   │
    │   └── analysis/
    │       └── run01/
    │           ├── WT/
    │           │   ├── r1/               # Single-replica reports (RMSD, Rg, RMSF)
    │           │   ├── r2/
    │           │   └── r3/
    │           ├── L119R/
    │           │   └── {r1, r2, r3}/
    │           └── compare/figures/      # Triplicate Mean ± SD overlay plots
    │
    └── atpmg/                            # Track 2: Holo state (ATP + Mg²⁺ complex)
        ├── structure/
        │   └── O76039_4BGQ/
        │       ├── donors/
        │       │   └── 8FP5/             # Isolated cofactor donor template folder
        │       │       └── 8fp5.cif      # Donor crystal archive containing ATP and Mg²⁺
        │       ├── raw/                  # Inherited target archive and mapping
        │       ├── models/
        │       │   ├── WT/
        │       │   │   ├── complex.pdb   # Assembled Protein + ATP + Mg²⁺ complex
        │       │   │   ├── complex.json  # Core RMSD & ATP/Mg fit audit trail
        │       │   │   └── ASSEMBLY_DONE
        │       │   └── L119R/
        │       │       └── complex.pdb
        │       └── STRUCTURE_DONE
        │
        ├── simulation/run01/{WT, L119R}/{01_prep, 02_eq, 03_md/{r1, r2, r3}}
        └── analysis/run01/{WT, L119R, compare/figures}
```

> 💡 **Why `UniProt_PDB` Folders (`O76039_4BGQ/`)?**
> Organizing structures by `<uniprot>_<pdb>` keeps coordinates, mappings, and models completely separated when comparing alternative templates (e.g. 4BGQ vs 8FP5) or modeling multiple proteins in the same project workspace.

> 💡 **Why Separate `donors/<template>/` Folders?**
> In holo simulations, the cofactor donor crystal structure (e.g. 8FP5) is kept in its own dedicated `donors/8FP5/` subfolder. This completely isolates the donor coordinates from the target protein files (`raw/`), preventing any confusion about which PDB provided the backbone vs which provided the ligands.

> 💡 **Why `model.pdb` in `01_prep/`?**
> When you run `cmda prep`, ColabMDA retrieves the modeled variant (`WT.pdb` or `complex.pdb`) and copies it into `01_prep/model.pdb`. This creates a self-contained, immutable starting snapshot for that specific simulation run (`run01`), guaranteeing full reproducibility even if upstream modeling files are modified later.

> 💡 **Why Shared Equilibration (`02_eq/`)?**
> The 24-stage equilibration ladder (EM, thermalization, density relaxation) only needs to be calculated once per variant. Replicas (`r1`, `r2`, `r3`) fork from `02_eq/` with distinct random velocity seeds, saving substantial GPU time while guaranteeing true statistical independence.

---

## 4. Simulation Workflow: Step-by-Step Execution

ColabMDA divides work logically into three clear phases:
* **Phase 1 (Local):** Fast CPU homology modeling and interactive PyMOL QA on your workstation (~30 seconds).
* **Phase 2 (Sync):** Pushing verified coordinates to Google Drive via `rclone copy` (~5 MB).
* **Phase 3 (Cloud):** Heavy GPU parameterization, 24-stage equilibration, and 100 ns OpenMM MD on Google Colab GPU.

---

### 💻 Phase 1: Local Modeling & PyMOL Visual Inspection (Local Workstation)
*Execute this on your local laptop or workstation (macOS, Linux, or Windows WSL). No GPU is needed.*

#### 4.1. Initialize Local Project Workspace
```bash
# 1. Create and navigate to your local project directory
# (Replace 'colabmda-example' with your custom project name if desired, e.g. cdkl5-dynamics)
mkdir -p ~/research/projects/colabmda-example
cd ~/research/projects/colabmda-example

# 2. Activate your local environment
conda activate colabmda-local
```

#### 4.2. Homology Modeling & Point Mutations (Modeller)
The modeling workflow enforces strict **Biological UniProt Numbering**, automated loop completion, and coordinate validation gates:

```bash
# 1. Fetch template coordinates (PDB 4BGQ) and map UniProt (O76039, residues 1–303)
cmda fetch --state apo --pdb-id 4BGQ --uniprot-id O76039 --range 1 303 --mutations L119R

# 2. Build full Wild-Type model (fills all 27 missing loop residues)
cmda model --state apo --variant WT

# 3. Build Point Mutant (L119R) on the verified WT backbone
cmda model --state apo --variant L119R
```

*What happens behind the scenes:*
* **`cmda fetch` Output Banner**: Displays experimental resolution, coordinate density, missing loop segments, and template sequence checks:
  ```text
  ==============================================================================
   ColabMDA Structural Mapping & Verification
  ==============================================================================
   UniProt Target:   O76039 (Residues 1–303, Total Length: 960 AA)
   Template PDB:     4BGQ (Chain A, Method: X-RAY DIFFRACTION, Resolution: 2.00 Å)
   Coordinate Gate:  276/303 observed residues resolved (91.1% coordinate density)
   Missing Loops:    27 missing residues identified across 4 loop segment(s):
                     - Loop 1: Res 1–8 (8 residues)
                     - Loop 2: Res 46–52 (7 residues)
                     - Loop 3: Res 163–173 (11 residues)
                     - Loop 4: Res 303–303 (1 residues)
   Template Check:   ✓ All 276 observed residues match Canonical UniProt WT
   Mutation Gate:    ✓ Target variant L119R verified on WT residue L119 (coordinate & UniProt match)
   Outputs Saved:    data/apo/structure/O76039_4BGQ/
  ==============================================================================
  ```
* **Engineered Template Mutations**: If crystallographers introduced mutations into the PDB template, `cmda fetch` reports them, and Modeller automatically reverts them to canonical UniProt WT during `cmda model WT`.
* **`cmda model WT`**: Invokes Modeller to build full coordinates for residues 1–303, closing missing loops and preserving author residue numbering.
* **`cmda model L119R`**: Introduces the point mutation onto the verified WT template and optimizes local sidechain rotamers.

#### 4.3. (Optional) Holo ATP/$\text{Mg}^{2+}$ Coordinate Assembly (Track 2)
To investigate catalytic active-site geometry and cofactor coordination, assemble ATP and $\text{Mg}^{2+}$ from donor PDB 8FP5 onto your models:

```bash
# 1. Assemble Holo Wild-Type from CDK2 donor (PDB 8FP5)
cmda assemble --state atpmg --variant WT --template 8FP5 --relax-sidechains

# 2. Assemble Holo Mutant (L119R) from the same donor pose
cmda assemble --state atpmg --variant L119R --template 8FP5 --relax-sidechains
```
*What happens:* Downloads `8fp5.cif` into `data/atpmg/structure/O76039_4BGQ/donors/8FP5/`, aligns kinase catalytic core residues (30–220), transfers ATP and $\text{Mg}^{2+}$, relaxes pocket sidechains, and writes `complex.pdb`.

#### 4.4. Visual Quality Assurance in Local PyMOL
Verify loop closure and mutated residue rotamers interactively:

```bash
# Inspect Wild-Type model in PyMOL:
cmda view --state apo --variant WT

# Or compare WT and Mutant side-by-side:
pymol data/apo/structure/O76039_4BGQ/models/WT/WT.pdb data/apo/structure/O76039_4BGQ/models/L119R/L119R.pdb
```
*Visual Check:* Confirm loop 163–173 is closed smoothly and Arg119 points into solvent with zero steric clashes.

---

### 🚀 Phase 2: Push Modeled Structures to Google Drive
*Once models are verified locally, transfer the lightweight `data/` directory (~5 MB) to Google Drive so Google Colab's GPU can access them.*

#### Option A: Terminal CLI via `rclone` (Recommended)
From your local project folder:
```bash
rclone copy data/ gdrive:research/projects/colabmda-example/data/ \
  --exclude ".git/**" --exclude ".~lock.*#" --fast-list -P
```

#### Option B: Google Drive for Desktop
If you worked inside your Google Drive desktop mount (`~/Google Drive/MyDrive/research/projects/colabmda-example/`), all files are already synced automatically in the background.

---

### ⚡ Phase 3: GPU Simulation & Analysis (Google Colab Cloud)
*Now connect to Google Colab to run heavy Amber parameterization and OpenMM molecular dynamics on an NVIDIA T4 / A100 GPU.*

#### 4.5. Connect to Colab GPU & Activate Environment
In your local terminal:
```bash
# Launch Colab session with T4 GPU and mount Google Drive
colab new -s colabmda-example --gpu T4
colab drivemount --session colabmda-example
colab console -s colabmda-example
```

Inside the interactive Colab console terminal:
```bash
# 1. Navigate to the synced project workspace on Google Drive
cd /content/drive/MyDrive/research/projects/colabmda-example

# 2. Activate the ColabMDA environment
source /root/miniforge3/etc/profile.d/conda.sh
conda activate colabmda
```

#### 4.6. Structure Preparation & 3D pKa Titration (PROPKA at pH 7.4)
Assign physiological protonation states at pH 7.4, solvate with OPC water, neutralize with 150 mM NaCl, and build Amber `ff19SB` topologies:

```bash
# Track 1: Apo state
cmda prep --state apo --variant WT
cmda prep --state apo --variant L119R

# Track 2: Holo state (includes Hu et al. B3 ATP and Li–Merz 12–6 Mg²⁺)
cmda prep --state atpmg --variant WT
cmda prep --state atpmg --variant L119R
```

*Key outputs in `data/<state>/simulation/run01/<variant>/01_prep/`:*
* `model.pdb`: Immutable copy of the input coordinate PDB for complete reproducibility.
* `system.prmtop`: Amber topology containing all physical parameters.
* `system.inpcrd`: Solvated, charge-neutralized initial coordinates.

#### 4.7. Simulation Execution: 24-Stage Equilibration & QC Gate
Run the 24-stage progressive heating and restraint release ladder (EM, NVT heating, NPT density relaxation) on CUDA GPU:

```bash
# 1. Initialize simulation directories
cmda stage --state apo --variant WT

# 2. Run full 24-stage equilibration
cmda sim eq --state apo --variant WT

# 3. Check equilibration stability QC gate (generates qc.png)
cmda sim qc --state apo --variant WT
```
*(For Holo state, run with `--state atpmg`).*

#### 4.8. Production MD Across Triplicates (`03_md/{r1, r2, r3}`)
Replicas fork from `02_eq/` with distinct random velocity seeds into `03_md/`:

```bash
# Replica 1 (100 ns production run, or 10 ns for testing)
cmda sim md --state apo --variant WT --replica r1 --total-ns 100

# Replica 2
cmda sim md --state apo --variant WT --replica r2 --total-ns 100

# Replica 3
cmda sim md --state apo --variant WT --replica r3 --total-ns 100
```
*(Repeat the same commands for mutant: `cmda stage --state apo --variant L119R`, `cmda sim eq ...`, and `cmda sim md ...`).*

> 🛡️ **Resume Safety:** Production MD saves checkpoints (`prod.chk`) at regular intervals directly to Google Drive. If Colab disconnects, simply re-run the exact same command to resume immediately with zero lost progress.

#### 4.9. Merge Trajectories & Autoimage PBC
Concatenate trajectory chunks, apply periodic boundary condition (PBC) autoimaging, center the protein, and strip solvent:

```bash
# Merge Wild-Type replicas
cmda merge --state apo --variant WT

# Merge Mutant replicas
cmda merge --state apo --variant L119R
```

#### 4.10. Trajectory Analysis & Direct Overlay Comparison
Calculate backbone metrics and overlay WT vs Mutant triplicate statistics:

```bash
# 1. Single-variant analysis (RMSD, Radius of Gyration, RMSF)
cmda analyze --state apo --variant WT
cmda analyze --state apo --variant L119R

# 2. Direct comparative overlay (WT vs L119R with Mean ± SD band)
cmda compare --state apo --variants WT,L119R
```

*Publication-grade figures generated in `data/apo/analysis/run01/compare/figures/`:*
* `rmsd_avg.png`: Backbone RMSD overlay across time (Mean ± SD band).
* `rg_avg.png`: Radius of gyration ($R_g$) compaction profiles.
* `rmsf_avg.png`: Per-residue flexibility plot highlighting dynamic differences.

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
