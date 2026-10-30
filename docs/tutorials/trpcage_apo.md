# Trp-Cage (1L2Y) Apo State: 5-Minute Quickstart & Smoke Test

This tutorial provides a lightweight, 5-minute simulation walkthrough using the 20-residue **Trp-cage miniprotein** (PDB **1L2Y**, NMR Model 1, residues **1–20**). It serves as a rapid "smoke test" to verify that your local environment, Google Drive synchronization, Colab GPU runtime, and simulation checkpointing are fully operational before launching full research production runs.

---

## 📁 Trp-Cage Directory Architecture Map

ColabMDA organizes the Trp-cage project files within a deterministic hierarchy:

```text
research/projects/colabmda-example/
├── scripts/                      # Fast runner scripts generated via staging
│   ├── 01_structure/01_fetch.py
│   └── 02_simulation/
│       ├── 01_prepare.sh
│       ├── 02_equilibrate.sh
│       ├── 03_production.sh
│       ├── 04_analysis.sh
│       └── 05_check.py
└── data/
    ├── 00_setup/example.json     # Staging metadata & integrity checksums
    ├── 01_structure/trpcage/     # 1L2Y.cif, WT.pdb (20 residues), input.json
    └── 02_simulation/<run>/
        ├── provenance.json       # Environment, GPU, and package manifest
        └── apo/WT/
            ├── 01_prep/          # system.prmtop, system.inpcrd, solvated.pdb
            ├── 02_eq/            # system.xml, npt.state.xml, equilibration_qc.json
            ├── 03_md/r1/         # prod_keep_ligand.dcd, prod_ca.dcd, prod.chk
            └── 04_analysis/      # rmsd.csv, rg.csv, rmsf.csv, figures/
```

---

## 💻 Phase 1: Local Setup & Structure Preparation (Local Machine)
*Execute this phase on your local laptop or workstation CPU. No GPU compute units are required.*

### 1.0. Local Machine Setup (One-Time)
Run this one-time setup in your terminal (Linux or Windows WSL):

```bash
# 1. Install cmda, OpenMM, AmberTools, PyMOL, and Modeller:
bash -c "$(curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/dev1/scripts/setup.sh)"

# 2. Activate environment:
source ~/miniforge3/etc/profile.d/conda.sh
conda activate cmda

# 3. Log in to Google Colab (one-time browser login):
colab sessions
```

### 1.1. Project Workspace & Example Staging
Create your project directory and stage the Trp-cage example:

```bash
mkdir -p ~/research/projects/colabmda-example
cd ~/research/projects/colabmda-example

# Stage Trp-cage runner scripts and verified baseline hashes:
cmda stage --example trpcage
```

### 1.2. Fetch Canonical 1L2Y Coordinates
Fetch the 20-residue NMR structure (Model 1, Chain A of 1L2Y):

```bash
python scripts/01_structure/01_fetch.py
```

*This verifies all 20 Cα atoms, confirms the deposited sequence (`NLYIQWLKDGGPSSGRPPPS`), and saves `data/01_structure/trpcage/WT.pdb`.*

### 1.3. Structure Preparation & Solvation
Solvate the miniprotein in OPC water and neutralize with counterions:

```bash
export COLABMDA_RUN_ID=run01
bash scripts/02_simulation/01_prepare.sh
```

*Key outputs generated in `data/02_simulation/run01/apo/WT/01_prep/`:*
* `model.pdb`: Input coordinate PDB.
* `solvated.pdb`: Solvated coordinate PDB with protein, water box, and counterion.
* `system.prmtop`: Amber topology containing all physical parameters.
* `system.inpcrd`: Solvated, charge-neutralized initial coordinates.

---

## 🚀 Phase 2: Launch Colab GPU & Stream Data to Drive
*Launch your Colab GPU instance and stream the prepared project directory to Google Drive.*

In your local terminal:

```bash
# 1. Launch a Colab GPU instance:
colab new -s cmda --gpu T4

# 2. Mount Google Drive:
colab drivemount --session cmda

# 3. Stream local data/ directly to Google Drive:
cmda push -s cmda
```

> [!TIP]
> **Data Destination on Google Drive:**
> * **Google Drive Filesystem Path:** `/content/drive/MyDrive/research/projects/colabmda-example/data/`
> * **Mechanism:** `cmda push` packages your local `data/` folder, streams it to your active Colab VM (`-s cmda`), and extracts it directly into your mounted Google Drive project directory.

---

## ⚡ Phase 3: Cloud GPU Simulation Execution (Google Colab)
*Connect to Google Colab to run the rapid smoke test on an NVIDIA GPU.*

### 3.1. Connect to Colab GPU Console
In your local terminal:
```bash
colab console -s cmda
```

Inside the Colab console terminal, run the cloud setup (~90 seconds) and navigate to the project directory:

```bash
# 1. Install cmda, OpenMM (CUDA/CPU), and AmberTools:
curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/dev1/scripts/setup.sh | bash

# 2. Activate environment & navigate to Drive project directory:
source /root/miniforge3/etc/profile.d/conda.sh
conda activate cmda
cd /content/drive/MyDrive/research/projects/colabmda-example
```

---

### 3.2. Fast Equilibration & Production MD (~2 Minutes)
Run quick 100 ps equilibration, 1 ns production MD, and automated trajectory analysis:

```bash
export COLABMDA_RUN_ID=run01

# 1. Run 100 ps equilibration:
bash scripts/02_simulation/02_equilibrate.sh

# 2. Run 1 ns production MD:
bash scripts/02_simulation/03_production.sh

# 3. Merge trajectories & generate RMSD / Rg / RMSF plots:
bash scripts/02_simulation/04_analysis.sh

# 4. Verify output integrity:
python scripts/02_simulation/05_check.py
```

*`05_check.py` validates all 100 trajectory frames, thermodynamic logs, nonfinite coordinate absence, and touches `RUN_DONE` upon passing.*

---

### 3.3. Checkpoint & Timeout Resume Test
To verify that simulation recovery works as expected, request an extension:

```bash
# Test extension request:
COLABMDA_EXAMPLE_NS=2 bash scripts/02_simulation/03_production.sh
```

*ColabMDA detects the existing `prod.chk` checkpoint on Google Drive and resumes dynamics from the last saved timestep.*

---

## 📥 Phase 4: Retrieve Results & Inspection (Local Machine)
*Switch back to your local workstation terminal.*

### 4.1. Pull Analyzed Data from Cloud
Download deliverables to your local machine:

```bash
cd ~/research/projects/colabmda-example
cmda pull -s cmda
```

### 4.2. Inspect Verification Reports
Verify that all deliverables arrived intact:

```bash
cat data/02_simulation/run01/summary.json
```

View the generated analysis plots in `data/02_simulation/run01/apo/WT/04_analysis/figures/`:
* `rmsd_vs_time.png` (Backbone RMSD over 1 ns)
* `rg_vs_time.png` (Radius of gyration compactness)
* `rmsf_per_residue.png` (Per-residue root mean square fluctuations for residues 1–20)

---

## Acknowledgements
- [OpenMM](https://openmm.org) & [PDBFixer](https://github.com/openmm/pdbfixer)
- [Modeller](https://salilab.org/modeller/)
- [MDAnalysis](https://www.mdanalysis.org) & [MDTraj](https://www.mdtraj.org)
- [PDB2PQR & PROPKA](https://server.poissonboltzmann.org/)
- [NumPy](https://numpy.org), [Matplotlib](https://matplotlib.org), [Biopython](https://biopython.org)
- [Google Colab](https://colab.research.google.com) & [Miniforge/Conda](https://github.com/conda-forge/miniforge)

---

## Citation

Initially, this repository was based on a GROMACS-on-Colab workflow, but was transitioned to OpenMM to improve performance and stability for automated simulation recovery on cloud platforms.

If you use this tool, please cite the underlying studies:

> Paul SK, Saddam M, Rahaman KA, Choi JG, Lee SS, Hasan M. **Molecular modeling, molecular dynamics simulation, and essential dynamics analysis of grancalcin: An upregulated biomarker in experimental autoimmune encephalomyelitis mice.** *Heliyon*. 2022 Oct 23;8(10):e11232. doi: [10.1016/j.heliyon.2022.e11232](https://doi.org/10.1016/j.heliyon.2022.e11232). PMID: 36340004; PMCID: PMC9626934.
>
> [![Cited by](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fapi.openalex.org%2Fworks%2Fhttps%3A%2F%2Fdoi.org%2F10.1016%2Fj.heliyon.2022.e11232&query=%24.cited_by_count&label=cited%20by&color=blue)](https://doi.org/10.1016/j.heliyon.2022.e11232)
>
> Paul SK, Saddam M, Tabassum N, Hasan M. **Molecular dynamics simulation of wild and mutant proteasome subunit beta type 8 (PSMB8) protein: Implications for restoration of inflammation in experimental autoimmune encephalomyelitis pathogenesis.** *Heliyon*. 2024 Dec 15;11(1):e41166. doi: [10.1016/j.heliyon.2024.e41166](https://doi.org/10.1016/j.heliyon.2024.e41166). PMID: 39802026; PMCID: PMC11719297.
>
> [![Cited by](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fapi.openalex.org%2Fworks%2Fhttps%3A%2F%2Fdoi.org%2F10.1016%2Fj.heliyon.2024.e41166&query=%24.cited_by_count&label=cited%20by&color=blue)](https://doi.org/10.1016/j.heliyon.2024.e41166)
