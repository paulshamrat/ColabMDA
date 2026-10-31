# Trp-Cage (1L2Y) Apo State: 5-Minute Quickstart & Smoke Test

This tutorial provides a lightweight, 5-minute simulation walkthrough using the 20-residue **Trp-cage miniprotein** (PDB **1L2Y**, NMR Model 1, residues **1–20**). It serves as a rapid "smoke test" to verify that your local environment, Google Drive synchronization, Colab GPU runtime, and simulation checkpointing are fully operational before launching full research production runs.

---

## 📁 Trp-Cage Directory Architecture Map

ColabMDA organizes the Trp-cage project files within a deterministic hierarchy across all four workflow phases:

```text
research/projects/colabmda-example/
├── scripts/                      # Staged fast runner scripts
│   ├── 01_structure/
│   │   └── 01_fetch.py           # Fetch & validate 1L2Y NMR coordinates
│   └── 02_simulation/
│       ├── common.sh             # Environment variables & run ID management
│       ├── 01_prepare.sh         # Solvation & AMBER parameterization (cmda prep)
│       ├── 02_equilibrate.sh     # Quick 100 ps heating & equilibration (cmda sim eq)
│       ├── 03_production.sh      # 1 ns checkpointed production MD (cmda sim md)
│       ├── 04_analysis.sh        # Trajectory merge, autoimaging, and analysis
│       └── 05_check.py           # Strict frame and thermodynamic audit (touches RUN_DONE)
└── data/
    ├── 00_setup/example.json     # Staging metadata & integrity checksums
    ├── 01_structure/trpcage/     # 1l2y.cif, WT.pdb (20 residues), input.json, STRUCTURE_DONE
    └── 02_simulation/<run>/
        ├── provenance.json       # Environment, GPU, and package manifest
        └── apo/WT/
            ├── 01_prep/          # system.prmtop, system.inpcrd, solvated.pdb, PREP_DONE
            ├── 02_eq/            # system.xml, npt.state.xml, equilibration_qc.json, EQ_DONE
            ├── 03_md/r1/         # md.dcd, md.log, prod.chk, MD_DONE, chunks/
            └── 04_analysis/      # rmsd.csv, rg.csv, rmsf.csv, figures/, ANALYSIS_DONE
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

> [!TIP]
> **Modeller Key:** Free for academic research ([register here](https://salilab.org/modeller/registration.html)). Enter it when prompted, or run:  
> `KEY_MODELLER="YOUR_KEY" bash -c "$(curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/dev1/scripts/setup.sh)"`

---

### 1.1. Project Workspace & Example Staging
Create and navigate to your project directory, then stage the Trp-cage runner scripts:

```bash
mkdir -p ~/research/projects/colabmda-example
cd ~/research/projects/colabmda-example

# Stage Trp-cage runner scripts and verified baseline hashes:
cmda stage --example trpcage
```

*What `cmda stage` generates:*
* `scripts/`: Staged execution scripts for structure fetching, preparation, equilibration, production, and analysis.
* `data/00_setup/example.json`: Checksum manifest recording script hashes to prevent accidental overwrites.

---

### 1.2. Structural Fetching & Coordinate Verification
Fetch the 20-residue NMR structure (Model 1, Chain A of 1L2Y):

```bash
# Fetch and validate coordinates using the staged script:
python scripts/01_structure/01_fetch.py

# (Optional) Direct CLI alternative:
# cmda fetch --state apo --pdb-id 1L2Y --chain A --range 1 20
```

*What `01_fetch.py` generates in `data/01_structure/trpcage/`:*

| File | Purpose | Scientific Role |
| :--- | :--- | :--- |
| **`1l2y.cif`** | Archive Copy | Raw crystallographic/NMR record directly from the RCSB PDB archive. |
| **`WT.pdb`** | Standard Model | 20-residue clean coordinates (`NLYIQWLKDGGPSSGRPPPS`) with NMR ensemble tags stripped. |
| **`input.json`** | Provenance Record | Machine-readable log containing parameters, timestamps, and SHA-256 hashes. |
| **`STRUCTURE_DONE`** | Gate Marker | Verification marker confirming all 20 Cα atoms and the sequence match exactly. |

#### Visual Inspection in PyMOL:
```bash
pymol data/01_structure/trpcage/WT.pdb
```

> [!TIP]
> **Trp-Cage Structural Verification:**
> * In PyMOL, enter: `show sticks, resn TRP or resn PRO`.
> * Confirm the central **Trp6** indole ring is enclosed within the hydrophobic "cage" formed by **Pro12**, **Pro17**, **Pro18**, and **Pro19**.
> * Verify that **Asp9** and **Arg16** form a stabilizing salt bridge across the $\alpha$-helix.

---

### 1.3. Structure Preparation & Solvation (`tleap` + PROPKA at pH 7.4)
Solvate the miniprotein with OPC water, neutralize net charge with counterions, and assign AMBER physical parameters:

```bash
export COLABMDA_RUN_ID=run01

# Execute preparation:
bash scripts/02_simulation/01_prepare.sh

# (Optional) Direct CLI equivalent:
# cmda prep --state apo --variant WT --receptor data/01_structure/trpcage/WT.pdb
```

*Key outputs generated in `data/02_simulation/run01/apo/WT/01_prep/`:*
* `model.pdb`: Exact copy of input coordinates for provenance.
* `solvated.pdb`: Solvated coordinate PDB (~1,800 OPC waters + neutralizing Cl⁻ counterion).
* `system.prmtop`: Amber topology containing all `ff19SB` physical force-field parameters.
* `system.inpcrd`: Solvated, charge-neutralized initial coordinates.
* `prep.json`: Parameter and atom-count receipt.
* `PREP_DONE`: Verification completion gate.

#### In-Context Reference: Explicit AMBER Protonation States (PROPKA 3 at pH 7.4)
ColabMDA automatically invokes **PROPKA 3** during `cmda prep` to evaluate 3D electrostatic environments and assign explicit AMBER titratable residue forms:

| Residue Type | Code | Physical Description & Protonation Criteria |
| :--- | :--- | :--- |
| **Histidine** | `HID` | Neutral Histidine protonated exclusively at $\text{N}_\delta$. |
| | `HIE` | Neutral Histidine protonated exclusively at $\text{N}_\epsilon$ (most common neutral form). |
| | `HIP` | Doubly protonated, positively charged Histidine (both nitrogens protonated, $pK_a > 7.4$). |
| **Acidic** | `ASP` / `GLU` | Negatively charged carboxylate forms (standard physiological state). |
| | `ASH` / `GLH` | Neutral protonated carboxylic acid forms (assigned when local microenvironment $pK_a > 7.4$). |
| **Basic** | `LYS` / `ARG` | Positively charged forms (Lys8 and Arg16 in Trp-cage are fully protonated). |
| | `LYN` | Neutral deprotonated Lysine (assigned when buried in hydrophobic core with $pK_a < 7.4$). |
| **Cysteine** | `CYS` | Free neutral thiol form. |
| | `CYX` | Disulfide-bridged Cysteine form. |

---

### 1.4. Visual Quality Assurance in Local PyMOL
Verify the solvated unit cell and counterions locally before uploading:

```bash
pymol data/02_simulation/run01/apo/WT/01_prep/solvated.pdb
```

Inside PyMOL:
1. `show cartoon` $\to$ Confirm the short N-terminal $\alpha$-helix (residues 2–8) and polyproline II tail (residues 12–19).
2. `show spheres, resn CL` $\to$ Confirm the neutralizing chloride counterion is placed in bulk water.
3. `show cell` $\to$ Confirm the truncated octahedron / cubic periodic unit cell padding ($\ge 10\text{ \AA}$ clearance).

---

## 🚀 Phase 2: Launch Colab GPU & Stream Data to Drive
*With models and solvated Amber parameters prepared locally, launch your Colab GPU instance and stream the project directory to Google Drive.*

In your local workstation terminal:

```bash
# 1. Launch a Colab GPU instance:
colab new -s cmda --gpu T4

# 2. Mount Google Drive (visit the printed authorization link once, then press Enter):
colab drivemount --session cmda

# 3. Stream local data/ directly to Google Drive:
cmda push -s cmda
```

> [!TIP]
> **Data Destination on Google Drive:**
> * **Google Drive Filesystem Path:** `/content/drive/MyDrive/research/projects/colabmda-example/data/`
> * **Mechanism:** `cmda push` packages your local `data/` folder, streams it to your active Colab VM (`-s cmda`), and extracts it directly into your mounted Google Drive project directory without manual file uploading.

---

## ⚡ Phase 3: Cloud GPU Simulation Execution (Google Colab)
*Connect to Google Colab to run the rapid smoke test on an NVIDIA GPU.*

### 3.1. Connect to Colab GPU Console
In your local workstation terminal:
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

### 3.2. Fast Equilibration (~30 Seconds)
Execute energy minimization and a fast 100 ps NPT equilibration smoke test:

```bash
export COLABMDA_RUN_ID=run01

# Run fast equilibration:
bash scripts/02_simulation/02_equilibrate.sh

# (Optional) Direct CLI equivalent:
# cmda sim --software-test --state apo --variant WT --stage eq --protocol quick --equil-time 100
```

*What `02_equilibrate.sh` generates in `data/02_simulation/run01/apo/WT/02_eq/`:*
* `system.xml`: Serialized OpenMM system definition.
* `npt.state.xml`: Final relaxed thermodynamic state (velocities, coordinates, box vectors).
* `equilibration_qc.json`: Fast thermodynamic summary (temperature, volume, density).
* `EQ_DONE`: Equilibration completion gate.

---

### 3.2.1. Intermediate Inspection Gate: Minimal Verification Checklist (Local Machine)
*Before committing compute to extended production runs, verify that equilibration completed stably in 4 minimal steps:*

#### 1. Pull Deliverables to Local Machine (Terminal)
Switch to your local workstation terminal:
```bash
cd ~/research/projects/colabmda-example

# Pull lightweight equilibration deliverables (~1 MB):
cmda pull -s cmda
```

#### 2. Verify Status in 1 Second (Terminal)
```bash
cat data/02_simulation/run01/apo/WT/02_eq/equilibration_qc.json
```
*Expected Output:* `T ≈ 300 K`, `Density ≈ 1.00–1.02 g/mL`, `Potential Energy: negative and stable`.

#### 3. Quick 3D Fold Check in PyMOL in 30 Seconds
Launch PyMOL with starting model and relaxed coordinates:
```bash
pymol data/02_simulation/run01/apo/WT/01_prep/model.pdb data/02_simulation/run01/apo/WT/01_prep/solvated.pdb
```
In the PyMOL right panel:
1. Next to **`(all)`**: Click **`H`** (Hide) $\to$ **`waters`**, then **`S`** (Show) $\to$ **`as`** $\to$ **`cartoon`**.
2. Confirm the 20-residue fold remains compact without steric distortion.

---

### 3.3. Production MD (~1 Minute on T4 GPU)
Run 1 ns of checkpointed production molecular dynamics:

```bash
export COLABMDA_RUN_ID=run01

# Run 1 ns production:
bash scripts/02_simulation/03_production.sh

# (Optional) Direct CLI equivalent:
# cmda sim --software-test --state apo --variant WT --stage md --replica r1 --total-ns 1 --checkpoint-ps 100 --traj-interval 10
```

#### In-Context Reference: Simulation Math & Storage Budgeting
Calculate expected trajectory frames:
$$\text{Expected Frames} = \frac{\text{Total Time (ps)}}{\text{Trajectory Saving Interval (ps)}}$$

* **Example:** A $1.0\text{ ns}$ ($1,000\text{ ps}$) simulation with a $10.0\text{ ps}$ saving interval:
  $$\text{Expected Frames} = \frac{1,000\text{ ps}}{10\text{ ps}} = 100\text{ frames}$$
* **Saving Every $10\text{ ps}$ (`--traj-interval 10`):** Provides high temporal resolution for fast conformational fluctuations while maintaining tiny storage footprints.
* **Stripped Trajectories (`cmda merge`):** Removing solvent water leaves only the 20 protein residues (304 atoms), reducing trajectory size to **<1 MB**.

#### 💡 Resuming After a Session Timeout
Production MD continuously saves checkpoints (`prod.chk`) directly to Google Drive every 100 ps. If your session expires:

1. Re-launch the session and re-mount Drive from your local workstation terminal:
   ```bash
   colab new -s cmda --gpu T4
   colab drivemount --session cmda
   colab console -s cmda
   ```
2. In the Colab console, re-activate the environment:
   ```bash
   [ ! -d "/root/miniforge3/envs/cmda" ] && curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/dev1/scripts/setup.sh | bash
   source /root/miniforge3/etc/profile.d/conda.sh
   conda activate cmda
   cd /content/drive/MyDrive/research/projects/colabmda-example
   ```
3. Re-run the production command (or request an extension):
   ```bash
   export COLABMDA_RUN_ID=run01
   COLABMDA_EXAMPLE_NS=2 bash scripts/02_simulation/03_production.sh
   ```
   *ColabMDA automatically detects the existing `prod.chk` checkpoint on Google Drive and seamlessly resumes dynamics.*

#### Trajectory Integrity Auditing (`cmda status`)
Verify chunk continuity, frame counts, and resume readiness:
```bash
cmda status --root . --state apo --variant WT --replica r1
```

<details>
<summary>🔍 <b>Click to view example status report</b></summary>

```text
[STATUS]
  Workdir          : data/02_simulation/run01/apo/WT/03_md/r1
  Chunks (DCD/log) : 1 / 1
  Topology File    : solvated.pdb
                     └─ 5878 atoms, 1860 residues
                     └─ (20 protein, 1839 water, 1 ion)
  Trajectory File  : md.dcd (100 frames)
  Log File         : md.log (100 frames)
  Frames (from logs): 100
  Sim time (ps/ns) : 1000.00 ps / 1.0000 ns
  Resume-ready     : YES (needs prod.chk, system.xml, solvated.pdb)
```

</details>

---

### 3.4. Trajectory Merge, PBC Autoimaging & Analysis
Post-process the trajectory by autoimaging periodic boundaries, stripping bulk water, and generating RMSD, Rg, and RMSF profiles:

```bash
export COLABMDA_RUN_ID=run01

# 1. Merge trajectory & calculate biophysical metrics:
bash scripts/02_simulation/04_analysis.sh

# 2. Run automated validation check:
python scripts/02_simulation/05_check.py
```

*What `04_analysis.sh` and `05_check.py` output:*
* `prod_keep_ligand.dcd`: Solvent-stripped, PBC autoimaged trajectory (304 atoms, 100 frames).
* `rmsd.csv` / `rmsd_vs_time.png`: Backbone Cα RMSD over 1 ns.
* `rg.csv` / `rg_vs_time.png`: Radius of gyration compactness.
* `rmsf.csv` / `rmsf_per_residue.png`: Per-residue fluctuation profile for residues 1–20.
* `RUN_DONE`: Final pipeline validation gate.

---

## 📥 Phase 4: Retrieve Results & PyMOL Visualization (Local Machine)
*Switch back to your local workstation terminal.*

### 4.1. Pull Analyzed Data from Cloud
Download deliverables to your local machine:

```bash
cd ~/research/projects/colabmda-example

# Pull deliverables:
cmda pull -s cmda

# To also pull the solvent-stripped trajectory:
cmda pull -s cmda --traj
```

---

### 4.2. In-Context Reference: Topology Pairing Guidelines
When visualizing or analyzing trajectories, the reference topology must match the atom count in the trajectory file:

| Trajectory File | Description | Matching Topology File | Atom Count | Intended Usage |
| :--- | :--- | :--- | :--- | :--- |
| `md.dcd` | Raw coordinate trajectory | `01_prep/solvated.pdb` or `system.prmtop` | ~5,800 atoms | Raw simulation checkpointing & resume. |
| `prod_keep_ligand.dcd` | Merged, PBC autoimaged, solvent-stripped | `03_md/r1/stripped.prmtop` or `WT.pdb` | 304 atoms | Standard biophysical analysis (`cmda analyze`, PyMOL). |
| `prod_ca.dcd` | Merged Cα backbone coordinates | `03_md/r1/ca.prmtop` | 20 atoms | Rapid PCA and essential dynamics. |

> [!WARNING]
> **Topology Pairing Requirement:**
> Loading `prod_keep_ligand.dcd` (304 atoms) with `solvated.pdb` (~5,800 atoms) in PyMOL or MDAnalysis will trigger an **atom mismatch error**. Always pair stripped trajectories with `stripped.prmtop` or `WT.pdb`.

---

### 4.3. Visualizing Trajectories in Local PyMOL
Animate the 1 ns Trp-cage trajectory in PyMOL:

```bash
# View the production trajectory using the cmda CLI:
cmda view --root . --state apo --variant WT --replica r1
```

Or run an automated PyMOL visualization script:

```bash
pymol data/01_structure/trpcage/WT.pdb data/02_simulation/run01/apo/WT/03_md/r1/prod_keep_ligand.dcd
```

<details>
<summary><b>Click to view custom PyMOL script configuration (<code>visualize_trpcage.pml</code>)</b></summary>

```text
# 1. Clear out memory & set white background
reinitialize
bg_color white

# 2. Load topology and trajectory
load data/01_structure/trpcage/WT.pdb, trpcage
load_traj data/02_simulation/run01/apo/WT/03_md/r1/prod_keep_ligand.dcd, trpcage

# 3. Fit backbone frame-by-frame
intra_fit trpcage and name CA

# 4. Display secondary structure cartoon
show cartoon, trpcage
color cyan, trpcage and name C*
util.cnc("trpcage")

# 5. Highlight central Trp6 and Polyproline cage
show sticks, trpcage and (resi 6 or resi 12 or resi 17-19) and not name N+C+O+H
color yellow, trpcage and resi 6 and name C*
color orange, trpcage and (resi 12 or resi 17-19) and name C*
set stick_radius, 0.22

# 6. Highlight Asp9 - Arg16 Salt Bridge
show sticks, trpcage and (resi 9 or resi 16) and not name N+C+O+H
distance salt_bridge, (trpcage and resi 9 and name OD*), (trpcage and resi 16 and name NH*), 4.0

# 7. Center camera & play movie
zoom trpcage, buffer=3
mplay
```

</details>

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

---

<p align="center">
  <img src="https://komarev.com/ghpvc/?username=paulshamrat&repo=ColabMDA&label=Visitors&color=blue&style=flat-square" alt="Visitors" />
</p>
