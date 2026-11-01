# Trp-Cage (1L2Y) Apo State: 5-Minute Quickstart & Smoke Test

This tutorial provides a lightweight, 5-minute simulation walkthrough using the 20-residue **Trp-cage miniprotein** (PDB **1L2Y**, NMR Model 1, residues **1–20**). It serves as a rapid "smoke test" to verify that your local environment, Google Drive synchronization, Colab GPU runtime, and simulation checkpointing are fully operational before launching full research production runs.

---

## 📁 Trp-Cage Directory Architecture Map

ColabMDA manages project files within a standardized hierarchy under `data/apo/` across all four workflow phases:

```text
research/projects/trpcage-example/
└── data/
    └── apo/
        ├── structure/
        │   └── 1L2Y/                     # Target: NMR PDB 1L2Y
        │       ├── raw/                  # Phase 1.2: 1l2y.cif, template.pdb, uniprot.fasta, residues.json
        │       ├── variants/             # Phase 1.3: Modeled coordinates (WT/WT.pdb)
        │       └── STRUCTURE_DONE        # Verification QC gate marker
        │
        ├── simulation/
        │   └── WT/
        │       ├── 01_prep/                  # Solvated Amber system
        │       │   ├── model.pdb             # Input coordinates
        │       │   ├── solvated.pdb          # Solvated visualization PDB
        │       │   ├── system.prmtop         # Amber parameter/topology file
        │       │   ├── system.inpcrd         # Initial orthogonal coordinates
        │       │   └── system.xml            # Serialized OpenMM system
        │       │
        │       ├── 02_eq/                    # Equilibration hub
        │       │   ├── eq.pdb                # Final relaxed coordinates
        │       │   ├── system.xml            # OpenMM system topology & parameters
        │       │   ├── state.xml             # OpenMM thermodynamic state snapshot
        │       │   └── stages/               # Granular restart checkpoints (state_01.xml ... state_24.xml)
        │       │
        │       └── 03_md/                    # Production replica (1 ns smoke test)
        │           └── r1/                   # Replica 1 hub
        │               ├── md.dcd            # Continuous merged trajectory (full system)
        │               ├── md.log            # Continuous merged thermodynamic log
        │               ├── state.xml         # Latest OpenMM simulation state snapshot
        │               ├── prod.chk          # OpenMM checkpoint binary for resume
        │               ├── stripped.prmtop   # Solvent-stripped topology for PyMOL & analysis
        │               ├── prod_keep_ligand.dcd # Continuous PBC-autoimaged stripped trajectory
        │               ├── ca.prmtop         # Cα backbone topology
        │               ├── prod_ca.dcd       # Cα backbone trajectory
        │               ├── MERGE_DONE        # Merge completion marker
        │               └── chunks/           # Raw crash-recovery chunks (0001.dcd, 0001.log, ...)
        │
        └── analysis/                         # Mirrored deliverables (Zero JSON policy)
            └── WT/
                ├── 02_eq/                    # qc.png, qc.txt (thermodynamic & structural convergence)
                └── 03_md/
                    └── r1/                   # rmsd.png, rg.png, rmsf.png, summary.txt, *.csv
```

---

## 💻 Phase 1: Local Modeling & Structure Preparation (Local Machine)
*Execute this phase on your local laptop or workstation CPU. No GPU compute units are required.*

### 1.0. Local Machine Setup (One-Time)
Run this one-time setup in your terminal on **Linux (Ubuntu)** or Windows via **[WSL2](https://learn.microsoft.com/en-us/windows/wsl/install)**:

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

### 1.1. Project Workspace
Create and navigate to your project directory:

```bash
mkdir -p ~/research/projects/trpcage-example
cd ~/research/projects/trpcage-example
```

---

### 1.2. Structural Fetching & Coordinate Verification
Fetch the 20-residue NMR structure (Model 1, Chain A of 1L2Y):

```bash
# Fetch and validate coordinates directly using cmda fetch:
cmda fetch --state apo --pdb-id 1L2Y --chain A

# (Optional) Visually inspect the cleaned template in local PyMOL:
pymol data/apo/structure/1L2Y/raw/template.pdb
```

*What `cmda fetch` generates in `data/apo/structure/1L2Y/`:*

| File | Purpose | Scientific Role |
| :--- | :--- | :--- |
| **`1l2y.cif`** | Archive Copy | Raw crystallographic/NMR record directly from the RCSB PDB archive. |
| **`template.pdb`** | Standard Model | 20-residue clean coordinates (`NLYIQWLKDGGPSSGRPPPS`) with NMR ensemble tags stripped. |
| **`uniprot.fasta`** | Sequence Record | Clean FASTA sequence used for coordinate mapping. |
| **`residues.json`** | Coordinate Metadata | Residue-by-residue mapping verifying 100% full coordinate density (20/20 resolved). |
| **`input.json`** | Provenance Record | Machine-readable log containing parameters, timestamps, and SHA-256 hashes. |
| **`STRUCTURE_DONE`** | Gate Marker | Verification marker confirming all 20 Cα atoms and the sequence match exactly. |

---

### 1.3. Homology Model & Coordinate Verification (`cmda model`)
Verify the canonical Wild-Type coordinates:

```bash
# Verify the canonical Wild-Type model:
cmda model --state apo --variant WT

# Inspect the modeled WT structure in PyMOL:
pymol data/apo/structure/1L2Y/variants/WT/WT.pdb
```

> [!TIP]
> **Trp-Cage Structural Verification:**
> * In PyMOL, enter: `show sticks, resn TRP or resn PRO`.
> * Confirm the central **Trp6** indole ring is enclosed within the hydrophobic "cage" formed by **Pro12**, **Pro17**, **Pro18**, and **Pro19**.
> * Verify that **Asp9** and **Arg16** form a stabilizing salt bridge across the $\alpha$-helix.

---

### 1.4. Structure Preparation & Solvation (`tleap` + PROPKA at pH 7.4)
Solvate the miniprotein with OPC water, neutralize net charge with counterions, and assign AMBER physical parameters:

```bash
# Solvate with OPC water, neutralize net charge, and assign AMBER parameters
cmda prep --state apo --variant WT

# Visually inspect the solvated system in PyMOL:
pymol data/apo/simulation/WT/01_prep/solvated.pdb
```

*Key outputs generated in `data/apo/simulation/WT/01_prep/`:*
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

### 1.5. Visual Quality Assurance in Local PyMOL
Verify the solvated unit cell and counterions locally before uploading:

```bash
# 1. Inspect modeled coordinates:
cmda view --state apo --variant WT

# 2. Inspect the fully solvated system (water box + neutralizing counterion):
pymol data/apo/simulation/WT/01_prep/solvated.pdb
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
> * **Google Drive Filesystem Path:** `/content/drive/MyDrive/research/projects/trpcage-example/data/`
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
cd /content/drive/MyDrive/research/projects/trpcage-example
```

---

### 3.2. Stage Simulation Hubs & Equilibration (~30 Seconds)
Stage simulation directories and execute the equilibration smoke test:

```bash
# 1. Initialize simulation directories:
cmda stage --state apo --variant WT

# 2. Run quick equilibration (automatically executes EM -> NVT -> NPT -> QC):
cmda sim eq --state apo --variant WT --protocol quick
```

*What `cmda sim eq` generates in `data/apo/simulation/WT/02_eq/`:*
* `eq.pdb`: Relaxed coordinate snapshot.
* `system.xml`: Serialized OpenMM system definition.
* `state.xml`: Final relaxed thermodynamic state (velocities, coordinates, box vectors).
* `EQ_DONE`: Equilibration completion gate.

---

### 3.2.1. Intermediate Inspection Gate: Minimal Verification Checklist (Local Machine)
*Before committing compute to extended production runs, verify that equilibration completed stably in 4 minimal steps:*

#### 1. Pull Deliverables to Local Machine (Terminal)
Switch to your local workstation terminal:
```bash
cd ~/research/projects/trpcage-example

# Pull lightweight equilibration deliverables (~1 MB):
cmda pull -s cmda --stage eq
```

#### 2. Verify Status in 1 Second (Terminal)
```bash
cat data/apo/analysis/WT/02_eq/qc.txt
```
*Expected Output:* `Status: PASSED`, `Cα RMSD: ~0.8 Å` (Confirms temperature $\approx 300\text{ K}$, density $\approx 1.01\text{ g/mL}$, fold intact).

#### 3. Glance at 6-Panel QC Plot in 5 Seconds
Open the unified publication plot in any image viewer:
* `data/apo/analysis/WT/02_eq/qc.png`

*Expected Visual:*
* **Panels A–D (Thermodynamics):** Temperature is flat at **$300\text{ K}$**; density is flat at **$\sim 1.01\text{ g/mL}$**; potential energy and volume are stable.
* **Panel E (Per-Residue C$\alpha$ Drift):** Core residues stay flat at **$< 0.5\text{ \AA}$**; termini remain below threshold ($< 1.0\text{ \AA}$).
* **Panel F (3D Backbone Alignment):** Solid cyan equilibrated fold matches semi-transparent salmon starting model ($\text{C}_\alpha\text{ RMSD} < 1.0\text{ \AA}$).

#### 4. Quick 3D Fold Check in PyMOL in 30 Seconds
Launch PyMOL with both starting and equilibrated coordinates:
```bash
pymol data/apo/simulation/WT/01_prep/model.pdb data/apo/simulation/WT/02_eq/eq.pdb
```
In the PyMOL right panel:
1. Next to **`(all)`**: Click **`H`** (Hide) $\to$ **`waters`**, then **`S`** (Show) $\to$ **`as`** $\to$ **`cartoon`**.
2. Next to **`eq`**: Click **`A`** (Action) $\to$ **`align`** $\to$ **`to molecule`** $\to$ choose **`model`**.

*Expected Visual:* Cyan equilibrated backbone (`eq`) overlays closely onto starting model (`model`). No unfolded helices or broken bonds.

---

### 3.3. Production MD (~1 Minute on T4 GPU)
Run 1 ns of checkpointed production molecular dynamics:

```bash
# Run 1 ns production simulation:
cmda sim md --state apo --variant WT --replica r1 --total-ns 1 --checkpoint-ps 100 --traj-interval 10
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
   cd /content/drive/MyDrive/research/projects/trpcage-example
   ```
3. Re-run the exact same `cmda sim md` command:
   ```bash
   cmda sim md --state apo --variant WT --replica r1 --total-ns 2
   ```
   *ColabMDA automatically detects the existing `prod.chk` checkpoint on Google Drive and seamlessly resumes dynamics.*

#### Trajectory Integrity Auditing (`cmda status`)
Verify chunk continuity, frame counts, and resume readiness:
```bash
cmda status --state apo --variant WT --replica r1
```

<details>
<summary>🔍 <b>Click to view example status report</b></summary>

```text
[STATUS]
  Workdir          : data/apo/simulation/WT/03_md/r1
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

### 3.4. Merge Trajectories & Autoimage PBC
Merge raw chunked trajectories into a single continuous file, re-center the peptide in the periodic unit cell (`autoimage`), and strip bulk solvent:

```bash
# Merge Wild-Type trajectory:
cmda merge --state apo --variant WT
```

*What `cmda merge` outputs in `data/apo/simulation/WT/03_md/r1/`:*
* `md.dcd`: Clean continuous trajectory merged from chunk files.
* `md.log`: Continuous merged thermodynamic log table.
* `state.xml`: Final OpenMM physical restart state.

---

### 3.5. Trajectory Analysis
Calculate backbone RMSD, Radius of Gyration, and per-residue RMSF profiles:

```bash
# Run biophysical analysis:
cmda analyze --state apo --variant WT
```

Outputs saved directly in `data/apo/analysis/WT/03_md/r1/`:
* `summary.txt` (Summary table with Mean ± SD for RMSD, Rg, Temp)
* `rmsd.png` (Backbone Cα RMSD over 1 ns)
* `rg.png` (Radius of Gyration compactness)
* `rmsf.png` (Per-residue root mean square fluctuations for residues 1–20)

---

## 📥 Phase 4: Retrieve Results & PyMOL Visualization (Local Machine)
*Switch back to your local workstation terminal.*

### 4.1. Pull Analyzed Data from Cloud
Download analyzed figures, summary tables, and reports:

```bash
cd ~/research/projects/trpcage-example

# Pull deliverables:
cmda pull -s cmda

# To also pull post-processed stripped trajectories for PyMOL inspection:
cmda pull -s cmda --traj
```

> [!TIP]
> **What `cmda pull` Does:**
> 1. Packages analysis tables, logs, summary manifests, and PNG plots from Google Drive.
> 2. Excludes heavy raw trajectory files (`chunks/*.dcd`, `.chk`) by default so the download completes quickly. Pass `--traj` to include stripped trajectories (`prod_keep_ligand.dcd`), or `--all` for full raw simulation chunks.
> 3. Extracts files non-destructively into your local `data/` folder.

---

### 4.2. In-Context Reference: Topology Pairing Guidelines
When visualizing or analyzing trajectories, the reference topology file used must match the coordinates in the trajectory file:

| Trajectory File | Description | Matching Topology File | Atom Count | Intended Usage |
| :--- | :--- | :--- | :--- | :--- |
| `md.dcd` | Raw chunked coordinates | `01_prep/solvated.pdb` or `system.prmtop` | ~5,800 atoms | Raw simulation checkpointing & resume. |
| `prod_keep_ligand.dcd` | Merged, PBC autoimaged, solvent-stripped | `03_md/r1/stripped.prmtop` or `WT.pdb` | 304 atoms | Standard biophysical analysis (`cmda analyze`, PyMOL). |
| `prod_ca.dcd` | Merged Cα backbone only | `03_md/r1/ca.prmtop` | 20 atoms | Rapid PCA and essential dynamics. |

> [!WARNING]
> **Topology Pairing Requirement:**
> Loading `prod_keep_ligand.dcd` (304 atoms) with `solvated.pdb` (~5,800 atoms) in PyMOL or MDAnalysis will trigger an **atom mismatch error**. Always pair stripped trajectories with `stripped.prmtop` or `WT.pdb`.

---

### 4.3. Visualizing Trajectories in Local PyMOL
View the production trajectory using the `cmda` CLI:

```bash
# 1. View production trajectory in local PyMOL:
cmda view --state apo --variant WT --replica r1

# 2. Or launch PyMOL manually:
pymol data/apo/structure/1L2Y/variants/WT/WT.pdb data/apo/simulation/WT/03_md/r1/prod_keep_ligand.dcd
```

<details>
<summary><b>Click to view custom PyMOL script configuration (<code>visualize_trpcage.pml</code>)</b></summary>

```text
# 1. Clear out memory & set white background
reinitialize
bg_color white

# 2. Load topology and trajectory
load data/apo/structure/1L2Y/variants/WT/WT.pdb, trpcage
load_traj data/apo/simulation/WT/03_md/r1/prod_keep_ligand.dcd, trpcage

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
