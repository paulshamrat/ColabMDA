# CDKL5 Apo State: Complete Step-by-Step Walkthrough

This tutorial provides the complete, end-to-end simulation workflow for the **Apo state** (unliganded protein in water) using human cyclin-dependent kinase-like 5 (**CDKL5**, UniProt **O76039**, PDB **4BGQ** Chain A, residues **1–303**) with Wild-Type (WT) and the disease-associated mutant **L119R** as the comparative pair.

---

Outputs are shown after each stage. The complete directory tree is collapsed at
the end of this tutorial. Paths are relative to the project root unless stated otherwise.

---

## 💻 Phase 1: Local Modeling & Structure Preparation (Local Machine)
*Execute this phase on your local laptop or workstation CPU. No GPU is required.*

### 1.0. Local Machine Setup (One-Time)
Run this one-time setup in your terminal on **Linux (Ubuntu)** or Windows via **[WSL2](https://learn.microsoft.com/en-us/windows/wsl/install)**:

```bash
# 1. Install required system tools (Ubuntu/WSL2):
sudo apt-get update
sudo apt-get install -y git curl wget ca-certificates

# 2. Install cmda, OpenMM, AmberTools, PyMOL, and Modeller:
bash -c "$(curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/main/scripts/setup.sh)"

# 3. Activate the environment in this terminal:
source ~/miniforge3/etc/profile.d/conda.sh
conda activate cmda

# 4. Log in to Google Colab (one-time browser login):
colab sessions
```

Run the activation commands in every new terminal or SSH session. The installer
runs in a child shell, so it cannot activate the environment in your current terminal.

When `colab sessions` prints an authorization URL:

1. Keep that terminal open and open the printed URL in a browser. Over SSH, use
   a browser on your own computer.
2. Sign in with the Google account you use for Colab and approve access.
3. Copy the authorization code Google displays.
4. Paste the code into the same terminal at `Enter the authorization code:` and
   press Enter. Keep the code private.

> [!TIP]
> **Modeller Key:** Free for academic research ([register here](https://salilab.org/modeller/registration.html)). Enter it when prompted, or run:
> `KEY_MODELLER="YOUR_KEY" bash -c "$(curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/main/scripts/setup.sh)"`

### 1.1. Project Workspace
Create and navigate to your project directory:

```bash
mkdir -p ~/research/projects/cmda-example
cd ~/research/projects/cmda-example
```

### 1.2. Structural Fetching & Biological UniProt Mapping
The `cmda fetch` command downloads experimental PDB coordinates, maps canonical UniProt residue numbering, and detects unresolved loop coordinates for subsequent modeling:

```bash
# Fetch CDKL5 kinase domain (residues 1–303)
cmda fetch --state apo --pdb-id 4BGQ --uniprot-id O76039 --range 1 303

# (Optional) Visually inspect the cleaned template in local PyMOL:
pymol data/apo/structure/O76039_4BGQ/raw/template.pdb
```

<details>
<summary><b>View stage directory tree</b></summary>

```text
data/apo/structure/O76039_4BGQ/
└── raw/
    ├── 4bgq.cif                                     # Experimental archive; PDB fallback is possible
    ├── template.pdb                                 # Cleaned template coordinates
    └── uniprot.fasta                                # Sequence used for modeling
```

Stage log: `logs/apo/structure/O76039_4BGQ/structure.log`.

</details>

<details>
<summary>🔍 <b>Click to view example structural verification banner</b></summary>

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
 Outputs Saved:    data/apo/structure/O76039_4BGQ/
==============================================================================
```

</details>

---

### 1.3. Homology Modeling & Point Mutations (Modeller)
Model the complete Wild-Type coordinates, closing all 27 missing loop residues, then introduce point mutations on the verified WT backbone:

```bash
# 1. Build full Wild-Type model (UniProt canonical, fills all 27 missing loop residues)
cmda model --state apo --variant WT

# 2. Visually inspect the modeled WT structure in PyMOL:
pymol data/apo/structure/O76039_4BGQ/variants/WT/WT.pdb

# 3. Build Point Mutant (L119R) on the verified WT backbone
cmda model --state apo --variant L119R

# 4. Compare WT and Mutant (L119R) side-by-side in PyMOL:
pymol data/apo/structure/O76039_4BGQ/variants/WT/WT.pdb data/apo/structure/O76039_4BGQ/variants/L119R/L119R.pdb
```

<details>
<summary><b>View stage directory tree</b></summary>

```text
data/apo/structure/O76039_4BGQ/variants/
├── WT/
│   └── WT.pdb                                       # Coordinates passed to preparation
└── L119R/
    └── L119R.pdb                                    # Coordinates passed to preparation
```

Stage logs: `logs/apo/structure/O76039_4BGQ/variants/WT/model.log`, `logs/apo/structure/O76039_4BGQ/variants/L119R/model.log`.

</details>

> [!TIP]
> **PyMOL Model Verification Checks**:
> * **Loop Closure Check**: In `WT.pdb`, confirm that residues 163–173 (activation loop) and 46–52 are modeled continuously as an unbroken cartoon ribbon.
> * **Mutation Rotamer Check**: In the WT vs. L119R comparison, enter `zoom resi 119` and `show sticks, resi 119` in PyMOL. Verify that Arg119 is oriented outward into solvent without steric clashes with neighboring helices.

*What happens behind the scenes:*
* **Engineered Template Mutations**: If crystallographers introduced mutations into the PDB template, `cmda fetch` detects them, and Modeller automatically reverts them to canonical UniProt WT during `cmda model WT`.
* **`cmda model WT`**: Invokes Modeller to build full coordinates for residues 1–303, closing missing loops and preserving author residue numbering.
* **`cmda model L119R`**: Introduces the point mutation onto the verified WT template and optimizes local sidechain rotamers.

---

### 1.4. Structure Preparation & Solvation (`tleap` + PROPKA at pH 7.4)
Prepare the solvated systems locally using CPU-bound AmberTools. By default, ColabMDA neutralizes system charge with counterions only:

```bash
# Solvate with OPC water, neutralize net charge, assign AMBER titratable protonation states
cmda prep --state apo --variant WT
cmda prep --state apo --variant L119R

# Visually inspect the solvated systems in PyMOL (WT and Mutant):
pymol data/apo/simulation/WT/prep/solvated.pdb
pymol data/apo/simulation/L119R/prep/solvated.pdb

# Optional: Add 150 mM physiological saline if specifically desired:
# cmda prep --state apo --variant WT --salt 0.15
```

<details>
<summary><b>View stage directory tree</b></summary>

```text
data/apo/simulation/
├── WT/
│   └── prep/
│       ├── model.pdb                                # Prepared solute coordinates
│       ├── solvated.pdb                             # Protein, water and ions
│       ├── system.prmtop                            # Amber topology passed to equilibration
│       └── system.inpcrd                            # Amber coordinates passed to equilibration
└── L119R/
    └── prep/
        ├── model.pdb                                # Prepared solute coordinates
        ├── solvated.pdb                             # Protein, water and ions
        ├── system.prmtop                            # Amber topology passed to equilibration
        └── system.inpcrd                            # Amber coordinates passed to equilibration
```

Stage logs: `logs/apo/simulation/WT/prep/prep.log`, `logs/apo/simulation/L119R/prep/prep.log`.

</details>

> [!TIP]
> **PyMOL Solvation & Neutrality Check**:
> Open `solvated.pdb` in PyMOL to verify box dimensions and ion distribution:
> * **Water Box Buffer**: Verify the water box provides at least ~14 Å clearance around the protein edges.
> * **Counterions**: Run `select ions, resn CL- or resn NA+` followed by `show spheres, ions` to ensure neutralizing counterions are solvated in bulk water without steric overlaps.

#### In-Context Reference: Explicit AMBER Protonation States (PROPKA 3 at pH 7.4)
Accurate protonation states are critical for biophysical fidelity. ColabMDA natively integrates **PDBFixer**, **PDB2PQR**, and **PROPKA 3** during `cmda prep` to evaluate 3D electrostatic microenvironments and assign explicit AMBER titratable residue forms:

| Residue Type | Code | Physical Description & Protonation Criteria |
| :--- | :--- | :--- |
| **Histidine** | `HID` | Neutral Histidine protonated exclusively at $\text{N}_\delta$. |
| | `HIE` | Neutral Histidine protonated exclusively at $\text{N}_\epsilon$ (most common neutral form). |
| | `HIP` | Doubly protonated, positively charged Histidine (both nitrogens protonated, $pK_a > 7.4$). |
| **Acidic** | `ASP` / `GLU` | Negatively charged carboxylate forms (standard physiological state). |
| | `ASH` / `GLH` | Neutral protonated carboxylic acid forms (assigned when local microenvironment $pK_a > 7.4$). |
| **Basic** | `LYS` / `ARG` | Positively charged forms (standard physiological state). |
| | `LYN` | Neutral deprotonated Lysine (assigned when buried in hydrophobic core with $pK_a < 7.4$). |
| **Cysteine** | `CYS` | Free neutral thiol form. |
| | `CYX` | Disulfide-bridged Cysteine form. |

*Audit Manifest:* Every assignment is permanently recorded.

---

### 1.5. Visual Quality Assurance in Local PyMOL
Verify coordinates and the solvated water box interactively before uploading:

```bash
# 1. Inspect modeled coordinates (WT and Mutant):
cmda view --state apo --variant WT
cmda view --state apo --variant L119R

# 2. Compare WT and Mutant side-by-side in PyMOL:
pymol data/apo/structure/O76039_4BGQ/variants/WT/WT.pdb data/apo/structure/O76039_4BGQ/variants/L119R/L119R.pdb

# 3. Inspect the fully solvated systems (water box + neutralizing counterions):
pymol data/apo/simulation/WT/prep/solvated.pdb
pymol data/apo/simulation/L119R/prep/solvated.pdb
```
*Visual Check:* Confirm loop 163–173 is closed smoothly, Arg119 points into solvent with zero steric clashes, and the OPC water box has adequate padding (~14 Å).

---

## 🚀 Phase 2: Launch Colab GPU & Stream Data to Drive
*Now that your models and solvated Amber parameters are built and verified locally, launch your Colab GPU instance and stream the coordinates to Google Drive.*

In your local terminal:

```bash
# 1. Launch a persistent Colab GPU instance
colab new -s cmda --gpu T4

# 2. Mount Google Drive (visit the printed OAuth URL once, then press Enter)
colab drivemount --session cmda

# 3. Stream local data/ directly to Google Drive in ONE single command:
cmda push -s cmda
```

<details>
<summary><b>View transferred directory tree</b></summary>

```text
/content/drive/MyDrive/research/projects/cmda-example/
└── data/apo/
    ├── structure/                                   # Transferred structural inputs and models
    └── simulation/                                  # Transferred prepared systems
```

</details>

> [!TIP]
> **Data Destination on Google Drive:**
> * **Google Drive Filesystem Path:** `/content/drive/MyDrive/research/projects/cmda-example/data/`
> * **Mechanism:** `cmda push` packages your local `data/` folder, streams it to your active Colab VM (`-s cmda`), and extracts it directly into your mounted Google Drive project directory.

---

## ⚡ Phase 3: Cloud GPU Simulation & Analysis (Google Colab)
*Connect to Google Colab to run OpenMM molecular dynamics simulations on an NVIDIA GPU.*

### 3.1. Connect to Colab GPU Terminal
In your local terminal:
```bash
colab ssh -s cmda
```

Inside the Colab terminal, run the cloud setup (the first installation can take several minutes) and navigate to the project directory:

```bash
# 1. Install cmda, OpenMM (CUDA/CPU), and AmberTools:
curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/main/scripts/setup.sh | bash

# 2. Activate environment & navigate to Drive project directory:
source /root/miniforge3/etc/profile.d/conda.sh
conda activate cmda
cd /content/drive/MyDrive/research/projects/cmda-example
```

---

### 3.2. Stage & Run 24-Stage Equilibration
Stage the simulation directories and execute the 24-stage progressive heating and restraint release ladder (EM, NVT heating, NPT density relaxation) on CUDA GPU:

```bash
# 1. Initialize simulation directories
cmda stage --state apo --variant WT
cmda stage --state apo --variant L119R

# 2. Run full 24-stage equilibration for WT (automatically executes EM -> NVT -> NPT -> QC)
cmda sim eq --state apo --variant WT

# 3. Run full 24-stage equilibration for L119R (automatically executes EM -> NVT -> NPT -> QC)
cmda sim eq --state apo --variant L119R

# (Optional) Re-run standalone QC gating without re-simulating:
# cmda sim qc --state apo --variant WT
# cmda sim qc --state apo --variant L119R
```

<details>
<summary><b>View stage directory tree</b></summary>

```text
data/apo/
├── simulation/
│   ├── WT/
│   │   └── eq/
│   │       ├── solvated.pdb                         # Topology coordinates used by OpenMM
│   │       ├── system.xml                           # Serialized OpenMM system
│   │       ├── em.state.xml                         # Minimized state
│   │       ├── nvt.state.xml                        # State after NVT equilibration
│   │       ├── npt.state.xml                        # Equilibrated state passed to production
│   │       ├── eq.pdb                               # Final equilibrated coordinates
│   │       └── stages/
│   │           ├── state_01.xml                     # First staged restart
│   │           └── state_24.xml                     # Last of 24 staged restart states
│   └── L119R/
│       └── eq/
│           ├── solvated.pdb                         # Topology coordinates used by OpenMM
│           ├── system.xml                           # Serialized OpenMM system
│           ├── em.state.xml                         # Minimized state
│           ├── nvt.state.xml                        # State after NVT equilibration
│           ├── npt.state.xml                        # Equilibrated state passed to production
│           ├── eq.pdb                               # Final equilibrated coordinates
│           └── stages/
│               ├── state_01.xml                     # First staged restart
│               └── state_24.xml                     # Last of 24 staged restart states
└── analysis/
    ├── WT/
    │   └── eq/
    │       ├── qc.txt                               # Equilibration diagnostics
    │       └── qc.png                               # Thermodynamic and structural QC figure
    └── L119R/
        └── eq/
            ├── qc.txt                               # Equilibration diagnostics
            └── qc.png                               # Thermodynamic and structural QC figure
```

</details>

---

### 3.2.1. Intermediate Inspection Gate: Minimal Verification Checklist (Local Machine)
*Before committing cloud compute to 100 ns production runs, verify that equilibration completed stably in 4 minimal steps:*

#### 1. Pull Deliverables to Local Machine (Terminal)
```bash
cd ~/research/projects/cmda-example

# Pull lightweight equilibration deliverables (~2 MB):
cmda pull -s cmda --stage eq

# Or pull only a specific variant:
# cmda pull -s cmda --stage eq --variant WT
```

<details>
<summary><b>View transferred directory tree</b></summary>

```text
data/apo/
├── simulation/                                      # Retrieved simulation deliverables
└── analysis/                                        # Retrieved equilibration QC
```

</details>

> [!TIP]
> **GPU vs. CPU Session for Retrieval:**
> * **Active Session:** If your GPU session is still running, pull directly.
> * **Session Ended / Timed Out:** You **do not** need a GPU to pull results. Simply launch an instant, free CPU runtime (`colab new -s cmda`), mount Drive (`colab drivemount --session cmda`), and run `cmda pull`.

#### 2. Verify Status in 1 Second (Terminal)
```bash
cat data/apo/analysis/WT/eq/qc.txt
cat data/apo/analysis/L119R/eq/qc.txt
```
*Expected Output:* `Status: PASSED`, `Cα RMSD: ~0.83 Å` (Confirms temperature $\approx 300\text{ K}$, density $\approx 1.01\text{ g/mL}$, fold intact).

#### 3. Glance at 6-Panel QC Plot in 5 Seconds
Open the unified 7-inch publication plot in any image viewer:
* `data/apo/analysis/WT/eq/qc.png`
* `data/apo/analysis/L119R/eq/qc.png`

*Expected Visual:*
* **Panels A–D (Thermodynamics):** Temperature is flat at **$300\text{ K}$**; density is flat at **$\sim 1.01\text{ g/mL}$**; potential energy and volume are stable.
* **Panel E (Per-Residue C$\alpha$ Drift):** Core secondary structures stay flat at **$< 0.5\text{ \AA}$**; loop peaks remain below threshold ($< 1.2\text{ \AA}$).
* **Panel F (3D Backbone Alignment):** Solid cyan equilibrated fold matches semi-transparent salmon starting model ($\text{C}_\alpha\text{ RMSD} < 1.5\text{ \AA}$).

#### 4. Quick 3D Fold Check in PyMOL in 30 Seconds
Launch PyMOL with both starting and equilibrated coordinates:
```bash
pymol data/apo/simulation/WT/prep/model.pdb data/apo/simulation/WT/eq/eq.pdb
```
In the PyMOL right panel:
1. Next to **`(all)`**: Click **`H`** (Hide) $\to$ **`waters`**, then **`S`** (Show) $\to$ **`as`** $\to$ **`cartoon`**
2. Next to **`eq`**: Click **`A`** (Action) $\to$ **`align`** $\to$ **`to molecule`** $\to$ choose **`model`**

*Expected Visual:* Cyan equilibrated backbone (`eq`) overlays closely onto starting model (`model`). No unfolded helices or broken loops ($\text{C}_\alpha\text{ RMSD} < 1.5\text{ \AA}$). Repeat for L119R.

<details>
<summary><b>🔍 Advanced PyMOL GUI Inspection (Periodic Box & Ions)</b></summary>

#### Fast PyMOL Inspection Commands
Type these commands directly into the PyMOL command prompt (`PyMOL>`):

```pml
# 1. Clean view & align backbones (salmon = starting model, cyan = equilibrated)
hide (resn HOH,WAT)
show cartoon
color salmon, model_
color cyan, eq
align eq, model_

# 2. Display the periodic unit cell boundary box wireframe
show cell

# 3. Display counter-ions and cofactors/ligands (if present)
show spheres, (resn NA,CL,MG)
show sticks, (organic or resn ATP)
```

> [!NOTE]
> **Why are waters seen outside the cyan box in `eq.pdb`?**
> OpenMM and AMBER integrate in continuous Cartesian space. Molecules diffusing across boundaries remain unwrapped in raw coordinate snapshots. In **Phase 4.2**, **`cmda merge`** automatically runs PBC autoimaging (folding all molecules back inside the unit cell) and solvent stripping for clean trajectory animation.
> *(Note: `eq.pdb` is a static 3D coordinate snapshot; multi-frame animation is recorded in Phase 3.3 and viewed via `cmda view`).*

</details>

*Once verified, proceed to launch Production MD (Section 3.3).*

---

### 3.3. Production MD Across Triplicates (`md/{r1, r2, r3}`)
Replicas fork from `eq/` with distinct random velocity seeds into `md/`:

```bash
# Wild-Type Triplicates (100 ns production runs)
cmda sim md --state apo --variant WT --replica r1 --total-ns 100
cmda sim md --state apo --variant WT --replica r2 --total-ns 100
cmda sim md --state apo --variant WT --replica r3 --total-ns 100

# Mutant Triplicates (100 ns production runs)
cmda sim md --state apo --variant L119R --replica r1 --total-ns 100
cmda sim md --state apo --variant L119R --replica r2 --total-ns 100
cmda sim md --state apo --variant L119R --replica r3 --total-ns 100
```

<details>
<summary><b>View stage directory tree</b></summary>

```text
data/apo/simulation/
├── WT/
│   └── md/
│       └── r1/
│           ├── system.xml                           # Copied OpenMM system
│           ├── solvated.pdb                         # Copied topology coordinates
│           ├── npt.state.xml                        # Starting state copied from equilibration
│           ├── prod.chk                             # Checkpoint used to resume production
│           └── chunks/
│               ├── 0001.dcd                         # First committed trajectory chunk
│               └── 0001.log                         # Matching committed thermodynamic log
└── L119R/
    └── md/
        └── r1/
            ├── system.xml                           # Copied OpenMM system
            ├── solvated.pdb                         # Copied topology coordinates
            ├── npt.state.xml                        # Starting state copied from equilibration
            ├── prod.chk                             # Checkpoint used to resume production
            └── chunks/
                ├── 0001.dcd                         # First committed trajectory chunk
                └── 0001.log                         # Matching committed thermodynamic log
```

</details>

Each replica has its own folder. Committed chunks are numbered `0001`, `0002`,
and so on; `.tmp` files belong to the chunk currently being written. Replicas
`r2` and `r3` use the same layout when run.

#### Simulation Parameters & Frame Saving Cadence

Production dynamics uses standard biophysical integration parameters optimized for cloud stability and storage efficiency:

| Parameter | Value | Biophysical Rationale & Details |
| :--- | :--- | :--- |
| **Timestep ($\Delta t$)** | **$2\text{ fs}$ ($0.002\text{ ps}$)** | Standard Langevin middle integrator with SHAKE/SETTLE constraints on all bonds involving hydrogen. |
| **Frame Saving Interval** | **Every $10\text{ ps}$ ($5,000\text{ steps}$)** | Writes coordinates to `.dcd` trajectory every 10 ps (`--traj-interval 10`), yielding high temporal resolution. |
| **Thermodynamic Log Interval** | **Every $10\text{ ps}$ ($5,000\text{ steps}$)** | Records instantaneous temperature, potential energy, kinetic energy, volume, and density to `.log`. |
| **Chunk Partitioning** | **$1\text{ ns}$ ($500,000\text{ steps}$)** | Simulations are partitioned into $1\text{ ns}$ self-contained chunks ($100\text{ frames}$ per chunk file). |
| **Checkpoint Cadence** | **Every $1\text{ ns}$** | Serializes binary state checkpoint `prod.chk` to Google Drive at the end of each chunk. |
| **Total Frames Generated** | **$10,000\text{ frames}$** | $\frac{100,000\text{ ps}}{10\text{ ps}} = 10,000\text{ frames}$ per $100\text{ ns}$ replica. |
| **Storage Footprint** | **~3.6 GB $\to$ <200 MB** | Raw solvated chunks total ~3.6 GB; stripped trajectories (`cmda merge`) reduce this by >90% to <200 MB. |

#### Downstream Handover

> **Downstream Handover:** These 100 raw 1-ns chunk files are handed over to the merge stage (`cmda merge`), which joins them, eliminates periodic boundary condition (PBC) box jumps, and strips solvent to produce `traj.dcd` for statistical analysis.

#### 💡 How to Resume After a Google Colab Timeout
Production MD writes periodic checkpoints (`prod.chk`) directly to Google Drive. If your session expires:

1. Re-launch the session and re-mount Drive from your local terminal:
   ```bash
   colab new -s cmda --gpu T4
   colab drivemount --session cmda
   colab ssh -s cmda
   ```
2. In the Colab terminal, re-activate the environment:
   ```bash
   [ ! -d "/root/miniforge3/envs/cmda" ] && curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/main/scripts/setup.sh | bash
   source /root/miniforge3/etc/profile.d/conda.sh
   conda activate cmda
   cd /content/drive/MyDrive/research/projects/cmda-example
   ```
3. Re-run the exact same `cmda sim md` command:
   ```bash
   cmda sim md --state apo --variant WT --replica r1 --total-ns 100
   ```
   *ColabMDA detects the existing `prod.chk` checkpoint on Google Drive and resumes dynamics from the last saved timestep.*

#### Trajectory Integrity Auditing (`cmda status`)
To verify chunk continuity, frame counts, and resume readiness before or after merging:
```bash
cmda status --state apo --variant WT --replica r1
```

<details>
<summary>🔍 <b>Click to view example status report</b></summary>

```text
[STATUS]
  Workdir          : data/apo/simulation/WT/md/r1
  Chunks (DCD/log) : 100 / 100
  Topology File    : solvated.pdb
                     └─ 27273 atoms, 8388 residues
                     └─ (169 protein, 8170 water, 49 ions)
  Trajectory File  : traj.dcd (10000 frames)
  Log File         : md.log (10000 frames)
  Frames (from logs): 10000
  Sim time (ps/ns) : 100000.00 ps / 100.0000 ns (from chunk names)
  Resume-ready     : YES (needs prod.chk, system.xml, solvated.pdb)
```

</details>

Alternatively, query frame counts directly in Python with MDAnalysis:
```bash
python3 -c "import MDAnalysis as mda; u = mda.Universe('data/apo/simulation/WT/md/r1/top.prmtop', 'data/apo/simulation/WT/md/r1/traj.dcd'); print('Frames:', len(u.trajectory))"
```

---

### 3.4. Merge Trajectories & Autoimage PBC
Merge raw 1-ns chunked trajectories into a single continuous file, re-center the protein in the periodic unit cell (`autoimage`), and strip bulk solvent:

```bash
# Merge Wild-Type replicas
cmda merge --state apo --variant WT

# Merge Mutant replicas
cmda merge --state apo --variant L119R
```

<details>
<summary><b>View stage directory tree</b></summary>

```text
data/apo/simulation/
├── WT/
│   └── md/
│       └── r1/
│           ├── traj.dcd                             # Autoimaged trajectory without bulk solvent
│           ├── top.prmtop                           # Matching topology; use with the trajectory above
│           ├── trajCa.dcd                           # Cα trajectory
│           ├── topCa.prmtop                         # Matching Cα topology
│           └── md.log                               # Copy of the merged production log
└── L119R/
    └── md/
        └── r1/
            ├── traj.dcd                             # Autoimaged trajectory without bulk solvent
            ├── top.prmtop                           # Matching topology; use with the trajectory above
            ├── trajCa.dcd                           # Cα trajectory
            ├── topCa.prmtop                         # Matching Cα topology
            └── md.log                               # Copy of the merged production log
```

Stage logs: `logs/apo/simulation/WT/md/r1/merge.log`, `logs/apo/simulation/L119R/md/r1/merge.log`.

</details>

---

### 3.5. Trajectory Analysis & Direct Overlay Comparison
Run biophysical trajectory analysis and generate publication-ready comparative overlay plots (WT vs. Mutant):

```bash
# 1. Single-variant analysis (RMSD, Radius of Gyration, RMSF)
cmda analyze --state apo --variant WT
cmda analyze --state apo --variant L119R
```

<details>
<summary><b>View stage directory tree</b></summary>

```text
data/apo/analysis/
├── WT/
│   └── r1/
│       ├── rmsd.csv                                 # RMSD values
│       ├── rg.csv                                   # Radius of gyration values
│       ├── rmsf.csv                                 # Per-residue fluctuations
│       └── figures/
│           ├── rmsd_vs_time.png                     # RMSD plot
│           ├── rg_vs_time.png                       # Radius of gyration plot
│           └── rmsf_per_residue.png                 # Fluctuation plot
└── L119R/
    └── r1/
        ├── rmsd.csv                                 # RMSD values
        ├── rg.csv                                   # Radius of gyration values
        ├── rmsf.csv                                 # Per-residue fluctuations
        └── figures/
            ├── rmsd_vs_time.png                     # RMSD plot
            ├── rg_vs_time.png                       # Radius of gyration plot
            └── rmsf_per_residue.png                 # Fluctuation plot
```

Stage logs: `logs/apo/analysis/WT/r1/analysis.log`, `logs/apo/analysis/L119R/r1/analysis.log`.

</details>

```bash
# 2. Direct comparative overlay (WT vs L119R with Mean ± SD band)
cmda compare --state apo --variants WT,L119R
```

<details>
<summary><b>View stage directory tree</b></summary>

```text
data/apo/analysis/compare/
├── rmsd.png                                         # Variant RMSD comparison
├── rg.png                                           # Variant compactness comparison
├── rmsf.png                                         # Variant fluctuation comparison
├── figures/
│   ├── rmsd_avg.png                                 # Copy of RMSD comparison
│   ├── rg_avg.png                                   # Copy of compactness comparison
│   └── rmsf_avg.png                                 # Copy of fluctuation comparison
└── summary.txt                                      # Comparison statistics
```

</details>

Outputs saved directly in `data/apo/analysis/compare/`:
* `summary.txt` (Executive summary table with Mean ± SD for RMSD, Rg, Temp)
* `rmsd.png` (300 DPI comparative time-series with Mean ± SD band)
* `rg.png` (Radius of Gyration compactness comparison)
* `rmsf.png` (Residue-by-residue fluctuation profiles highlighting loop flexibility)

---

## 📥 Phase 4: Retrieve Results in ONE Command & PyMOL Visualization
*Switch back to your local workstation terminal.*

### 4.1. Pull Analyzed Data from Cloud
Download all analyzed figures, summary tables, and reports:

```bash
cmda pull -s cmda

# To also pull post-processed stripped trajectories for PyMOL inspection:
# cmda pull -s cmda --traj
```

<details>
<summary><b>View transferred directory tree</b></summary>

```text
data/apo/
├── simulation/                                      # Retrieved simulation deliverables
└── analysis/                                        # Retrieved QC, numerical tables and figures
```

</details>

> [!TIP]
> **What `cmda pull` Does Automatically:**
> 1. Packages analysis CSVs, logs, summary manifests, and 300 DPI PNG plots from Google Drive.
> 2. Excludes heavy raw trajectory files (`chunks/*.dcd`, `.chk`) by default so the download finishes in seconds. Pass `--traj` to include stripped trajectories (`traj.dcd`), or `--all` for full raw simulation chunks.
> 3. Extracts everything non-destructively into your local `data/` folder.
> 4. Cleans up temporary cloud archives.

### 4.2. In-Context Reference: Topology Pairing Guidelines

When visualizing or analyzing trajectories, the reference topology file used must match the coordinates in the trajectory file:

| Trajectory File | Description | Matching Topology File | Atom Count | Intended Usage |
| :--- | :--- | :--- | :--- | :--- |
| `prod_raw.dcd` | Raw chunked coordinates | `prep/solvated.pdb` or `system.prmtop` | ~27,000 atoms | Raw simulation checkpointing & resume. |
| `traj.dcd` | Merged, PBC autoimaged, water-stripped | `md/r1/top.prmtop` | ~2,600 atoms | Standard biophysical analysis (`cmda analyze`, `cmda compare`). |
| `trajCa.dcd` | Merged $C_\alpha$ backbone only | `md/r1/topCa.prmtop` | ~300 atoms | Rapid PCA and essential dynamics. |

> [!WARNING]
> **Critical PyMOL Topology Warning:**
> If you load `traj.dcd` (~2,600 atoms) with `solvated.pdb` (~27,000 atoms) in PyMOL or MDAnalysis, the software will crash with a **fatal atom mismatch error**. Always pair stripped trajectories with `top.prmtop`.

### 4.3. Visualizing Trajectories in Local PyMOL

Run these commands in the project directory on a computer with desktop PyMOL.
For a Colab run, first download the processed trajectories using `cmda pull --traj`.
Open one variant at a time; close its PyMOL window before running the next command.

```bash
cmda view --state apo --variant WT --replica r1 --frame all
cmda view --state apo --variant L119R --replica r1 --frame all
```

PyMOL opens and starts playing all saved trajectory frames. Use its playback
controls to pause, play or step through frames. `--frame all` displays the movie;
`--frame 1` displays only the first frame.

The command reads the matching `md/r1/top.prmtop` and `md/r1/traj.dcd` for the
selected state and variant. It generates `view.pdb` and `visualize.pml` in the replica's analysis folder.
Display provenance and verified completion are recorded in the matching
`logs/` folder as `visualize.log`. The PDB has the same atoms as
the processed trajectory, so loading is handled for you.

<details>
<summary><b>Manual loading in PyMOL</b></summary>

After the command has generated the viewing files, these commands can also be
entered in PyMOL from the project directory:

```pml
load data/apo/analysis/WT/r1/view.pdb, protein
load_traj data/apo/simulation/WT/md/r1/traj.dcd, protein, state=1
hide everything
show cartoon, protein
mplay
```

For another variant or replica, change both paths together. Use `solvated.pdb`
only with the original full-system chunks, whose atom count includes water and
salt. The processed `traj.dcd` matches `top.prmtop` and the generated `view.pdb`.

</details>

<details>
<summary><b>Existing runs with older filenames</b></summary>

After production has finished, rename existing processed files without rerunning
or rewriting the trajectories:

```bash
cmda merge --state apo --variant WT --replica r1 --rename-files
```

Then move legacy completion records and visual provenance into `logs/`:

```bash
cmda logs --migrate
```

Repeat the rename command for each completed variant and replica. It validates matching
atom and frame counts before renaming, rejects conflicting files, and leaves raw
chunks, equilibration states and restart checkpoints in place. Stage completion
and metadata live in verified logs; no separate completion flags are created. Older names remain
readable when working with historical runs.

</details>

## Complete project directory structure

<details>
<summary><b>View the complete workflow directory tree</b></summary>

The tree covers the workflow outputs; MODELLER and preparation tools may also
leave auxiliary files in their owning folders. Raw archives may use `.pdb`
instead of `.cif` when the archive fallback is used. The `r2`, `r3` and
`aggregate` entries apply after all three replicas have been run and analyzed.

```text
~/research/projects/cmda-example/
├── data/
│   └── apo/
│       ├── structure/
│       │   └── O76039_4BGQ/
│       │       ├── raw/
│       │       │   ├── 4bgq.cif                     # Experimental archive; PDB fallback is possible
│       │       │   ├── template.pdb                 # Cleaned template coordinates
│       │       │   └── uniprot.fasta                # Sequence used for modeling
│       │       └── variants/
│       │           ├── WT/
│       │           │   └── WT.pdb                   # Coordinates passed to preparation
│       │           └── L119R/
│       │               └── L119R.pdb                # Coordinates passed to preparation
│       ├── simulation/
│       │   ├── WT/
│       │   │   ├── prep/
│       │   │   │   ├── model.pdb                    # Prepared solute coordinates
│       │   │   │   ├── solvated.pdb                 # Protein, water and ions
│       │   │   │   ├── system.prmtop                # Amber topology passed to equilibration
│       │   │   │   └── system.inpcrd                # Amber coordinates passed to equilibration
│       │   │   ├── eq/
│       │   │   │   ├── solvated.pdb                 # Topology coordinates used by OpenMM
│       │   │   │   ├── system.xml                   # Serialized OpenMM system
│       │   │   │   ├── em.state.xml                 # Minimized state
│       │   │   │   ├── nvt.state.xml                # State after NVT equilibration
│       │   │   │   ├── npt.state.xml                # Equilibrated state passed to production
│       │   │   │   ├── eq.pdb                       # Final equilibrated coordinates
│       │   │   │   └── stages/
│       │   │   │       ├── state_01.xml             # First staged restart
│       │   │   │       └── state_24.xml             # Last of 24 staged restart states
│       │   │   └── md/
│       │   │       ├── r1/
│       │   │       │   ├── system.xml               # Copied OpenMM system
│       │   │       │   ├── solvated.pdb             # Copied topology coordinates
│       │   │       │   ├── npt.state.xml            # Starting state copied from equilibration
│       │   │       │   ├── prod.chk                 # Checkpoint used to resume production
│       │   │       │   ├── chunks/
│       │   │       │   │   ├── 0001.dcd             # First committed trajectory chunk
│       │   │       │   │   └── 0001.log             # Matching committed thermodynamic log
│       │   │       │   ├── traj.dcd                 # Autoimaged trajectory without bulk solvent
│       │   │       │   ├── top.prmtop               # Matching topology; use with the trajectory above
│       │   │       │   ├── trajCa.dcd               # Cα trajectory
│       │   │       │   ├── topCa.prmtop             # Matching Cα topology
│       │   │       │   └── md.log                   # Copy of the merged production log
│       │   │       ├── r2/
│       │   │       │   ├── system.xml               # Copied OpenMM system
│       │   │       │   ├── solvated.pdb             # Copied topology coordinates
│       │   │       │   ├── npt.state.xml            # Starting state copied from equilibration
│       │   │       │   ├── prod.chk                 # Checkpoint used to resume production
│       │   │       │   ├── chunks/
│       │   │       │   │   ├── 0001.dcd             # First committed trajectory chunk
│       │   │       │   │   └── 0001.log             # Matching committed thermodynamic log
│       │   │       │   ├── traj.dcd                 # Autoimaged trajectory without bulk solvent
│       │   │       │   ├── top.prmtop               # Matching topology; use with the trajectory above
│       │   │       │   ├── trajCa.dcd               # Cα trajectory
│       │   │       │   ├── topCa.prmtop             # Matching Cα topology
│       │   │       │   └── md.log                   # Copy of the merged production log
│       │   │       └── r3/
│       │   │           ├── system.xml               # Copied OpenMM system
│       │   │           ├── solvated.pdb             # Copied topology coordinates
│       │   │           ├── npt.state.xml            # Starting state copied from equilibration
│       │   │           ├── prod.chk                 # Checkpoint used to resume production
│       │   │           ├── chunks/
│       │   │           │   ├── 0001.dcd             # First committed trajectory chunk
│       │   │           │   └── 0001.log             # Matching committed thermodynamic log
│       │   │           ├── traj.dcd                 # Autoimaged trajectory without bulk solvent
│       │   │           ├── top.prmtop               # Matching topology; use with the trajectory above
│       │   │           ├── trajCa.dcd               # Cα trajectory
│       │   │           ├── topCa.prmtop             # Matching Cα topology
│       │   │           └── md.log                   # Copy of the merged production log
│       │   └── L119R/
│       │       ├── prep/
│       │       │   ├── model.pdb                    # Prepared solute coordinates
│       │       │   ├── solvated.pdb                 # Protein, water and ions
│       │       │   ├── system.prmtop                # Amber topology passed to equilibration
│       │       │   └── system.inpcrd                # Amber coordinates passed to equilibration
│       │       ├── eq/
│       │       │   ├── solvated.pdb                 # Topology coordinates used by OpenMM
│       │       │   ├── system.xml                   # Serialized OpenMM system
│       │       │   ├── em.state.xml                 # Minimized state
│       │       │   ├── nvt.state.xml                # State after NVT equilibration
│       │       │   ├── npt.state.xml                # Equilibrated state passed to production
│       │       │   ├── eq.pdb                       # Final equilibrated coordinates
│       │       │   └── stages/
│       │       │       ├── state_01.xml             # First staged restart
│       │       │       └── state_24.xml             # Last of 24 staged restart states
│       │       └── md/
│       │           ├── r1/
│       │           │   ├── system.xml               # Copied OpenMM system
│       │           │   ├── solvated.pdb             # Copied topology coordinates
│       │           │   ├── npt.state.xml            # Starting state copied from equilibration
│       │           │   ├── prod.chk                 # Checkpoint used to resume production
│       │           │   ├── chunks/
│       │           │   │   ├── 0001.dcd             # First committed trajectory chunk
│       │           │   │   └── 0001.log             # Matching committed thermodynamic log
│       │           │   ├── traj.dcd                 # Autoimaged trajectory without bulk solvent
│       │           │   ├── top.prmtop               # Matching topology; use with the trajectory above
│       │           │   ├── trajCa.dcd               # Cα trajectory
│       │           │   ├── topCa.prmtop             # Matching Cα topology
│       │           │   └── md.log                   # Copy of the merged production log
│       │           ├── r2/
│       │           │   ├── system.xml               # Copied OpenMM system
│       │           │   ├── solvated.pdb             # Copied topology coordinates
│       │           │   ├── npt.state.xml            # Starting state copied from equilibration
│       │           │   ├── prod.chk                 # Checkpoint used to resume production
│       │           │   ├── chunks/
│       │           │   │   ├── 0001.dcd             # First committed trajectory chunk
│       │           │   │   └── 0001.log             # Matching committed thermodynamic log
│       │           │   ├── traj.dcd                 # Autoimaged trajectory without bulk solvent
│       │           │   ├── top.prmtop               # Matching topology; use with the trajectory above
│       │           │   ├── trajCa.dcd               # Cα trajectory
│       │           │   ├── topCa.prmtop             # Matching Cα topology
│       │           │   └── md.log                   # Copy of the merged production log
│       │           └── r3/
│       │               ├── system.xml               # Copied OpenMM system
│       │               ├── solvated.pdb             # Copied topology coordinates
│       │               ├── npt.state.xml            # Starting state copied from equilibration
│       │               ├── prod.chk                 # Checkpoint used to resume production
│       │               ├── chunks/
│       │               │   ├── 0001.dcd             # First committed trajectory chunk
│       │               │   └── 0001.log             # Matching committed thermodynamic log
│       │               ├── traj.dcd                 # Autoimaged trajectory without bulk solvent
│       │               ├── top.prmtop               # Matching topology; use with the trajectory above
│       │               ├── trajCa.dcd               # Cα trajectory
│       │               ├── topCa.prmtop             # Matching Cα topology
│       │               └── md.log                   # Copy of the merged production log
│       └── analysis/
│           ├── WT/
│           │   ├── eq/
│           │   │   ├── qc.txt                       # Equilibration diagnostics
│           │   │   └── qc.png                       # Thermodynamic and structural QC figure
│           │   ├── r1/
│           │   │   ├── rmsd.csv                     # RMSD values
│           │   │   ├── rg.csv                       # Radius of gyration values
│           │   │   ├── rmsf.csv                     # Per-residue fluctuations
│           │   │   └── figures/
│           │   │       ├── rmsd_vs_time.png         # RMSD plot
│           │   │       ├── rg_vs_time.png           # Radius of gyration plot
│           │   │       └── rmsf_per_residue.png     # Fluctuation plot
│           │   ├── r2/
│           │   │   ├── rmsd.csv                     # RMSD values
│           │   │   ├── rg.csv                       # Radius of gyration values
│           │   │   ├── rmsf.csv                     # Per-residue fluctuations
│           │   │   └── figures/
│           │   │       ├── rmsd_vs_time.png         # RMSD plot
│           │   │       ├── rg_vs_time.png           # Radius of gyration plot
│           │   │       └── rmsf_per_residue.png     # Fluctuation plot
│           │   ├── r3/
│           │   │   ├── rmsd.csv                     # RMSD values
│           │   │   ├── rg.csv                       # Radius of gyration values
│           │   │   ├── rmsf.csv                     # Per-residue fluctuations
│           │   │   └── figures/
│           │   │       ├── rmsd_vs_time.png         # RMSD plot
│           │   │       ├── rg_vs_time.png           # Radius of gyration plot
│           │   │       └── rmsf_per_residue.png     # Fluctuation plot
│           │   └── aggregate/
│           │       ├── rmsd_avg.csv                 # Replica mean and standard deviation
│           │       ├── rg_avg.csv                   # Replica mean and standard deviation
│           │       ├── rmsf_avg.csv                 # Replica mean and standard deviation
│           │       └── figures/
│           │           ├── rmsd_aggregate.png       # Replica RMSD summary
│           │           ├── rg_aggregate.png         # Replica compactness summary
│           │           └── rmsf_aggregate.png       # Replica fluctuation summary
│           ├── L119R/
│           │   ├── eq/
│           │   │   ├── qc.txt                       # Equilibration diagnostics
│           │   │   └── qc.png                       # Thermodynamic and structural QC figure
│           │   ├── r1/
│           │   │   ├── rmsd.csv                     # RMSD values
│           │   │   ├── rg.csv                       # Radius of gyration values
│           │   │   ├── rmsf.csv                     # Per-residue fluctuations
│           │   │   └── figures/
│           │   │       ├── rmsd_vs_time.png         # RMSD plot
│           │   │       ├── rg_vs_time.png           # Radius of gyration plot
│           │   │       └── rmsf_per_residue.png     # Fluctuation plot
│           │   ├── r2/
│           │   │   ├── rmsd.csv                     # RMSD values
│           │   │   ├── rg.csv                       # Radius of gyration values
│           │   │   ├── rmsf.csv                     # Per-residue fluctuations
│           │   │   └── figures/
│           │   │       ├── rmsd_vs_time.png         # RMSD plot
│           │   │       ├── rg_vs_time.png           # Radius of gyration plot
│           │   │       └── rmsf_per_residue.png     # Fluctuation plot
│           │   ├── r3/
│           │   │   ├── rmsd.csv                     # RMSD values
│           │   │   ├── rg.csv                       # Radius of gyration values
│           │   │   ├── rmsf.csv                     # Per-residue fluctuations
│           │   │   └── figures/
│           │   │       ├── rmsd_vs_time.png         # RMSD plot
│           │   │       ├── rg_vs_time.png           # Radius of gyration plot
│           │   │       └── rmsf_per_residue.png     # Fluctuation plot
│           │   └── aggregate/
│           │       ├── rmsd_avg.csv                 # Replica mean and standard deviation
│           │       ├── rg_avg.csv                   # Replica mean and standard deviation
│           │       ├── rmsf_avg.csv                 # Replica mean and standard deviation
│           │       └── figures/
│           │           ├── rmsd_aggregate.png       # Replica RMSD summary
│           │           ├── rg_aggregate.png         # Replica compactness summary
│           │           └── rmsf_aggregate.png       # Replica fluctuation summary
│           └── compare/
│               ├── rmsd.png                         # Variant RMSD comparison
│               ├── rg.png                           # Variant compactness comparison
│               ├── rmsf.png                         # Variant fluctuation comparison
│               ├── figures/
│               │   ├── rmsd_avg.png                 # Copy of RMSD comparison
│               │   ├── rg_avg.png                   # Copy of compactness comparison
│               │   └── rmsf_avg.png                 # Copy of fluctuation comparison
│               └── summary.txt                      # Comparison statistics
└── logs/                                            # Execution logs mirror the data folders
    └── apo/
        ├── structure/
        │   └── O76039_4BGQ/
        │       ├── structure.log                    # Stage metadata and verified completion
        │       └── variants/
        │           ├── WT/
        │           │   └── model.log                # Stage metadata and verified completion
        │           └── L119R/
        │               └── model.log                # Stage metadata and verified completion
        ├── simulation/
        │   ├── WT/
        │   │   ├── prep/
        │   │   │   └── prep.log                     # Stage metadata and verified completion
        │   │   └── md/
        │   │       ├── r1/
        │   │       │   └── merge.log                # Stage metadata and verified completion
        │   │       ├── r2/
        │   │       │   └── merge.log                # Stage metadata and verified completion
        │   │       └── r3/
        │   │           └── merge.log                # Stage metadata and verified completion
        │   └── L119R/
        │       ├── prep/
        │       │   └── prep.log                     # Stage metadata and verified completion
        │       └── md/
        │           ├── r1/
        │           │   └── merge.log                # Stage metadata and verified completion
        │           ├── r2/
        │           │   └── merge.log                # Stage metadata and verified completion
        │           └── r3/
        │               └── merge.log                # Stage metadata and verified completion
        └── analysis/
            ├── WT/
            │   ├── r1/
            │   │   ├── analysis.log                 # Stage metadata and verified completion
            │   │   └── visualize.log                # Movie inputs and display provenance
            │   ├── r2/
            │   │   ├── analysis.log                 # Stage metadata and verified completion
            │   │   └── visualize.log                # Movie inputs and display provenance
            │   ├── r3/
            │   │   ├── analysis.log                 # Stage metadata and verified completion
            │   │   └── visualize.log                # Movie inputs and display provenance
            │   └── aggregate/
            │       └── aggregate.log                # Stage metadata and verified completion
            └── L119R/
                ├── r1/
                │   ├── analysis.log                 # Stage metadata and verified completion
                │   └── visualize.log                # Movie inputs and display provenance
                ├── r2/
                │   ├── analysis.log                 # Stage metadata and verified completion
                │   └── visualize.log                # Movie inputs and display provenance
                ├── r3/
                │   ├── analysis.log                 # Stage metadata and verified completion
                │   └── visualize.log                # Movie inputs and display provenance
                └── aggregate/
                    └── aggregate.log                # Stage metadata and verified completion
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
