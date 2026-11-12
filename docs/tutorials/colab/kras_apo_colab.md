# KRAS (4LDJ) Apo State: Google Colab Walkthrough

This tutorial provides a streamlined simulation walkthrough for the **Apo state** of the **KRAS GTPase** (UniProt **P01116**, PDB **4LDJ** Chain A, residues **1–169**), comparing the Wild-Type (WT) protein against the oncogenic point mutant **G12D**.

Everything is executed on a **Google Colab cloud GPU via the `colab` CLI**, with all outputs, checkpoints, and figures continuously persisted inside your **Google Drive**.

---

Outputs are shown after each stage. The complete directory tree is collapsed at
the end of this tutorial. Paths are relative to the project root unless stated otherwise.

---

## 1. Connect to Google Colab

### 1.1. Local Setup
Run this one-time setup in your local terminal on **Linux (Ubuntu)** or Windows via **[WSL2](https://learn.microsoft.com/en-us/windows/wsl/install)**:

```bash
# 1. Install required system tools (Ubuntu/WSL2):
sudo apt-get update
sudo apt-get install -y git curl wget ca-certificates

# 2. Install cmda, OpenMM, AmberTools, PyMOL, and google-colab-cli:
bash -c "$(curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/main/scripts/setup.sh)"

# 3. Activate the environment in this terminal:
source ~/miniforge3/etc/profile.d/conda.sh
conda activate cmda

# 4. Log in to Google Colab (one-time browser authorization):
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

### 1.2. Launch Cloud GPU
Once logged in, provision a Google Colab GPU instance and enter the remote console:

```bash
# 1. Launch a Colab GPU instance:
colab new -s cmda --gpu T4

# 2. Mount Google Drive (visit the printed link once, then press Enter):
colab drivemount --session cmda

# 3. Enter the Colab Console:
colab console -s cmda
```

---

## 2. Cloud Environment Setup

Inside the active **Colab console terminal**, provide your Modeller license key, run the bootstrap installer, and navigate to Google Drive:

```bash
# 1. Export Modeller key & run automated installer:
export KEY_MODELLER="YOUR_MODELLER_KEY"
curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/main/scripts/setup.sh | bash

# 2. Activate environment:
source /root/miniforge3/etc/profile.d/conda.sh
conda activate cmda

# 3. Navigate to persistent project folder on Google Drive:
mkdir -p /content/drive/MyDrive/research/projects/kras-example
cd /content/drive/MyDrive/research/projects/kras-example
```

> [!NOTE]
> Setting `KEY_MODELLER` tells `setup.sh` to install Modeller from the official `salilab` channel and inject your license. Academic keys are free at [salilab.org/modeller/registration.html](https://salilab.org/modeller/registration.html).

---

## 3. Modeling & Mutations

Fetch experimental coordinates from the PDB archive, align to canonical UniProt indexing, and model the full WT domain and G12D point mutation:

```bash
# 1. Download 4LDJ (Chain A, residues 1–169) & verify canonical mapping:
cmda fetch --state apo --pdb-id 4LDJ --uniprot-id P01116 --chain A --range 1 169
```

<details>
<summary><b>View stage directory tree</b></summary>

```text
data/apo/structure/P01116_4LDJ/
├── raw/
│   ├── 4ldj.cif                                     # Experimental archive; PDB fallback is possible
│   ├── template.pdb                                 # Cleaned template coordinates
│   └── uniprot.fasta                                # Sequence used for modeling
└── STRUCTURE_DONE                                   # Completed fetch and mapping
```

</details>

```bash
# 2. Model canonical Wild-Type coordinates (closes missing crystallographic loops):
cmda model --state apo --variant WT

# 3. Build Point Mutant (G12D) on the verified WT backbone:
cmda model --state apo --variant G12D
```

<details>
<summary><b>View stage directory tree</b></summary>

```text
data/apo/structure/P01116_4LDJ/variants/
├── WT/
│   ├── WT.pdb                                       # Coordinates passed to preparation
│   └── MODEL_DONE                                   # Completed model verification
└── G12D/
    ├── G12D.pdb                                     # Coordinates passed to preparation
    └── MODEL_DONE                                   # Completed model verification
```

</details>

> **Downstream Handover:** These validated coordinate models (`WT.pdb`, `G12D.pdb`) are handed over to **Section 4 (`cmda prep`)** for PROPKA protonation and Amber OPC water box solvation.

---

## 4. Structure Preparation & Solvation

ColabMDA runs **PROPKA 3** to calculate 3D electrostatic microenvironments, assign explicit AMBER titratable protonation states, solvate with OPC 4-point water, and neutralize net charge:

```bash
cmda prep --state apo --variant WT
cmda prep --state apo --variant G12D
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
│       ├── system.inpcrd                            # Amber coordinates passed to equilibration
│       └── PREP_DONE                                # Verified preparation outputs
└── G12D/
    └── prep/
        ├── model.pdb                                # Prepared solute coordinates
        ├── solvated.pdb                             # Protein, water and ions
        ├── system.prmtop                            # Amber topology passed to equilibration
        ├── system.inpcrd                            # Amber coordinates passed to equilibration
        └── PREP_DONE                                # Verified preparation outputs
```

</details>

| Residue | State | Physical Criteria |
| :--- | :--- | :--- |
| **Histidine** | `HID` / `HIE` / `HIP` | Neutral $\text{N}_\delta$, neutral $\text{N}_\epsilon$ (default), or positively charged doubly protonated ($pK_a > 7.4$). |
| **Acidic** | `ASP` / `GLU` | Standard carboxylate anions (`ASH`/`GLH` assigned if buried with $pK_a > 7.4$). |
| **Basic** | `LYS` / `ARG` | Positively charged physiological cations (`LYN` if neutral). |
| **Cysteine** | `CYS` / `CYX` | Free neutral thiol (`CYS`) or disulfide bridge (`CYX`). |

> **Downstream Handover:** These solvated parameters are loaded by OpenMM for GPU equilibration in **Section 5 (`cmda sim eq`)**.

---

## 5. Equilibration & QC

Initialize simulation directories and run the progressive 24-stage equilibration relaxation ladder. This heats the system from 50 K to 300 K and gradually releases positional restraints:

```bash
# 1. Initialize simulation directories:
cmda stage --state apo --variant WT
cmda stage --state apo --variant G12D

# 2. Execute 24-stage progressive equilibration (EM -> NVT -> NPT -> QC):
cmda sim eq --state apo --variant WT
cmda sim eq --state apo --variant G12D
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
│   └── G12D/
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
    └── G12D/
        └── eq/
            ├── qc.txt                               # Equilibration diagnostics
            └── qc.png                               # Thermodynamic and structural QC figure
```

</details>

### Stage Schedule (~1.96 ns Total)

| Stages | Operation | Duration | Timestep & Restraints |
| :--- | :--- | :--- | :--- |
| **01–03** | Restrained Minimization | 5,000 steps each | $k = 100, 50, 25\text{ kcal/mol/\AA}^2$ |
| **04–09** | NVT Progressive Heating ($50\text{ K} \to 300\text{ K}$) | **$6 \times 10\text{ ps}$ ($60\text{ ps}$)** | $1\text{ fs timestep}$, $k = 25$ |
| **10–14** | Post-Heating Minimization | 5,000 / 1,000 steps | $k = 25 \to 1\text{ kcal/mol/\AA}^2$ |
| **15–18** | NVT Restraint Release at $300\text{ K}$ | **$4 \times 100\text{ ps}$ ($400\text{ ps}$)** | $2\text{ fs timestep}$, $k = 10, 5, 2, 1$ |
| **19–22** | NPT Density Relaxation ($300\text{ K}, 1\text{ bar}$) | **$4 \times 100\text{ ps}$ ($400\text{ ps}$)** | $2\text{ fs timestep}$, $k = 10, 5, 2, 1$ |
| **23** | Unrestrained NPT Relaxation | **$100\text{ ps}$** | $2\text{ fs timestep}$, unrestrained ($k = 0$) |
| **24** | Final Unrestrained NPT Stabilization | **$1,000\text{ ps}$ ($1.0\text{ ns}$)** | $2\text{ fs timestep}$, unrestrained ($k = 0$) |

> [!TIP]
> **Fast Smoke Test Alternative:** For a 2-minute quick check rather than full 24-stage production relaxation, append `--protocol quick` (runs 100 ps NVT + 100 ps NPT).

### QC Verification
```bash
cat data/apo/analysis/WT/eq/qc.txt
cat data/apo/analysis/G12D/eq/qc.txt
```
*Expected Output:* `Status: PASSED`, `Cα RMSD: ~0.83 Å` (Confirms temperature $\approx 300\text{ K}$, density $\approx 1.01\text{ g/mL}$).

> **Downstream Handover:** The final thermodynamic state (`eq/npt.state.xml`) and topology serve as the physical starting point for independent production replicas in **Section 6 (`cmda sim md`)**, forking with distinct Maxwell-Boltzmann velocity seeds.

---

## 6. Production MD (100 ns)

### Why Triplicates?
Molecular dynamics is inherently stochastic. A single trajectory can easily become trapped in a metastable conformational sub-state. Running **triplicates** ($N=3$ independent trajectories: `r1`, `r2`, `r3`) ensures conformational sampling convergence, allowing you to compute true **Mean ± SD** statistical error bands for RMSD, Radius of Gyration, and residue fluctuations.

> 💡 **Why Shared `eq/`?**
> Equilibration only needs to be run **once** per system. Replicas (`r1`, `r2`, `r3`) fork directly from the final relaxed coordinates in `eq/npt.state.xml` using distinct Maxwell-Boltzmann random velocity seeds. This eliminates redundant equilibration, saving hours of GPU time on Google Colab while guaranteeing strict statistical independence.

### 6.1. Launch Production Runs (`r1`)
Launch production dynamics for the primary replica (`r1`) of Wild-Type and G12D:

```bash
# 1. Run Wild-Type replica 1 (100 ns):
cmda sim md --state apo --variant WT --replica r1 --total-ns 100

# 2. Run G12D mutant replica 1 (100 ns):
cmda sim md --state apo --variant G12D --replica r1 --total-ns 100
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
└── G12D/
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

> [!IMPORTANT]
> **💡 How to Resume After a Session Disconnection or Timeout:**
> Google Colab GPU sessions have runtime limits (typically 4–12 hours) before the instance is recycled. Because a recycled Colab instance starts with a fresh filesystem, **you must reinstall the environment on the new instance**. Your simulation progress and checkpoints remain safe on Google Drive.
>
> If your session times out or disconnects:
> 1. **Launch a new session & mount Drive (Local Terminal):**
>    ```bash
>    colab new -s cmda --gpu T4
>    colab drivemount --session cmda
>    colab console -s cmda
>    ```
> 2. **Reinstall environment & navigate to Drive (Colab Console, ~90 seconds):**
>    ```bash
>    curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/main/scripts/setup.sh | bash
>    source /root/miniforge3/etc/profile.d/conda.sh
>    conda activate cmda
>    cd /content/drive/MyDrive/research/projects/kras-example
>    ```
> 3. **Re-run the exact same command:**
>    ```bash
>    # Re-run whichever variant was running:
>    cmda sim md --state apo --variant WT --replica r1 --total-ns 100
>    # or for mutant:
>    cmda sim md --state apo --variant G12D --replica r1 --total-ns 100
>    ```
> ColabMDA automatically detects `prod.chk` and existing completed 1-ns trajectory chunks on Google Drive, and **resumes dynamics from the exact last saved timestep without losing any progress**.

### 6.2. Simulation Parameters & Frame Saving Cadence

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

### 6.3. Downstream Handover

> **Downstream Handover:** These 100 raw 1-ns chunk files are handed over to **Section 7 (`cmda merge`)**, which joins them, eliminates periodic boundary condition (PBC) box jumps, and strips solvent to produce `prod_keep_ligand.dcd` for statistical analysis.

<details>
<summary><b>💡 (Optional) How to Launch Additional Replicas (<code>r2</code>, <code>r3</code>) Later</b></summary>

You do **not** need to run these now—proceed directly to **Section 7 (Merge & Analysis)** using `r1`.

When you are ready to produce final publication error bands (Mean ± SD), you can launch the remaining independent replicas:

```bash
# Additional Wild-Type replicas (run when ready):
# cmda sim md --state apo --variant WT --replica r2 --total-ns 100
# cmda sim md --state apo --variant WT --replica r3 --total-ns 100

# Additional G12D mutant replicas (run when ready):
# cmda sim md --state apo --variant G12D --replica r2 --total-ns 100
# cmda sim md --state apo --variant G12D --replica r3 --total-ns 100
```

</details>

---

## 7. Merge & Analysis

Concatenate raw 1-ns chunks, apply periodic boundary condition (PBC) autoimaging, strip bulk water, audit frame continuity, and generate comparative figures:

```bash
# 1. Merge chunks & autoimage PBC (produces autoimaged, stripped trajectories):
cmda merge --state apo --variant WT
cmda merge --state apo --variant G12D
```

<details>
<summary><b>View stage directory tree</b></summary>

```text
data/apo/simulation/
├── WT/
│   └── md/
│       └── r1/
│           ├── prod_keep_ligand.dcd                 # Autoimaged trajectory without bulk solvent
│           ├── stripped.prmtop                      # Matching topology; use with the trajectory above
│           ├── prod_ca.dcd                          # Cα trajectory
│           ├── ca.prmtop                            # Matching Cα topology
│           ├── prod_full.log                        # Merged production log
│           ├── md.log                               # Copy of the merged production log
│           └── MERGE_DONE                           # Verified merged outputs
└── G12D/
    └── md/
        └── r1/
            ├── prod_keep_ligand.dcd                 # Autoimaged trajectory without bulk solvent
            ├── stripped.prmtop                      # Matching topology; use with the trajectory above
            ├── prod_ca.dcd                          # Cα trajectory
            ├── ca.prmtop                            # Matching Cα topology
            ├── prod_full.log                        # Merged production log
            ├── md.log                               # Copy of the merged production log
            └── MERGE_DONE                           # Verified merged outputs
```

</details>

```bash
# 2. Audit frame counts & chunk continuity:
cmda status --state apo --variant WT --replica r1
cmda status --state apo --variant G12D --replica r1
```

```bash
# 3. Calculate single-variant biophysical metrics (RMSD, Rg, RMSF):
cmda analyze --state apo --variant WT
cmda analyze --state apo --variant G12D
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
│       ├── figures/
│       │   ├── rmsd_vs_time.png                     # RMSD plot
│       │   ├── rg_vs_time.png                       # Radius of gyration plot
│       │   └── rmsf_per_residue.png                 # Fluctuation plot
│       └── ANALYSIS_DONE                            # Verified replica tables and figures
└── G12D/
    └── r1/
        ├── rmsd.csv                                 # RMSD values
        ├── rg.csv                                   # Radius of gyration values
        ├── rmsf.csv                                 # Per-residue fluctuations
        ├── figures/
        │   ├── rmsd_vs_time.png                     # RMSD plot
        │   ├── rg_vs_time.png                       # Radius of gyration plot
        │   └── rmsf_per_residue.png                 # Fluctuation plot
        └── ANALYSIS_DONE                            # Verified replica tables and figures
```

</details>

```bash
# 4. Generate direct comparative overlay figures (WT vs G12D with Mean ± SD bands):
cmda compare --state apo --variants WT,G12D
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

> **Downstream Handover:** These publication figures and statistical summaries are retrieved to your local machine in **Section 8 (`cmda pull`)** for final inspection and publication rendering.

---

## 8. Pull Results & Visualization

Exit the Colab console (`exit`), then switch back to your local workstation terminal:

### 8.1. Pull Data
Download analyzed figures, summary tables, and trajectories:

```bash
mkdir -p ~/research/projects/kras-example
cd ~/research/projects/kras-example

# Pull lightweight figures, logs, and summary tables (~2 MB):
cmda pull -s cmda

# To also pull merged trajectories for local PyMOL inspection:
cmda pull -s cmda --traj
```

<details>
<summary><b>View transferred directory tree</b></summary>

```text
data/apo/
├── simulation/             # Retrieved simulation deliverables
└── analysis/               # Retrieved QC, numerical tables and figures
```

</details>

> [!TIP]
> **What `cmda pull` Does:**
> 1. Packages analysis CSVs, logs, summary manifests, and PNG plots from Google Drive.
> 2. Excludes heavy raw simulation checkpoints (`chunks/*.dcd`, `.chk`) by default so the download completes quickly. Passing `--traj` pulls lightweight solvent-stripped trajectories (`prod_keep_ligand.dcd`), while `--all` downloads full raw chunks.
> 3. Extracts files non-destructively into your local `data/` folder.

> **Downstream Handover:** With deliverables pulled to your local machine, proceed to **Section 8.2** for PyMOL trajectory inspection and structural verification.

---

### 8.2. PyMOL Trajectory Visualization

To inspect trajectories, open the stripped topology (`stripped.prmtop`) paired with the stripped trajectory (`prod_keep_ligand.dcd`).

#### Option A: Quick Launch via CLI
The built-in `cmda view` command automatically generates a PyMOL script with backbone fitting and launches the GUI:

```bash
# View Wild-Type trajectory:
cmda view --state apo --variant WT --replica r1

# View Mutant trajectory (highlighting mutated residue 12):
cmda view --state apo --variant G12D --replica r1 --resi 12
```

#### Option B: Launch Directly with Native PyMOL
Open the topology and trajectory directly from your terminal:

```bash
pymol data/apo/simulation/WT/md/r1/stripped.prmtop data/apo/simulation/WT/md/r1/prod_keep_ligand.dcd
```

#### Option C: Interactive PyMOL GUI (Desktop App)
If opening PyMOL from your desktop application:
1. **Load Topology:** Go to **File → Open...** (or drag and drop) and select `data/apo/simulation/WT/md/r1/stripped.prmtop` (or `WT.pdb`).
2. **Load Trajectory:** Go to **File → Open...** and select `data/apo/simulation/WT/md/r1/prod_keep_ligand.dcd`. Ensure **Load into object** is set to the already-loaded protein object.
3. **Play & Format:** Use the right-hand panel buttons (`S` Show → `as cartoon`, `C` Color → `cyan`) and click `▶` on the bottom-right movie controls to play.

#### Inside PyMOL (Quick Commands):
Enter these commands in the PyMOL command bar to clean up and play the simulation:

```text
# 1. Fit backbone frame-by-frame (eliminates tumbling & drift):
intra_fit all and name CA

# 2. Show clean secondary structure cartoon:
show cartoon, all
color cyan, all and name C*

# 3. Highlight residue 12 as yellow sticks:
show sticks, resi 12 and not name N+C+O+H
color yellow, resi 12

# 4. Play trajectory:
# Use the ▶ button in the bottom-right movie controls, or Left/Right arrow keys.
```

#### Render Comparative Snapshot Grids:
Render 3×8 transition snapshot grids comparing structural conformations across time:

```bash
cmda snapshots --state apo
```

<details>
<summary><b>Click to view PyMOL visualization script (<code>visualize.pml</code>)</b></summary>

```text
# 1. Clear out memory
reinitialize
bg_color white

# 2. Load stripped topology and trajectory
load stripped.prmtop, kras
load_traj prod_keep_ligand.dcd, kras

# 3. Fit backbone frame-by-frame
intra_fit kras and name CA

# 4. Generate secondary structure cartoon
dss kras
cartoon automatic, kras
show cartoon, kras

# 5. Apply cyan color scheme
color cyan, kras and name C*
util.cnc("kras")

# 6. Highlight mutated residue 12 side chain sticks
show sticks, kras and resi 12 and not name N+C+O+H
color yellow, kras and resi 12
set stick_radius, 0.25

# 7. Focus camera onto protein
zoom kras and polymer, buffer=4
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

## Complete project directory structure

<details>
<summary><b>View the complete workflow directory tree</b></summary>

The tree covers the workflow outputs; MODELLER and preparation tools may also
leave auxiliary files in their owning folders. Raw archives may use `.pdb`
instead of `.cif` when the archive fallback is used. The `r2`, `r3` and
`aggregate` entries apply after all three replicas have been run and analyzed.

```text
/content/drive/MyDrive/research/projects/kras-example/
└── data/
    └── apo/
        ├── structure/
        │   └── P01116_4LDJ/
        │       ├── raw/
        │       │   ├── 4ldj.cif                     # Experimental archive; PDB fallback is possible
        │       │   ├── template.pdb                 # Cleaned template coordinates
        │       │   └── uniprot.fasta                # Sequence used for modeling
        │       ├── STRUCTURE_DONE                   # Completed fetch and mapping
        │       └── variants/
        │           ├── WT/
        │           │   ├── WT.pdb                   # Coordinates passed to preparation
        │           │   └── MODEL_DONE               # Completed model verification
        │           └── G12D/
        │               ├── G12D.pdb                 # Coordinates passed to preparation
        │               └── MODEL_DONE               # Completed model verification
        ├── simulation/
        │   ├── WT/
        │   │   ├── prep/
        │   │   │   ├── model.pdb                    # Prepared solute coordinates
        │   │   │   ├── solvated.pdb                 # Protein, water and ions
        │   │   │   ├── system.prmtop                # Amber topology passed to equilibration
        │   │   │   ├── system.inpcrd                # Amber coordinates passed to equilibration
        │   │   │   └── PREP_DONE                    # Verified preparation outputs
        │   │   ├── eq/
        │   │   │   ├── solvated.pdb                 # Topology coordinates used by OpenMM
        │   │   │   ├── system.xml                   # Serialized OpenMM system
        │   │   │   ├── em.state.xml                 # Minimized state
        │   │   │   ├── nvt.state.xml                # State after NVT equilibration
        │   │   │   ├── npt.state.xml                # Equilibrated state passed to production
        │   │   │   ├── eq.pdb                       # Final equilibrated coordinates
        │   │   │   └── stages/
        │   │   │       ├── state_01.xml             # First staged restart
        │   │   │       └── state_24.xml             # Last of 24 staged restart states
        │   │   └── md/
        │   │       ├── r1/
        │   │       │   ├── system.xml               # Copied OpenMM system
        │   │       │   ├── solvated.pdb             # Copied topology coordinates
        │   │       │   ├── npt.state.xml            # Starting state copied from equilibration
        │   │       │   ├── prod.chk                 # Checkpoint used to resume production
        │   │       │   ├── chunks/
        │   │       │   │   ├── 0001.dcd             # First committed trajectory chunk
        │   │       │   │   └── 0001.log             # Matching committed thermodynamic log
        │   │       │   ├── prod_keep_ligand.dcd     # Autoimaged trajectory without bulk solvent
        │   │       │   ├── stripped.prmtop          # Matching topology; use with the trajectory above
        │   │       │   ├── prod_ca.dcd              # Cα trajectory
        │   │       │   ├── ca.prmtop                # Matching Cα topology
        │   │       │   ├── prod_full.log            # Merged production log
        │   │       │   ├── md.log                   # Copy of the merged production log
        │   │       │   └── MERGE_DONE               # Verified merged outputs
        │   │       ├── r2/
        │   │       │   ├── system.xml               # Copied OpenMM system
        │   │       │   ├── solvated.pdb             # Copied topology coordinates
        │   │       │   ├── npt.state.xml            # Starting state copied from equilibration
        │   │       │   ├── prod.chk                 # Checkpoint used to resume production
        │   │       │   ├── chunks/
        │   │       │   │   ├── 0001.dcd             # First committed trajectory chunk
        │   │       │   │   └── 0001.log             # Matching committed thermodynamic log
        │   │       │   ├── prod_keep_ligand.dcd     # Autoimaged trajectory without bulk solvent
        │   │       │   ├── stripped.prmtop          # Matching topology; use with the trajectory above
        │   │       │   ├── prod_ca.dcd              # Cα trajectory
        │   │       │   ├── ca.prmtop                # Matching Cα topology
        │   │       │   ├── prod_full.log            # Merged production log
        │   │       │   ├── md.log                   # Copy of the merged production log
        │   │       │   └── MERGE_DONE               # Verified merged outputs
        │   │       └── r3/
        │   │           ├── system.xml               # Copied OpenMM system
        │   │           ├── solvated.pdb             # Copied topology coordinates
        │   │           ├── npt.state.xml            # Starting state copied from equilibration
        │   │           ├── prod.chk                 # Checkpoint used to resume production
        │   │           ├── chunks/
        │   │           │   ├── 0001.dcd             # First committed trajectory chunk
        │   │           │   └── 0001.log             # Matching committed thermodynamic log
        │   │           ├── prod_keep_ligand.dcd     # Autoimaged trajectory without bulk solvent
        │   │           ├── stripped.prmtop          # Matching topology; use with the trajectory above
        │   │           ├── prod_ca.dcd              # Cα trajectory
        │   │           ├── ca.prmtop                # Matching Cα topology
        │   │           ├── prod_full.log            # Merged production log
        │   │           ├── md.log                   # Copy of the merged production log
        │   │           └── MERGE_DONE               # Verified merged outputs
        │   └── G12D/
        │       ├── prep/
        │       │   ├── model.pdb                    # Prepared solute coordinates
        │       │   ├── solvated.pdb                 # Protein, water and ions
        │       │   ├── system.prmtop                # Amber topology passed to equilibration
        │       │   ├── system.inpcrd                # Amber coordinates passed to equilibration
        │       │   └── PREP_DONE                    # Verified preparation outputs
        │       ├── eq/
        │       │   ├── solvated.pdb                 # Topology coordinates used by OpenMM
        │       │   ├── system.xml                   # Serialized OpenMM system
        │       │   ├── em.state.xml                 # Minimized state
        │       │   ├── nvt.state.xml                # State after NVT equilibration
        │       │   ├── npt.state.xml                # Equilibrated state passed to production
        │       │   ├── eq.pdb                       # Final equilibrated coordinates
        │       │   └── stages/
        │       │       ├── state_01.xml             # First staged restart
        │       │       └── state_24.xml             # Last of 24 staged restart states
        │       └── md/
        │           ├── r1/
        │           │   ├── system.xml               # Copied OpenMM system
        │           │   ├── solvated.pdb             # Copied topology coordinates
        │           │   ├── npt.state.xml            # Starting state copied from equilibration
        │           │   ├── prod.chk                 # Checkpoint used to resume production
        │           │   ├── chunks/
        │           │   │   ├── 0001.dcd             # First committed trajectory chunk
        │           │   │   └── 0001.log             # Matching committed thermodynamic log
        │           │   ├── prod_keep_ligand.dcd     # Autoimaged trajectory without bulk solvent
        │           │   ├── stripped.prmtop          # Matching topology; use with the trajectory above
        │           │   ├── prod_ca.dcd              # Cα trajectory
        │           │   ├── ca.prmtop                # Matching Cα topology
        │           │   ├── prod_full.log            # Merged production log
        │           │   ├── md.log                   # Copy of the merged production log
        │           │   └── MERGE_DONE               # Verified merged outputs
        │           ├── r2/
        │           │   ├── system.xml               # Copied OpenMM system
        │           │   ├── solvated.pdb             # Copied topology coordinates
        │           │   ├── npt.state.xml            # Starting state copied from equilibration
        │           │   ├── prod.chk                 # Checkpoint used to resume production
        │           │   ├── chunks/
        │           │   │   ├── 0001.dcd             # First committed trajectory chunk
        │           │   │   └── 0001.log             # Matching committed thermodynamic log
        │           │   ├── prod_keep_ligand.dcd     # Autoimaged trajectory without bulk solvent
        │           │   ├── stripped.prmtop          # Matching topology; use with the trajectory above
        │           │   ├── prod_ca.dcd              # Cα trajectory
        │           │   ├── ca.prmtop                # Matching Cα topology
        │           │   ├── prod_full.log            # Merged production log
        │           │   ├── md.log                   # Copy of the merged production log
        │           │   └── MERGE_DONE               # Verified merged outputs
        │           └── r3/
        │               ├── system.xml               # Copied OpenMM system
        │               ├── solvated.pdb             # Copied topology coordinates
        │               ├── npt.state.xml            # Starting state copied from equilibration
        │               ├── prod.chk                 # Checkpoint used to resume production
        │               ├── chunks/
        │               │   ├── 0001.dcd             # First committed trajectory chunk
        │               │   └── 0001.log             # Matching committed thermodynamic log
        │               ├── prod_keep_ligand.dcd     # Autoimaged trajectory without bulk solvent
        │               ├── stripped.prmtop          # Matching topology; use with the trajectory above
        │               ├── prod_ca.dcd              # Cα trajectory
        │               ├── ca.prmtop                # Matching Cα topology
        │               ├── prod_full.log            # Merged production log
        │               ├── md.log                   # Copy of the merged production log
        │               └── MERGE_DONE               # Verified merged outputs
        └── analysis/
            ├── WT/
            │   ├── eq/
            │   │   ├── qc.txt                       # Equilibration diagnostics
            │   │   └── qc.png                       # Thermodynamic and structural QC figure
            │   ├── r1/
            │   │   ├── rmsd.csv                     # RMSD values
            │   │   ├── rg.csv                       # Radius of gyration values
            │   │   ├── rmsf.csv                     # Per-residue fluctuations
            │   │   ├── figures/
            │   │   │   ├── rmsd_vs_time.png         # RMSD plot
            │   │   │   ├── rg_vs_time.png           # Radius of gyration plot
            │   │   │   └── rmsf_per_residue.png     # Fluctuation plot
            │   │   └── ANALYSIS_DONE                # Verified replica tables and figures
            │   ├── r2/
            │   │   ├── rmsd.csv                     # RMSD values
            │   │   ├── rg.csv                       # Radius of gyration values
            │   │   ├── rmsf.csv                     # Per-residue fluctuations
            │   │   ├── figures/
            │   │   │   ├── rmsd_vs_time.png         # RMSD plot
            │   │   │   ├── rg_vs_time.png           # Radius of gyration plot
            │   │   │   └── rmsf_per_residue.png     # Fluctuation plot
            │   │   └── ANALYSIS_DONE                # Verified replica tables and figures
            │   ├── r3/
            │   │   ├── rmsd.csv                     # RMSD values
            │   │   ├── rg.csv                       # Radius of gyration values
            │   │   ├── rmsf.csv                     # Per-residue fluctuations
            │   │   ├── figures/
            │   │   │   ├── rmsd_vs_time.png         # RMSD plot
            │   │   │   ├── rg_vs_time.png           # Radius of gyration plot
            │   │   │   └── rmsf_per_residue.png     # Fluctuation plot
            │   │   └── ANALYSIS_DONE                # Verified replica tables and figures
            │   └── aggregate/
            │       ├── rmsd_avg.csv                 # Replica mean and standard deviation
            │       ├── rg_avg.csv                   # Replica mean and standard deviation
            │       ├── rmsf_avg.csv                 # Replica mean and standard deviation
            │       ├── figures/
            │       │   ├── rmsd_aggregate.png       # Replica RMSD summary
            │       │   ├── rg_aggregate.png         # Replica compactness summary
            │       │   └── rmsf_aggregate.png       # Replica fluctuation summary
            │       └── AGGREGATE_DONE               # Verified multi-replica summary
            ├── G12D/
            │   ├── eq/
            │   │   ├── qc.txt                       # Equilibration diagnostics
            │   │   └── qc.png                       # Thermodynamic and structural QC figure
            │   ├── r1/
            │   │   ├── rmsd.csv                     # RMSD values
            │   │   ├── rg.csv                       # Radius of gyration values
            │   │   ├── rmsf.csv                     # Per-residue fluctuations
            │   │   ├── figures/
            │   │   │   ├── rmsd_vs_time.png         # RMSD plot
            │   │   │   ├── rg_vs_time.png           # Radius of gyration plot
            │   │   │   └── rmsf_per_residue.png     # Fluctuation plot
            │   │   └── ANALYSIS_DONE                # Verified replica tables and figures
            │   ├── r2/
            │   │   ├── rmsd.csv                     # RMSD values
            │   │   ├── rg.csv                       # Radius of gyration values
            │   │   ├── rmsf.csv                     # Per-residue fluctuations
            │   │   ├── figures/
            │   │   │   ├── rmsd_vs_time.png         # RMSD plot
            │   │   │   ├── rg_vs_time.png           # Radius of gyration plot
            │   │   │   └── rmsf_per_residue.png     # Fluctuation plot
            │   │   └── ANALYSIS_DONE                # Verified replica tables and figures
            │   ├── r3/
            │   │   ├── rmsd.csv                     # RMSD values
            │   │   ├── rg.csv                       # Radius of gyration values
            │   │   ├── rmsf.csv                     # Per-residue fluctuations
            │   │   ├── figures/
            │   │   │   ├── rmsd_vs_time.png         # RMSD plot
            │   │   │   ├── rg_vs_time.png           # Radius of gyration plot
            │   │   │   └── rmsf_per_residue.png     # Fluctuation plot
            │   │   └── ANALYSIS_DONE                # Verified replica tables and figures
            │   └── aggregate/
            │       ├── rmsd_avg.csv                 # Replica mean and standard deviation
            │       ├── rg_avg.csv                   # Replica mean and standard deviation
            │       ├── rmsf_avg.csv                 # Replica mean and standard deviation
            │       ├── figures/
            │       │   ├── rmsd_aggregate.png       # Replica RMSD summary
            │       │   ├── rg_aggregate.png         # Replica compactness summary
            │       │   └── rmsf_aggregate.png       # Replica fluctuation summary
            │       └── AGGREGATE_DONE               # Verified multi-replica summary
            └── compare/
                ├── rmsd.png                         # Variant RMSD comparison
                ├── rg.png                           # Variant compactness comparison
                ├── rmsf.png                         # Variant fluctuation comparison
                ├── figures/
                │   ├── rmsd_avg.png                 # Copy of RMSD comparison
                │   ├── rg_avg.png                   # Copy of compactness comparison
                │   └── rmsf_avg.png                 # Copy of fluctuation comparison
                └── summary.txt                      # Comparison statistics
```

</details>
