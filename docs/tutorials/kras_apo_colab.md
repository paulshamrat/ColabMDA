# KRAS (4LDJ) Apo State: Google Colab Walkthrough

This tutorial provides a streamlined simulation walkthrough for the **Apo state** of the **KRAS GTPase** (UniProt **P01116**, PDB **4LDJ** Chain A, residues **1–169**), comparing the Wild-Type (WT) protein against the oncogenic point mutant **G12D**.

Everything is executed on a **Google Colab cloud GPU via the `colab` CLI**, with all outputs, checkpoints, and figures continuously persisted inside your **Google Drive**.

---

## Directory Structure

All files reside inside your mounted Google Drive under a standardized hierarchy:

```text
/content/drive/MyDrive/research/projects/kras-example/
└── data/
    └── apo/
        ├── structure/
        │   └── P01116_4LDJ/                      # Canonical UniProt P01116 + Template PDB 4LDJ
        │       ├── raw/                          # Untouched experimental inputs (cmda fetch)
        │       │   ├── 4ldj.cif                  # Raw crystallographic archive record
        │       │   ├── template.pdb              # Cleaned monomer (waters/buffers stripped)
        │       │   └── uniprot.fasta             # Canonical sequence for biological residue numbering
        │       │
        │       ├── variants/                     # Modeled coordinates (cmda model)
        │       │   ├── WT/
        │       │   │   ├── WT.pdb                # Canonical modeled Wild-Type coordinate
        │       │   │   ├── target.B99990001.pdb  # Modeller build
        │       │   │   ├── target.B99990001_with_cryst.pdb # Model with crystal headers
        │       │   │   ├── target.ini            # Starting coordinate mesh
        │       │   │   ├── target.rsr            # Spatial restraints
        │       │   │   ├── target.sch            # Optimization schedule
        │       │   │   ├── target.D00000001      # Restraint progress log
        │       │   │   ├── target.V99990001      # Restraint violation profile
        │       │   │   ├── target.ali            # PIR target sequence alignment
        │       │   │   ├── alignment.ali         # PIR template pairwise alignment
        │       │   │   ├── P01116.fasta          # Sequence reference
        │       │   │   ├── 4ldj_clean.pdb        # Sanitized template
        │       │   │   ├── 4ldj_orig.pdb         # Raw template backup
        │       │   │   ├── 4ldj.pdb              # Isolated chain coordinate
        │       │   │   ├── pipeline.log          # Modeller build log
        │       │   │   └── MODEL_DONE            # Modeling completion gate
        │       │   │
        │       │   └── G12D/                     # Oncogenic point mutant
        │       │       ├── G12D.pdb              # Mutated coordinate
        │       │       ├── WT_G12D.pdb           # Direct mutation coordinate
        │       │       ├── WT_G12D_with_cryst.pdb # Mutant with crystal headers
        │       │       ├── mutate.log            # Rotamer sampling & energy log
        │       │       └── MODEL_DONE            # Mutant completion gate
        │       │
        │       └── STRUCTURE_DONE                # Verification QC gate marker
        │
        ├── simulation/
        │   ├── WT/
        │   │   ├── prep/              # Solvated Amber system
        │   │   │   ├── model.pdb         # Input receptor coordinates
        │   │   │   ├── solvated.pdb      # Solvated visualization PDB (protein + OPC water box + ions)
        │   │   │   ├── system.prmtop     # Amber ff19SB parameter/topology file
        │   │   │   ├── system.inpcrd     # Initial orthogonal coordinates
        │   │   │   └── system.xml        # Serialized OpenMM system
        │   │   │
        │   │   ├── eq/                # 24-stage equilibration
        │   │   │   ├── eq.pdb            # Final relaxed coordinates
        │   │   │   ├── system.xml        # OpenMM system topology & parameters
        │   │   │   ├── state.xml         # OpenMM thermodynamic state snapshot
        │   │   │   └── stages/           # Granular restart checkpoints (state_01.xml ... state_24.xml)
        │   │   │
        │   │   └── md/                # Production replicas (100 ns)
        │   │       ├── r1/               # Production replica 1
        │   │       │   ├── md.dcd        # Continuous merged trajectory (full system)
        │   │       │   ├── md.log        # Continuous merged thermodynamic log
        │   │       │   ├── state.xml     # Latest OpenMM simulation state snapshot
        │   │       │   ├── prod.chk      # OpenMM checkpoint binary for resume
        │   │       │   ├── stripped.prmtop      # Solvent-stripped topology for PyMOL & analysis
        │   │       │   ├── prod_keep_ligand.dcd # Continuous PBC-autoimaged stripped trajectory
        │   │       │   ├── ca.prmtop     # Cα backbone topology
        │   │       │   ├── prod_ca.dcd   # Cα backbone trajectory
        │   │       │   ├── MERGE_DONE    # Merge completion marker
        │   │       │   └── chunks/       # Raw crash-recovery chunks (0001.dcd, 0001.log, ...)
        │   │       ├── r2/               # Independent velocity seed (triplicate)
        │   │       └── r3/               # Independent velocity seed (triplicate)
        │   │
        │   └── G12D/
        │       ├── prep/
        │       ├── eq/
        │       └── md/{r1, r2, r3}/
        │
        └── analysis/                     # Mirrored deliverables (Zero JSON policy)
            ├── WT/
            │   ├── eq/                # Equilibration QC deliverables (qc.png, qc.txt)
            │   └── md/
            │       ├── r1/               # Per-replica metrics (rmsd.png, rg.png, rmsf.png, summary.txt)
            │       ├── r2/
            │       └── r3/
            ├── G12D/
            │   ├── eq/                # qc.png, qc.txt
            │   └── md/{r1, r2, r3}/
            └── compare/                  # Flat cross-variant comparison (rmsd.png, rg.png, rmsf.png, summary.txt)
```

---

## 1. Connect to Google Colab

### 1.1. Local Setup
Run this one-time setup in your local terminal on **Linux (Ubuntu)** or Windows via **[WSL2](https://learn.microsoft.com/en-us/windows/wsl/install)**:

```bash
# 1. Install cmda, OpenMM, AmberTools, PyMOL, and google-colab-cli:
bash -c "$(curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/dev1/scripts/setup.sh)"

# 2. Activate environment:
source ~/miniforge3/etc/profile.d/conda.sh
conda activate cmda

# 3. Log in to Google Colab (one-time browser authorization):
colab sessions
```

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
curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/dev1/scripts/setup.sh | bash

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

# 2. Model canonical Wild-Type coordinates (closes missing crystallographic loops):
cmda model --state apo --variant WT

# 3. Build Point Mutant (G12D) on the verified WT backbone:
cmda model --state apo --variant G12D
```

---

## 4. Structure Preparation & Solvation

ColabMDA runs **PROPKA 3** to calculate 3D electrostatic microenvironments, assign explicit AMBER titratable protonation states, solvate with OPC 4-point water, and neutralize net charge:

```bash
cmda prep --state apo --variant WT
cmda prep --state apo --variant G12D
```

| Residue | State | Physical Criteria |
| :--- | :--- | :--- |
| **Histidine** | `HID` / `HIE` / `HIP` | Neutral $\text{N}_\delta$, neutral $\text{N}_\epsilon$ (default), or positively charged doubly protonated ($pK_a > 7.4$). |
| **Acidic** | `ASP` / `GLU` | Standard carboxylate anions (`ASH`/`GLH` assigned if buried with $pK_a > 7.4$). |
| **Basic** | `LYS` / `ARG` | Positively charged physiological cations (`LYN` if neutral). |
| **Cysteine** | `CYS` / `CYX` | Free neutral thiol (`CYS`) or disulfide bridge (`CYX`). |

*Outputs in `data/apo/simulation/<variant>/prep/`:* `solvated.pdb`, `system.prmtop`, `system.inpcrd`, and `system.xml`.

---

## 5. Equilibration & QC

Initialize simulation directories and run the progressive 24-stage VarMDyn relaxation ladder. This heats the system from 50 K to 300 K and gradually releases positional restraints:

```bash
# 1. Initialize simulation directories:
cmda stage --state apo --variant WT
cmda stage --state apo --variant G12D

# 2. Execute 24-stage progressive equilibration (EM -> NVT -> NPT -> QC):
cmda sim eq --state apo --variant WT
cmda sim eq --state apo --variant G12D
```

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

---

## 6. Production MD (100 ns)

### Why Triplicates?
Molecular dynamics is inherently stochastic. A single trajectory can easily become trapped in a metastable conformational sub-state. Running **triplicates** ($N=3$ independent trajectories: `r1`, `r2`, `r3`) ensures conformational sampling convergence, allowing you to compute true **Mean ± SD** statistical error bands for RMSD, Radius of Gyration, and residue fluctuations.

> 💡 **Why Shared `eq/`?**
> Equilibration only needs to be run **once** per system. Replicas (`r1`, `r2`, `r3`) fork directly from the final relaxed coordinates in `eq/state.xml` using distinct Maxwell-Boltzmann random velocity seeds. This eliminates redundant equilibration, saving hours of GPU time on Google Colab while guaranteeing strict statistical independence.

### 6.1. Launch Production Runs (`r1`)
Launch production dynamics for the primary replica (`r1`) of Wild-Type and G12D:

```bash
# 1. Run Wild-Type replica 1 (100 ns):
cmda sim md --state apo --variant WT --replica r1 --total-ns 100

# 2. Run G12D mutant replica 1 (100 ns):
cmda sim md --state apo --variant G12D --replica r1 --total-ns 100
```

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
>    curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/dev1/scripts/setup.sh | bash
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

### 6.3. Generated Output Files & Downstream Handover

During simulation, ColabMDA populates `data/apo/simulation/<variant>/md/<replica>/`:

```text
simulation/WT/md/r1/
├── chunk_001.dcd ... chunk_100.dcd   # 1-ns raw trajectory chunks (100 frames each, 10,000 frames total)
├── chunk_001.log ... chunk_100.log   # 1-ns thermodynamic logs (temperature, volume, density, energy)
├── state.xml                         # Latest thermodynamic state snapshot
└── prod.chk                          # Checkpoint binary for auto-resuming after timeouts
```

> **Downstream Handover:** These 100 raw 1-ns chunk files are handed over to **Section 7 (`cmda merge`)**, which concatenates them into a continuous unified trajectory (`md.dcd`), eliminates periodic boundary condition (PBC) box jumps, and strips solvent to produce `prod_clean.dcd` for statistical analysis.

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
# 1. Merge chunks & autoimage PBC (produces continuous md.dcd and stripped trajectories):
cmda merge --state apo --variant WT
cmda merge --state apo --variant G12D

# 2. Audit frame counts & chunk continuity:
cmda status --state apo --variant WT --replica r1
cmda status --state apo --variant G12D --replica r1

# 3. Calculate single-variant biophysical metrics (RMSD, Rg, RMSF):
cmda analyze --state apo --variant WT
cmda analyze --state apo --variant G12D

# 4. Generate direct comparative overlay figures (WT vs G12D with Mean ± SD bands):
cmda compare --state apo --variants WT,G12D
```

*Key Deliverables Generated in `data/apo/analysis/compare/`:*
* `summary.txt`: Executive statistical table with Mean ± SD for RMSD, Radius of Gyration, and Temperature.
* `rmsd.png`: Comparative backbone RMSD time-series.
* `rg.png`: Comparative compactness profiles.
* `rmsf.png`: Per-residue fluctuation profiles highlighting Switch I (residues 30–38) and Switch II (residues 59–67) loop flexibility.

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

> [!TIP]
> **What `cmda pull` Does:**
> 1. Packages analysis CSVs, logs, summary manifests, and PNG plots from Google Drive.
> 2. Excludes heavy raw simulation checkpoints (`chunks/*.dcd`, `.chk`) by default so the download completes quickly. Passing `--traj` pulls lightweight solvent-stripped trajectories (`prod_keep_ligand.dcd`), while `--all` downloads full raw chunks.
> 3. Extracts files non-destructively into your local `data/` folder.

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
