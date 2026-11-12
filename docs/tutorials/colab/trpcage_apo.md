# Trp-Cage (1L2Y) Apo State: Google Colab Walkthrough

This tutorial provides a lightweight, 5-minute simulation walkthrough for the **Apo state** of the 20-residue **Trp-cage miniprotein** (PDB **1L2Y**, NMR Model 1, residues **1–20**). It serves as a rapid benchmark to verify that your environment, Google Drive synchronization, Colab GPU runtime, and simulation checkpointing are fully operational.

Everything is executed on a **Google Colab cloud GPU via the `colab` CLI**, with all outputs, checkpoints, and figures continuously persisted inside your **Google Drive**.

---

## Directory Structure

All files reside inside your mounted Google Drive under a standardized hierarchy:

```text
/content/drive/MyDrive/research/projects/trpcage-example/
└── data/
    └── apo/
        ├── structure/
        │   └── 1L2Y/                             # Target: NMR PDB 1L2Y
        │       ├── raw/                          # Untouched experimental inputs (cmda fetch)
        │       │   ├── 1l2y.cif                  # Raw PDB archive record
        │       │   ├── template.pdb              # Cleaned monomer (waters/buffers stripped)
        │       │   └── uniprot.fasta             # Target sequence for biological numbering
        │       │
        │       ├── variants/                     # Modeled coordinates (cmda model)
        │       │   └── WT/
        │       │       ├── WT.pdb                # Canonical modeled Wild-Type coordinate
        │       │       ├── target.B99990001.pdb  # Modeller build
        │       │       ├── target.B99990001_with_cryst.pdb # Model with crystal headers
        │       │       ├── target.ini            # Starting coordinate mesh
        │       │       ├── target.rsr            # Spatial restraints
        │       │       ├── target.sch            # Optimization schedule
        │       │       ├── target.D00000001      # Restraint progress log
        │       │       ├── target.V99990001      # Restraint violation profile
        │       │       ├── target.ali            # PIR target sequence alignment
        │       │       ├── alignment.ali         # PIR template pairwise alignment
        │       │       ├── 1l2y_clean.pdb        # Sanitized template
        │       │       ├── 1l2y_orig.pdb         # Raw template backup
        │       │       ├── 1l2y.pdb              # Isolated coordinate
        │       │       ├── pipeline.log          # Modeller build log
        │       │       └── MODEL_DONE            # Modeling completion gate
        │       │
        │       └── STRUCTURE_DONE                # Verification QC gate marker
        │
        ├── simulation/
        │   └── WT/
        │       ├── prep/                         # Solvated Amber system
        │       │   ├── model.pdb                 # Input receptor coordinates
        │       │   ├── solvated.pdb              # Solvated visualization PDB (protein + OPC water box + ions)
        │       │   ├── system.prmtop             # Amber ff19SB parameter/topology file
        │       │   ├── system.inpcrd             # Initial orthogonal coordinates
        │       │   └── system.xml                # Serialized OpenMM system
        │       │
        │       ├── eq/                           # Equilibration directory
        │       │   ├── eq.pdb                    # Final relaxed coordinates
        │       │   ├── system.xml                # OpenMM system topology & parameters
        │       │   ├── state.xml                 # OpenMM thermodynamic state snapshot
        │       │   └── stages/                   # Granular restart checkpoints
        │       │
        │       └── md/                           # Production replicas (1 ns)
        │           └── r1/                       # Production replica 1
        │               ├── md.dcd                # Continuous merged trajectory (full system)
        │               ├── md.log                # Continuous merged thermodynamic log
        │               ├── state.xml             # Latest OpenMM simulation state snapshot
        │               ├── prod.chk              # OpenMM checkpoint binary for resume
        │               ├── stripped.prmtop       # Solvent-stripped topology for PyMOL & analysis
        │               ├── prod_keep_ligand.dcd  # Continuous PBC-autoimaged stripped trajectory
        │               ├── ca.prmtop             # Cα backbone topology
        │               ├── prod_ca.dcd           # Cα backbone trajectory
        │               ├── MERGE_DONE            # Merge completion marker
        │               └── chunks/               # Raw crash-recovery chunks
        │
        └── analysis/                             # Mirrored deliverables (Zero JSON policy)
            └── WT/
                ├── eq/                           # Equilibration QC deliverables (qc.png, qc.txt)
                └── r1/                           # Per-replica metrics (figures/, rmsd.csv, rg.csv, rmsf.csv)
                    ├── figures/
                    │   ├── rmsd_vs_time.png
                    │   ├── rg_vs_time.png
                    │   └── rmsf_per_residue.png
                    ├── rmsd.csv
                    ├── rg.csv
                    ├── rmsf.csv
                    └── ANALYSIS_DONE
```

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
mkdir -p /content/drive/MyDrive/research/projects/trpcage-example
cd /content/drive/MyDrive/research/projects/trpcage-example
```

> [!NOTE]
> Academic Modeller keys are free at [salilab.org/modeller/registration.html](https://salilab.org/modeller/registration.html).

---

## 3. Modeling & Structure Verification

Fetch the 20-residue NMR coordinates (`1L2Y`, Model 1, Chain A) directly from the PDB archive and generate the canonical modeled Wild-Type structure:

```bash
# 1. Download 1L2Y (Chain A, residues 1–20) & verify coordinates:
cmda fetch --state apo --pdb-id 1L2Y --chain A

# 2. Model canonical Wild-Type coordinates:
cmda model --state apo --variant WT
```

*What `cmda fetch` generates in `data/apo/structure/1L2Y/raw/`:*

| File | Purpose | Scientific Role |
| :--- | :--- | :--- |
| **`1l2y.cif`** | Archive Copy | Raw crystallographic/NMR record directly from the RCSB PDB archive. |
| **`template.pdb`** | Modeller Input | Stripped 20-residue monomer coordinates (`NLYIQWLKDGGPSSGRPPPS`). |
| **`uniprot.fasta`** | Sequence Target | Canonical FASTA sequence for residue numbering. |
| **`STRUCTURE_DONE`** | Gate Marker | Verification marker confirming all 20 Cα atoms and coordinates passed QC gates. |

*What `cmda model` outputs in `data/apo/structure/1L2Y/variants/WT/`:*
* `WT.pdb`: Canonical modeled Wild-Type coordinates.
* `MODEL_DONE`: Modeling completion gate marker.

> **Downstream Handover:** These validated coordinates (`WT.pdb`) are handed over to **Section 4 (`cmda prep`)** for PROPKA protonation and Amber OPC water box solvation.

---

## 4. Structure Preparation & Solvation

ColabMDA runs **PROPKA 3** at physiological pH 7.4 to evaluate 3D electrostatic microenvironments, assign explicit AMBER titratable protonation states, solvate with OPC 4-point water, and neutralize net charge:

```bash
cmda prep --state apo --variant WT
```

| Residue | State | Physical Criteria |
| :--- | :--- | :--- |
| **Histidine** | `HID` / `HIE` / `HIP` | Neutral $\text{N}_\delta$, neutral $\text{N}_\epsilon$, or doubly protonated positive ($pK_a > 7.4$). |
| **Acidic** | `ASP` / `GLU` | Standard carboxylate anions (`ASH`/`GLH` assigned if buried with $pK_a > 7.4$). |
| **Basic** | `LYS` / `ARG` | Fully protonated cations (Lys8 and Arg16 in Trp-cage). |
| **Cysteine** | `CYS` / `CYX` | Free neutral thiol (`CYS`) or disulfide bridge (`CYX`). |

*What `cmda prep` outputs in `data/apo/simulation/WT/prep/`:*

```text
prep/
├── model.pdb        # Input receptor coordinates copied for reproducibility
├── solvated.pdb     # Solvated coordinates (protein, ~1,800 OPC waters, neutralizing Cl⁻ ion)
├── system.prmtop    # Amber ff19SB topology containing all force field parameters
├── system.inpcrd    # Solvated, charge-neutralized initial coordinates
└── system.xml       # Serialized OpenMM system (physical forces, constraints, box vectors)
```

> **Downstream Handover:** These solvated parameters are loaded by OpenMM for GPU equilibration in **Section 5 (`cmda sim eq`)**.

---

## 5. Equilibration & QC

Initialize simulation directories and run equilibration. For this 5-minute quickstart, use `--protocol quick` (runs restrained energy minimization, 100 ps NVT heating, and 100 ps NPT density relaxation in ~30 seconds):

```bash
# 1. Initialize simulation directories:
cmda stage --state apo --variant WT

# 2. Execute quick equilibration (~30 seconds on GPU):
cmda sim eq --state apo --variant WT --protocol quick
```

> [!TIP]
> **Full 24-Stage Progressive Ladder:** Omitting `--protocol quick` executes the complete 24-stage heating and restraint release ladder (`--protocol staged`).

### QC Verification
Verify thermodynamic and structural convergence:

```bash
cat data/apo/analysis/WT/eq/qc.txt
```
*Expected Output:* `Status: PASSED`, `Cα RMSD: ~0.8 Å` (Confirms temperature $\approx 300\text{ K}$, density $\approx 1.01\text{ g/mL}$, fold intact).

*What `cmda sim eq` outputs:*
* In `data/apo/simulation/WT/eq/`: `eq.pdb`, `system.xml`, `state.xml`, `stages/`.
* In `data/apo/analysis/WT/eq/`: `qc.png` (6-panel publication QC plot) and `qc.txt` (pass/fail diagnostic receipt).

> **Downstream Handover:** The final thermodynamic state (`eq/state.xml`) and topology serve as the physical starting point for independent production MD in **Section 6 (`cmda sim md`)**.

---

## 6. Checkpointed Production MD (1 ns)

Run 1 ns of checkpointed production molecular dynamics on the Colab GPU (~1 minute on T4):

```bash
cmda sim md --state apo --variant WT --replica r1 --total-ns 1 --checkpoint-ps 100 --traj-interval 10
```

> [!IMPORTANT]
> **💡 How to Resume After a Session Timeout:**
> If your Colab session expires, checkpoints remain safe on Google Drive:
> 1. Launch a new session & mount Drive from your local terminal:
>    ```bash
>    colab new -s cmda --gpu T4
>    colab drivemount --session cmda
>    colab console -s cmda
>    ```
> 2. Re-activate environment & navigate to Drive in Colab console (~90 seconds):
>    ```bash
>    [ ! -d "/root/miniforge3/envs/cmda" ] && curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/main/scripts/setup.sh | bash
>    source /root/miniforge3/etc/profile.d/conda.sh
>    conda activate cmda
>    cd /content/drive/MyDrive/research/projects/trpcage-example
>    ```
> 3. Re-run the exact same `cmda sim md` command:
>    ```bash
>    cmda sim md --state apo --variant WT --replica r1 --total-ns 1 --checkpoint-ps 100 --traj-interval 10
>    ```
>    *ColabMDA automatically detects `prod.chk` on Google Drive and seamlessly resumes dynamics.*

### Simulation Parameters & Frame Saving Cadence

| Parameter | Value | Biophysical Rationale & Details |
| :--- | :--- | :--- |
| **Timestep ($\Delta t$)** | **$2\text{ fs}$ ($0.002\text{ ps}$)** | Standard Langevin middle integrator with SHAKE/SETTLE constraints on all bonds involving hydrogen. |
| **Frame Saving Interval** | **Every $10\text{ ps}$ ($5,000\text{ steps}$)** | Writes coordinates to `.dcd` trajectory every 10 ps (`--traj-interval 10`). |
| **Thermodynamic Log Interval** | **Every $10\text{ ps}$ ($5,000\text{ steps}$)** | Records instantaneous temperature, potential energy, kinetic energy, volume, and density to `.log`. |
| **Checkpoint Cadence** | **Every $100\text{ ps}$** | Serializes binary state checkpoint `prod.chk` to Google Drive every 100 ps. |
| **Total Frames Generated** | **$100\text{ frames}$** | $\frac{1,000\text{ ps}}{10\text{ ps}} = 100\text{ frames}$ for a fast $1\text{ ns}$ benchmark run. |
| **Storage Footprint** | **~35 MB $\to$ <1 MB** | Stripped trajectories (`cmda merge`) reduce trajectory size to <1 MB. |

---

## 7. Merge & Analysis

Concatenate raw trajectory chunks, apply periodic boundary condition (PBC) autoimaging, strip bulk water, audit frame continuity, and calculate biophysical metrics:

```bash
# 1. Merge chunks & autoimage PBC (produces continuous md.dcd and stripped trajectories):
cmda sim merge --state apo --variant WT --replica r1

# 2. Audit frame counts & chunk continuity:
cmda status --state apo --variant WT --replica r1

# 3. Calculate biophysical metrics (RMSD, Rg, RMSF):
cmda analyze --state apo --variant WT --replica r1
```

*What `cmda sim merge` outputs in `data/apo/simulation/WT/md/r1/`:*
* `md.dcd`: Continuous concatenated solvated trajectory.
* `md.log`: Continuous concatenated thermodynamic log table.
* `prod_keep_ligand.dcd`: PBC-autoimaged, solvent-stripped trajectory (<1 MB).
* `stripped.prmtop`: Amber topology matching stripped trajectory (304 atoms).
* `ca.prmtop`: Cα-only topology (20 atoms).
* `prod_ca.dcd`: Cα-only trajectory.
* `MERGE_DONE`: Merge verification gate marker.

*What `cmda analyze` outputs in `data/apo/analysis/WT/r1/`:*
* `figures/rmsd_vs_time.png`: Backbone Cα RMSD over 1 ns.
* `figures/rg_vs_time.png`: Radius of Gyration compactness.
* `figures/rmsf_per_residue.png`: Per-residue root mean square fluctuations for residues 1–20.
* `rmsd.csv`, `rg.csv`, `rmsf.csv`: Raw numeric data files.
* `ANALYSIS_DONE`: Analysis completion marker.

---

## 8. Download Deliverables & PyMOL Visualization

### 8.1. Download Deliverables to Local Machine
In your **local terminal**, download the analysis deliverables and stripped trajectory using the `colab` CLI:

```bash
mkdir -p ~/research/projects/trpcage-example/data/apo

# Download analysis figures and summaries:
colab download -s cmda /content/drive/MyDrive/research/projects/trpcage-example/data/apo/analysis ~/research/projects/trpcage-example/data/apo/analysis

# Download solvent-stripped trajectory and coordinates:
colab download -s cmda /content/drive/MyDrive/research/projects/trpcage-example/data/apo/simulation/WT/md/r1/prod_keep_ligand.dcd ~/research/projects/trpcage-example/data/apo/
colab download -s cmda /content/drive/MyDrive/research/projects/trpcage-example/data/apo/structure/1L2Y/variants/WT/WT.pdb ~/research/projects/trpcage-example/data/apo/
```

### 8.2. PyMOL Trajectory Visualization
Inspect the production trajectory in PyMOL locally:

```bash
cd ~/research/projects/trpcage-example/data/apo
pymol WT.pdb prod_keep_ligand.dcd
```

In the PyMOL console:
```text
intra_fit (all and name CA)
show cartoon
color cyan, all and name C*
show sticks, resi 6 or resi 12 or resi 17-19
color yellow, resi 6 and name C*
color orange, (resi 12 or resi 17-19) and name C*
mplay
```

*Confirm the central **Trp6** indole ring remains securely packed inside the hydrophobic cage formed by **Pro12**, **Pro17**, **Pro18**, and **Pro19** throughout the simulation.*

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
