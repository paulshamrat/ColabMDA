# Trp-Cage (1L2Y) Apo State: Google Colab Walkthrough

This tutorial provides a lightweight, 5-minute simulation walkthrough for the **Apo state** of the 20-residue **Trp-cage miniprotein** (PDB **1L2Y**, NMR Model 1, residues **1–20**). It serves as a rapid benchmark to verify that your environment, Google Drive synchronization, Colab GPU runtime, and simulation checkpointing are fully operational.

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
```

<details>
<summary><b>View stage directory tree</b></summary>

```text
data/apo/structure/1L2Y/
├── raw/
│   ├── 1l2y.cif                                     # Experimental archive; PDB fallback is possible
│   ├── template.pdb                                 # Cleaned template coordinates
│   └── uniprot.fasta                                # Sequence used for modeling
└── STRUCTURE_DONE                                   # Completed fetch and mapping
```

</details>

```bash
# 2. Model canonical Wild-Type coordinates:
cmda model --state apo --variant WT
```

<details>
<summary><b>View stage directory tree</b></summary>

```text
data/apo/structure/1L2Y/variants/WT/
├── WT.pdb                                           # Coordinates passed to preparation
└── MODEL_DONE                                       # Completed model verification
```

</details>

> **Downstream Handover:** These validated coordinates (`WT.pdb`) are handed over to **Section 4 (`cmda prep`)** for PROPKA protonation and Amber OPC water box solvation.

---

## 4. Structure Preparation & Solvation

ColabMDA runs **PROPKA 3** at physiological pH 7.4 to evaluate 3D electrostatic microenvironments, assign explicit AMBER titratable protonation states, solvate with OPC 4-point water, and neutralize net charge:

```bash
cmda prep --state apo --variant WT
```

<details>
<summary><b>View stage directory tree</b></summary>

```text
data/apo/simulation/WT/prep/
├── model.pdb                                        # Prepared solute coordinates
├── solvated.pdb                                     # Protein, water and ions
├── system.prmtop                                    # Amber topology passed to equilibration
├── system.inpcrd                                    # Amber coordinates passed to equilibration
└── PREP_DONE                                        # Verified preparation outputs
```

</details>

| Residue | State | Physical Criteria |
| :--- | :--- | :--- |
| **Histidine** | `HID` / `HIE` / `HIP` | Neutral $\text{N}_\delta$, neutral $\text{N}_\epsilon$, or doubly protonated positive ($pK_a > 7.4$). |
| **Acidic** | `ASP` / `GLU` | Standard carboxylate anions (`ASH`/`GLH` assigned if buried with $pK_a > 7.4$). |
| **Basic** | `LYS` / `ARG` | Fully protonated cations (Lys8 and Arg16 in Trp-cage). |
| **Cysteine** | `CYS` / `CYX` | Free neutral thiol (`CYS`) or disulfide bridge (`CYX`). |

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

<details>
<summary><b>View stage directory tree</b></summary>

```text
data/apo/
├── simulation/
│   └── WT/
│       └── eq/
│           ├── solvated.pdb                         # Topology coordinates used by OpenMM
│           ├── system.xml                           # Serialized OpenMM system
│           ├── em.state.xml                         # Minimized state
│           ├── nvt.state.xml                        # State after NVT equilibration
│           ├── npt.state.xml                        # Equilibrated state passed to production
│           └── eq.pdb                               # Final equilibrated coordinates
└── analysis/
    └── WT/
        └── eq/
            ├── qc.txt                               # Equilibration diagnostics
            └── qc.png                               # Thermodynamic and structural QC figure
```

</details>

> [!TIP]
> **Full 24-Stage Progressive Ladder:** Omitting `--protocol quick` executes the complete 24-stage heating and restraint release ladder (`--protocol staged`).

### QC Verification
Verify thermodynamic and structural convergence:

```bash
cat data/apo/analysis/WT/eq/qc.txt
```
*Expected Output:* `Status: PASSED`, `Cα RMSD: ~0.8 Å` (Confirms temperature $\approx 300\text{ K}$, density $\approx 1.01\text{ g/mL}$, fold intact).

> **Downstream Handover:** The final thermodynamic state (`eq/npt.state.xml`) and topology serve as the physical starting point for independent production MD in **Section 6 (`cmda sim md`)**.

---

## 6. Checkpointed Production MD (1 ns)

Run 1 ns of checkpointed production molecular dynamics on the Colab GPU (~1 minute on T4):

```bash
cmda sim md --state apo --variant WT --replica r1 --total-ns 1 --checkpoint-ps 100 --traj-interval 10
```

<details>
<summary><b>View stage directory tree</b></summary>

```text
data/apo/simulation/WT/md/r1/
├── system.xml                                       # Copied OpenMM system
├── solvated.pdb                                     # Copied topology coordinates
├── npt.state.xml                                    # Starting state copied from equilibration
├── prod.chk                                         # Checkpoint used to resume production
└── chunks/
    ├── 0001.dcd                                     # First committed trajectory chunk
    └── 0001.log                                     # Matching committed thermodynamic log
```

</details>

Each replica has its own folder. Committed chunks are numbered `0001`, `0002`,
and so on; `.tmp` files belong to the chunk currently being written.

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
# 1. Merge chunks & autoimage PBC (produces autoimaged, stripped trajectories):
cmda merge --state apo --variant WT --replica r1
```

<details>
<summary><b>View stage directory tree</b></summary>

```text
data/apo/simulation/WT/md/r1/
├── traj.dcd                             # Autoimaged trajectory without bulk solvent
├── top.prmtop                                  # Matching topology; use with the trajectory above
├── trajCa.dcd                                      # Cα trajectory
├── topCa.prmtop                                        # Matching Cα topology
├── md.log                                           # Copy of the merged production log
└── MERGE_DONE                                       # Verified merged outputs
```

</details>

```bash
# 2. Audit frame counts & chunk continuity:
cmda status --state apo --variant WT --replica r1
```

```bash
# 3. Calculate biophysical metrics (RMSD, Rg, RMSF):
cmda analyze --state apo --variant WT --replica r1
```

<details>
<summary><b>View stage directory tree</b></summary>

```text
data/apo/analysis/WT/r1/
├── rmsd.csv                                         # RMSD values
├── rg.csv                                           # Radius of gyration values
├── rmsf.csv                                         # Per-residue fluctuations
├── figures/
│   ├── rmsd_vs_time.png                             # RMSD plot
│   ├── rg_vs_time.png                               # Radius of gyration plot
│   └── rmsf_per_residue.png                         # Fluctuation plot
└── ANALYSIS_DONE                                    # Verified replica tables and figures
```

</details>

---

## 8. Download Deliverables & PyMOL Visualization

### 8.1. Download Deliverables to Local Machine
In your **local terminal**, download the analysis deliverables and stripped trajectory using the `colab` CLI:

```bash
mkdir -p ~/research/projects/trpcage-example/data/apo

# Download analysis figures and summaries:
colab download -s cmda /content/drive/MyDrive/research/projects/trpcage-example/data/apo/analysis ~/research/projects/trpcage-example/data/apo/analysis

# Download solvent-stripped trajectory and coordinates:
colab download -s cmda /content/drive/MyDrive/research/projects/trpcage-example/data/apo/simulation/WT/md/r1/traj.dcd ~/research/projects/trpcage-example/data/apo/
colab download -s cmda /content/drive/MyDrive/research/projects/trpcage-example/data/apo/structure/1L2Y/variants/WT/WT.pdb ~/research/projects/trpcage-example/data/apo/
```

<details>
<summary><b>View transferred directory tree</b></summary>

```text
~/research/projects/trpcage-example/data/apo/
├── analysis/               # Downloaded analysis tables and figures
├── WT.pdb                  # Model for local visualization
└── traj.dcd     # Downloaded stripped trajectory
```

</details>

### 8.2. PyMOL Trajectory Visualization

Run these commands in the project directory on a computer with desktop PyMOL.
For a Colab run, first download the processed trajectories using `cmda pull --traj`.
Open one variant at a time; close its PyMOL window before running the next command.

```bash
cmda view --state apo --variant WT --replica r1 --frame all --resi 1
```

PyMOL opens and starts playing all saved trajectory frames. Use its playback
controls to pause, play or step through frames. `--frame all` displays the movie;
`--frame 1` displays only the first frame.

The command reads the matching `md/r1/top.prmtop` and `md/r1/traj.dcd` for the
selected state and variant. It generates `view.pdb`, `visualize.pml` and
`visualize.log` in that replica's analysis folder. The PDB has the same atoms as
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

Repeat for each completed variant and replica. The command validates matching
atom and frame counts before renaming, rejects conflicting files, and leaves raw
chunks, equilibration states and restart checkpoints in place. Older names remain
readable when working with historical runs.

</details>

## Complete project directory structure

<details>
<summary><b>View the complete workflow directory tree</b></summary>

The tree covers the workflow outputs; MODELLER and preparation tools may also
leave auxiliary files in their owning folders. Raw archives may use `.pdb`
instead of `.cif` when the archive fallback is used.

```text
/content/drive/MyDrive/research/projects/trpcage-example/
└── data/
    └── apo/
        ├── structure/
        │   └── 1L2Y/
        │       ├── raw/
        │       │   ├── 1l2y.cif                     # Experimental archive; PDB fallback is possible
        │       │   ├── template.pdb                 # Cleaned template coordinates
        │       │   └── uniprot.fasta                # Sequence used for modeling
        │       ├── STRUCTURE_DONE                   # Completed fetch and mapping
        │       └── variants/
        │           └── WT/
        │               ├── WT.pdb                   # Coordinates passed to preparation
        │               └── MODEL_DONE               # Completed model verification
        ├── simulation/
        │   └── WT/
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
        │       │   └── eq.pdb                       # Final equilibrated coordinates
        │       └── md/
        │           └── r1/
        │               ├── system.xml               # Copied OpenMM system
        │               ├── solvated.pdb             # Copied topology coordinates
        │               ├── npt.state.xml            # Starting state copied from equilibration
        │               ├── prod.chk                 # Checkpoint used to resume production
        │               ├── chunks/
        │               │   ├── 0001.dcd             # First committed trajectory chunk
        │               │   └── 0001.log             # Matching committed thermodynamic log
        │               ├── traj.dcd     # Autoimaged trajectory without bulk solvent
        │               ├── top.prmtop          # Matching topology; use with the trajectory above
        │               ├── trajCa.dcd              # Cα trajectory
        │               ├── topCa.prmtop                # Matching Cα topology
        │               ├── md.log                   # Copy of the merged production log
        │               └── MERGE_DONE               # Verified merged outputs
        └── analysis/
            └── WT/
                ├── eq/
                │   ├── qc.txt                       # Equilibration diagnostics
                │   └── qc.png                       # Thermodynamic and structural QC figure
                └── r1/
                    ├── rmsd.csv                     # RMSD values
                    ├── rg.csv                       # Radius of gyration values
                    ├── rmsf.csv                     # Per-residue fluctuations
                    ├── figures/
                    │   ├── rmsd_vs_time.png         # RMSD plot
                    │   ├── rg_vs_time.png           # Radius of gyration plot
                    │   └── rmsf_per_residue.png     # Fluctuation plot
                    └── ANALYSIS_DONE                # Verified replica tables and figures
```

</details>

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
