# Protein-in-Water Simulation Guide

This guide describes how to run Molecular Dynamics (MD) simulations for protein systems in water using **ColabMDA**.

> 💡 **Where to Run the Commands:**
> *   **Staging, Running, Merging, and Analysis (Sections 1-4):** These commands are typically run on **Google Colab** (either natively in the **Colab Terminal** [⋮ → Terminal] or inside notebook cells by prefixing the commands with an exclamation mark `!`).
> *   **Trajectory Visualization (Section 5):** PyMOL trajectory viewing and snapshot generation are run **locally** on your workstation/laptop (which has a graphical display).
>
> For details on setting up these environments, see the [Installation & Setup Guide](installation.md).

---

## 1. Simulation Workflow Overview

The standard simulation workflow consists of:
1. **Model:** Build structures (Wild-Type & Mutants) in `structures/` using Modeller.
2. **Stage:** Initialize the simulation directory structure under `simulations/`.
3. **Run:** Execute the resume-safe production MD (EM, NVT, NPT, MD stages).
4. **Merge:** Concatenate raw trajectory chunks, align/wrap boundaries, and center the protein.
5. **Analyze:** Compute and plot RMSD, Radius of Gyration ($R_g$), and RMSF.

---

## 2. Step-by-Step Execution

### 2.1. Build Structures (WT and Mutants)
**Environment:** `modeller_env`

Build your protein coordinates with **Biological Numbering** and automated QC:

```bash
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate modeller_env

# 1. Build Wild-Type KRAS (Starting at Residue 1 matching UniProt indexing)
colabmda modeller build --pdb-id 4ldj --uniprot-id P01116 --chain A --range 1 169 --uniprot-numbering --outdir structures/4ldj/wt

# 2. Create Mutant (G12D) structure from Wild-Type template
colabmda modeller mutate --pdb-in structures/4ldj/wt/target.B99990001_with_cryst.pdb --chain A --mut G12D --outdir-mut structures/4ldj/mutants/4ldj_G12D
```

### 2.2. Setup and Run MD
**Environment:** `openmm_env`

```bash
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate openmm_env

# 1. Stage simulation directories
colabmda openmm stage --pdb-file structures/4ldj/wt/target.B99990001_with_cryst.pdb --name 4ldj_wt --replica r1
colabmda openmm stage --pdb-file structures/4ldj/mutants/4ldj_G12D/target.B99990001_G12D.pdb --name 4ldj_G12D --replica r1

# 2. Start the Production Run (e.g., 10ns production)
colabmda openmm run --name 4ldj_wt --replica r1 --total-ns 10.0 --traj-interval 10 --equil-time 1000 --checkpoint-ps 1000
colabmda openmm run --name 4ldj_G12D --replica r1 --total-ns 10.0 --traj-interval 10 --equil-time 1000 --checkpoint-ps 1000
```

> 💡 **Storage Tip:** For a typical system (e.g., KRAS in water, ~30,000 atoms), a 100ns run at high resolution (1ps) can produce over **36GB** of data. On a free 15GB Google Drive, we recommend using **`--traj-interval 10`** to reduce this to ~3.6GB. Always calculate your storage needs based on your specific system size before starting long runs.

> **Note:** The `run` command includes an **Automated Stability Gate**. It automatically analyzes equilibration logs and aborts if the system hasn't stabilized, saving GPU time.

### 2.3. Merge and Center
Combine trajectory chunks into a single DCD, apply periodic boundary condition (PBC) correction, and center the protein using the robust MDAnalysis-based engine (`--mda`).

```bash
colabmda openmm merge --pdb-dir simulations/4ldj_wt/r1 --center --wrap
colabmda openmm merge --pdb-dir simulations/4ldj_G12D/r1 --center --wrap
```

#### Output Files:
After a standard merge, the following files are created in the simulation replica folder:
*   `prod_full.dcd`: The final, concatenated, and centered/wrapped trajectory file.
*   `prod_full.log`: The consolidated log file containing energy and temperature statistics.

---

## 3. Post-Merge Guidelines & Verification (FAQ)

### 3.1. 🔍 Verifying Trajectory Frames
To quickly check the status, simulation time, and number of frames in your merged files, run:
```bash
colabmda openmm status --pdb-dir simulations/4ldj_wt/r1
```

This will print a comprehensive status report including topology stats, chunks, and exact frame counts:
```text
[STATUS]
  Workdir          : /path/to/simulations/4ldj_wt/r1
  Chunks (DCD/log) : 100 / 100
  Topology File    : solvated.pdb
                     └─ 27273 atoms, 8388 residues
                     └─ (169 protein, 8170 water, 49 ions)
  Trajectory File  : prod_full.dcd (10000 frames)
  Log File         : prod_full.log (10000 frames)
  Frames (from logs): 10000
```

### 3.2. ⏱️ Understanding Simulation Time vs. Frame Counts
To calculate how much simulation time (in nanoseconds) your trajectory represents or how many frames you should expect, use this quick guide:

*   **Total Simulation Time**: Controlled by `--total-ns` (e.g., `100.0` ns = `100,000` ps).
*   **Frame Saving Frequency**: Controlled by `--traj-interval` in picoseconds (default is `10.0` ps = `0.01` ns).
*   **Calculating Expected Frames**:
    $$\text{Expected Frames} = \frac{\text{Total Time (ps)}}{\text{Trajectory Saving Interval (ps)}}$$
    *   *Example:* If you run a `100.0` ns simulation with a `10.0` ps saving interval:
        $$\text{Expected Frames} = \frac{100,000\text{ ps}}{10\text{ ps}} = 10,000\text{ frames}$$

*   **Effect of Striding on Merged Trajectories**:
    If you merge chunks using a stride (e.g., `--stride 10` for lightweight local viewing):
    $$\text{Merged Frames} = \frac{\text{Total Frames}}{\text{Stride}} = \frac{10,000}{10} = 1,000\text{ frames}$$

### 3.3. 🧬 Reference Topology Guidelines
When merging trajectories, the reference topology file used and the resulting output files depend on your merge options:

| Merge Mode | Command Flags | Reference Topology | Generated Topology File | Subsequent Command Usage |
| :--- | :--- | :--- | :--- | :--- |
| **Standard Merge** | `colabmda openmm merge --center --wrap` | `solvated.pdb` | *None* (only `prod_full.dcd` is written) | Use `solvated.pdb` for analysis/visualization. |
| **MDAnalysis Merge** | `colabmda openmm merge --mda --center --wrap` | `solvated.pdb` | `prod_full.pdb` (all atoms) | Use `prod_full.pdb` (or `solvated.pdb`). |
| **Protein-Only Merge** | `colabmda openmm merge --protein-only` | `solvated.pdb` | `prod_full.pdb` (protein atoms only) | **Must** use `prod_full.pdb` (since `prod_full.dcd` contains only protein coordinates). |

---

## 4. Analysis & Comparison

### 4.1. Single System Analysis
```bash
colabmda openmm analysis --pdb-id 4ldj_wt
colabmda openmm analysis --pdb-id 4ldj_G12D
```

### 4.2. WT vs Mutant Comparison
```bash
colabmda openmm compare --series "WT=analysis/single/4ldj_wt/r1,analysis/single/4ldj_wt/r2" --series "G12D=analysis/single/4ldj_G12D/r1,analysis/single/4ldj_G12D/r2" --outdir analysis/compare/wt_vs_g12d_avg
```

---

## 5. Trajectory Visualization in PyMOL

### 5.1. Local Setup and Installation
It is highly recommended to create a dedicated conda environment containing PyMOL and ColabMDA:

```bash
conda create -y -n colabmda_env -c conda-forge python=3.11 pymol-open-source
conda activate colabmda_env
python3 -m pip install --upgrade "git+https://github.com/paulshamrat/ColabMDA.git@pwpl"
```

### 5.2. Running the Visualization in PyMOL
Simply run the view command from your simulation replica folder:
```bash
colabmda openmm view --pdb-dir simulations/4ldj_wt/r1
```

### 5.3. Generating Comparative Snapshot Grids
To render publication-quality structural snapshot grids:

```bash
colabmda openmm snapshots
```
