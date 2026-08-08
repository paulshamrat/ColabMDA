# Protein-in-Water Simulation Guide

This guide describes how to run Molecular Dynamics (MD) simulations for protein systems in water using **ColabMDA**.

> 💡 **Where to Run the Commands:**
> *   **Staging, Running, Merging, and Analysis (Sections 1-4):** These commands are typically run on **Google Colab** (either natively in the **Colab Terminal** [⋮ → Terminal] or inside notebook cells by prefixing the commands with an exclamation mark `!`).
> *   **Trajectory Visualization (Section 5):** PyMOL trajectory viewing and snapshot generation are run **locally** on your workstation/laptop (which has a graphical display).
>
> For details on setup, see the [Installation Guide](installation.md). For the scientific
> definitions of EM/NVT/NPT/production timing, checkpoint continuation, and independent
> replicas, see [Simulation Protocol, Timing, and Replica Semantics](simulation_protocol.md).

---

## 1. Simulation Workflow Overview

The standard simulation workflow consists of:
1. **Model:** Build structures (Wild-Type & Mutants) in `data/str/` using Modeller.
2. **Stage:** Initialize the simulation directory structure under `data/sim/`.
3. **Run:** Execute the resume-safe production MD (EM, NVT, NPT, MD stages).
4. **Merge:** Concatenate raw trajectory chunks, align/wrap boundaries, and center the protein.
5. **Analyze:** Compute and plot RMSD, Radius of Gyration ($R_g$), and RMSF.

---

## 2. Step-by-Step Execution

### 2.1. Build & Prepare Structures (Modeller + PDB2PQR / PROPKA Titration)
**Environments:** `modeller_env` (for structure building) & `openmm_env` (for pKa titration)

#### 1. Build Coordinates (Biological Numbering)
```bash
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate modeller_env

# Build Wild-Type KRAS (Starting at Residue 1 matching UniProt indexing)
colabmda modeller build --pdb-id 4ldj --uniprot-id P01116 --chain A --range 1 169 --uniprot-numbering --outdir data/str/4ldj/proteins/4ldj_wt

# Create Mutant (G12D) structure from Wild-Type template
colabmda modeller mutate --pdb-in data/str/4ldj/proteins/4ldj_wt/4ldj_wt.pdb --chain A --mut G12D --outdir-mut data/str/4ldj/proteins/4ldj_g12d
```

#### 2. Protonate & Titrate Structures at pH 7.4 (PDB2PQR / PROPKA)
```bash
conda activate openmm_env

# Perform 3D pKa titration at pH 7.4 for WT and Mutant structures in data/str/
colabmda openmm prep --pdb-file data/str/4ldj/proteins/4ldj_wt/4ldj_wt.pdb --outdir data/str/4ldj/proteins/4ldj_wt/prep --name 4ldj_wt --ph 7.4
colabmda openmm prep --pdb-file data/str/4ldj/proteins/4ldj_g12d/4ldj_g12d.pdb --outdir data/str/4ldj/proteins/4ldj_g12d/prep --name 4ldj_g12d --ph 7.4
```
*Under the Hood:* PDB2PQR and PROPKA evaluate local 3D hydrogen-bonding networks and assign explicit AMBER titration forms (`HID`, `HIE`, `HIP`, `ASH`, `GLH`, `LYN`, `CYX`). An audit log `protonation_summary.json` is generated alongside each cleaned starting PDB in `data/str/`.

### 2.2. Stage Simulation Workspace
**Environment:** `openmm_env`

Stage the WT and mutant simulation workspaces in `data/sim/` using the PROPKA-protonated starting structures from `data/str/`:
```bash
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate openmm_env

# Stage simulation directories under data/sim/ using the cleaned, titratable starting structures
colabmda openmm stage --pdb-file data/str/4ldj/proteins/4ldj_wt/prep/4ldj_wt_cleaned.pdb --name 4ldj_wt --replica r1
colabmda openmm stage --pdb-file data/str/4ldj/proteins/4ldj_g12d/prep/4ldj_g12d_cleaned.pdb --name 4ldj_g12d --replica r1
```

> 💡 **Tip: Unified One-Command Execution (Alternative)**
> If you prefer to run the entire simulation sequentially (EM, NVT, NPT, QC, and MD) using a single, resume-safe command instead of executing each step manually, you can run:
> ```bash
> colabmda openmm run --name 4ldj_wt --replica r1 --total-ns 10.0 --traj-interval 10 --checkpoint-ps 1000
> ```

### 2.3. Energy Minimization (EM)
**Environment:** `openmm_env`

Stages the water molecules, parameterizes the protein-water system, and minimizes potential energy to resolve steric clashes:
```bash
colabmda openmm em --name 4ldj_wt --workdir data/sim/4ldj_wt/r1
```
*Under the Hood:* Parameterizes the protein with **AMBER ff19SB** and **OPC** water, neutralizes the system, and uses PME with a 1.0 nm cutoff and 0.8-1.0 nm switching. `--padding-nm auto` normally selects 1.0 nm for KRAS; use `1.4` for conservative 14 A padding. Outputs include `solvated.pdb`, `system.xml`, `parameterization.json`, and portable `em.state.xml`.

### 2.4. NVT Equilibration
**Environment:** `openmm_env`

Relaxes solvent molecules and adjusts temperature to 300 K while keeping protein heavy atoms restrained:
```bash
colabmda openmm nvt --name 4ldj_wt --workdir data/sim/4ldj_wt/r1
```
*Under the Hood:* Runs the default staged NVT profile: restrained minimization, gradual 50-300 K heating, post-heating minimization, and a 10/5/2/1 restraint-release ladder at 300 K. See the [stage table](simulation_protocol.md#default-staged-equilibration-profile).

### 2.5. NPT Equilibration
**Environment:** `openmm_env`

Relaxes the simulation box volume and adjusts system density under pressure:
```bash
colabmda openmm npt --name 4ldj_wt --workdir data/sim/4ldj_wt/r1
```
*Under the Hood:* Runs stages 19-24 at 300 K and **1 bar**: four restrained 100 ps NPT stages, 100 ps unrestrained NPT, then 1 ns unrestrained NPT.

> 💡 **Random Seed Note:** The `--seed` option is completely optional. If omitted, ColabMDA automatically generates a secure, cryptographically random seed (`fresh_seed`). For multi-replica runs (`r1`, `r2`, `r3`), ColabMDA automatically derives a unique random seed per replica (`derive_replica_seed`). Specify `--seed <integer>` only when you wish to force exact numerical reproducibility for a specific run.

### 2.6. Equilibration Validation & QC Check
**Environment:** `openmm_env`

Plots and validates the stability of thermodynamic parameters (temperature, density, potential energy) before starting production MD:
```bash
colabmda openmm check-equil --name 4ldj_wt --workdir data/sim/4ldj_wt/r1
```

### 2.7. Production Molecular Dynamics (MD)
**Environment:** `openmm_env`

Runs the production trajectory in chunked, resume-safe segments:
```bash
colabmda openmm md --name 4ldj_wt --workdir data/sim/4ldj_wt/r1 --total-ns 10.0 --traj-interval 10 --checkpoint-ps 1000
```
*Parameters:* `--total-ns` (total duration in ns), `--traj-interval` (coordinate saving frequency in ps), `--checkpoint-ps` (duration of each chunk before writing checkpoint to enable safe restarts). All bonds involving hydrogen atoms are constrained using the SHAKE-like algorithm (`HBonds`), allowing a stable 2 fs timestep. Long-range electrostatics are handled via Particle Mesh Ewald (PME) with a $1.0\text{ nm}$ cutoff.

> 💡 **Force field options:** Examples use the default `--protein-ff ff19SB --water-model opc`. For compatibility tests, you can choose alternatives such as `--protein-ff ff14SB --water-model tip3p`, but keep Amber19/OPC for the standard KRAS examples unless you intentionally want a different model.

> 💡 **Multi-Replica Acceleration & Shared Equilibration (`r1`, `r2`, `r3`, ...):**
> To run multiple independent replicas without repeating CPU/GPU-intensive solvation, minimization, and equilibration steps, `ColabMDA` supports **Shared-Equilibration Branching**:
>
> 1. **Option A: Run Equilibration Once (`--equil-only`):**
>    ```bash
>    colabmda openmm run --name 4ldj_wt --equil-only
>    ```
>    This runs EM, NVT, and NPT to produce `system.xml`, `solvated.pdb`, and `npt.state.xml`, then stops.
>
> 2. **Option B: Branch Production Replicas Instantly:**
>    ```bash
>    colabmda openmm run --name 4ldj_wt --replica r1 --total-ns 100.0
>    colabmda openmm run --name 4ldj_wt --replica r2 --total-ns 100.0
>    colabmda openmm run --name 4ldj_wt --replica r3 --total-ns 100.0
>    ```
>    Replicas auto-discover the shared equilibration state in `r1/`, `equil/`, or the parent folder (`../`), copy `system.xml` and `npt.state.xml`, assign a **unique per-replica random seed** (`derive_replica_seed`), and immediately launch production MD without wasting time re-minimizing or re-heating!

> 💡 **Storage Tip:** For a typical system (e.g., KRAS in water, ~30,000 atoms), a 100ns run at high resolution (1ps) can produce over **36GB** of data. On a free 15GB Google Drive, we recommend using **`--traj-interval 10`** to reduce this to ~3.6GB. Always calculate your storage needs based on your specific system size before starting long runs.

### 2.8. Merge and Center
**Environment:** `openmm_env`

Combine trajectory chunks into a single DCD, apply periodic boundary condition (PBC) correction, and center the protein using the robust MDAnalysis-based engine (`--mda`):
```bash
colabmda openmm merge --pdb-dir data/sim/4ldj_wt/r1 --center --wrap
colabmda openmm merge --pdb-dir data/sim/4ldj_g12d/r1 --center --wrap
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
colabmda openmm status --pdb-dir data/sim/4ldj_wt/r1
```

This will print a comprehensive status report including topology stats, chunks, and exact frame counts:
```text
[STATUS]
  Workdir          : /path/to/data/sim/4ldj_wt/r1
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
colabmda openmm analysis --pdb-id 4ldj_g12d
```

### 4.2. WT vs Mutant Comparison
```bash
colabmda openmm compare --series "WT=data/analysis/4ldj_wt/r1,data/analysis/4ldj_wt/r2" --series "G12D=data/analysis/4ldj_g12d/r1,data/analysis/4ldj_g12d/r2" --outdir data/analysis/compare/wt_vs_g12d_avg
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
colabmda openmm view --pdb-dir data/sim/4ldj_wt/r1
```

### 5.3. Generating Comparative Snapshot Grids
To render publication-quality structural snapshot grids:

```bash
colabmda openmm snapshots
```
