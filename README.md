# ColabMDA


**ColabMDA** is a specialized tool that lets you run high-quality Molecular Dynamics simulations on Google Colab without the fear of losing your work. The biggest problem with Colab is that it disconnects, often destroying hours of simulation data. **ColabMDA fixes this** with a "resume-safe" system that automatically saves your progress to Google Drive. If your session expires, you can resume exactly where you left off with one simple command. From modeling protein mutations to generating publication-ready analysis, ColabMDA handles the complex setup for you, making it the easiest way to get high-quality MD results using free cloud GPUs.

📖 **Full Documentation:** Visit our official manual at [colabmda.readthedocs.io](https://colabmda.readthedocs.io/)

## 🛠 Project Information
| Category | Details |
| :--- | :--- |
| **Release** | [![GitHub tag](https://img.shields.io/github/v/tag/paulshamrat/ColabMDA)](https://github.com/paulshamrat/ColabMDA/tags) |
| **Availability** | [![GitHub](https://img.shields.io/badge/GitHub-Repository-blue?logo=github)](https://github.com/paulshamrat/ColabMDA) |
| **Documentation** | [![Documentation Status](https://readthedocs.org/projects/colabmda/badge/?version=latest)](https://colabmda.readthedocs.io/en/latest/?badge=latest) |
| **Workflows** | [![Python CI](https://github.com/paulshamrat/ColabMDA/actions/workflows/python-test.yml/badge.svg?branch=main)](https://github.com/paulshamrat/ColabMDA/actions/workflows/python-test.yml) |
| **Issues** | [![GitHub issues](https://img.shields.io/github/issues/paulshamrat/ColabMDA)](https://github.com/paulshamrat/ColabMDA/issues) |
| **License** | [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) |
| **Style / Lint** | [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black) [![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff) |
| **Dependencies** | `OpenMM`, `Modeller`, `MDAnalysis`, `MDTraj` |
| **Platform** | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/paulshamrat/ColabMDA/blob/main/notebooks/05-colabmd-simulation-2024.ipynb) `Linux` `HPC (SLURM)` |
| **Structure** | [Source Layout](#-project-structure) |

### 📂 Project Structure
ColabMDA is organized into clear, functional modules:
*   `src/colabmda/modeller/`: Homology modeling engine (Biological numbering supported).
*   `src/colabmda/openmm_pw/`: OpenMM simulation engines (Modular EM/NVT/NPT/MD).
*   `envs/`: Automated installation scripts for scientific environments.
*   `scripts/`: Quick-start bootstrap scripts for Google Colab.
*   `notebooks/`: Ready-to-use Colab notebooks for simulation and analysis.

---

## 1. Quick Start (Google Colab Installation)

> 💡 **Terminal Access:** All bash commands should be run in the **Colab Terminal** (Open via the **⋮** menu -> **Terminal**).

### 1.1. Setup Colab Runtime & Drive
Before starting, ensure your environment is ready:
1. **Enable GPU:** Go to `Runtime` -> `Change runtime type` and select **T4 GPU**.
2. **Verify GPU:** Run `!nvidia-smi` in a cell to confirm GPU access.
3. **Mount Drive:** Click the **Folder icon** 📂 in the left sidebar, then click the **Drive icon** (Mount Drive), or run the code block below:

```ipython
from google.colab import drive
drive.mount('/content/drive')
!nvidia-smi
```

### 1.2. Environment & Package Installation (Required)
Run the following in the **Colab Terminal** (⋮ → Terminal). 
*Estimated time: ~3–5 minutes.*

```bash
# 1. Install the core scientific environment
cd /content
curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/main/scripts/bootstrap_colab_openmm_gpu.sh -o bootstrap_colab_openmm_gpu.sh
WITH_MODELLER=1 bash bootstrap_colab_openmm_gpu.sh latest

# 2. Install ColabMDA package
python3 -m pip install --upgrade "git+https://github.com/paulshamrat/ColabMDA.git@main"
```

> [!IMPORTANT]
> **Modeller License Prompt:** During Step 1, the script will **pause** and ask you to `Enter your Modeller License Key`. You must paste your key and press **Enter** to proceed. The installation will not complete without it.
>
> 🔑 **Get a Free License:** If you don't have one, register at [salilab.org/modeller/registration.html](https://salilab.org/modeller/registration.html) (Academic licenses are free and sent instantly via email).

---

### 💡 Tip: How to Resume After a Timeout
If your Google Colab session expires:
1. Re-run **Required Steps 1.2** to reinstall the environment.
2. Run the **exact same `colabmda openmm run` command** you used before.
3. The tool will automatically detect your `.chk` files and resume from where it left off.

---

## 2. Manual / Detailed Installation (Advanced)

> ⚠️ **Note:** This section is for local workstations, HPC, or advanced manual setups. For standard Google Colab runs, please use the **Quick Start (Section 1)** above.

<details>
<summary><b>🛠️ A. Manual Terminal Installation (Step-by-Step)</b></summary>

In the Colab Terminal (⋮ → Terminal), run each step one at a time:

```bash
# Step 1: Download & install Miniforge (Conda)
wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh && \
  bash /tmp/miniforge.sh -b -p "$HOME/miniforge3"

# Step 2: Initialize Conda in this shell
export PATH="$HOME/miniforge3/bin:$PATH" && source "$HOME/miniforge3/etc/profile.d/conda.sh"

# Step 3: Install Mamba into the base environment
conda install -y -n base -c conda-forge mamba

# Step 4: Install CUDA-enabled OpenMM and OpenMMTools
mamba install -y -c conda-forge cudatoolkit=11.8 openmm openmmtools

# Step 5: Install PDBFixer (conda, fallback to pip)
conda install -y -c conda-forge pdbfixer || pip install pdbfixer

# Step 6: Install MDAnalysis, MDTraj, NumPy, Matplotlib, and Biopython
mamba install -y -c conda-forge mdanalysis mdtraj numpy matplotlib biopython

# Step 7: Verify installations
python3 - << 'EOF'
from openmm import Platform; print("OpenMM platforms:", [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())])
import MDAnalysis, mdtraj, Bio; print("MDAnalysis:", MDAnalysis.__version__, "MDTraj:", mdtraj.__version__, "Biopython:", Bio.__version__)
EOF
```
</details>

<details>
<summary><b>📜 B. Alternative: Script-based Installation</b></summary>

```bash
cd /content
curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/main/scripts/install_colabmda_release.sh -o install_colabmda_release.sh
bash install_colabmda_release.sh latest /content/colabmda
```
</details>

<details>
<summary><b>🧬 C. Modeller CPU Environment Setup</b></summary>

```bash
cd /content/drive/MyDrive/openmm/ColabMDA
bash envs/install_modeller_env.sh
```
</details>

</details>

<details>
<summary><b>💻 D. Local Workstation Setup (Laptop/Desktop)</b></summary>

Beyond the Cloud ☁️: ColabMDA works on any Linux system with an NVIDIA GPU. Use the provided `environment.yml` to create a production-ready environment:
```bash
mamba env create -f environment.yml
conda activate colabmda
```
</details>

<details>
<summary><b>🏢 E. HPC Usage (SLURM)</b></summary>

You can easily incorporate ColabMDA into SLURM batch scripts. Since it processes trajectories in chunks, it is highly efficient for long-running jobs on cluster partitions with time limits.
</details>


---

## 3. Simulation Workflow

### Pipeline at a Glance:
1. **Model:** Build structures (WT & Mutants) in `structures/`
2. **Stage:** Initialize the simulation folder in `simulations/`
3. **Run:** Execute the MD simulation (Resume-safe)
4. **Merge:** Combine trajectory chunks into a final file
5. **Analyze:** Generate RMSD, Rg, and RMSF plots

### 3.1. Build Structures (WT and Mutants)
**Environment:** `modeller_env`

The build workflow now includes **Biological Numbering** and **Automated Quality Control**.

#### New Features:
*   **`--uniprot-numbering`**: Physically re-numbers the PDB residues to match the UniProt biological index (Best Practice).
*   **Automatic Alignment Summary**: Displays a full sequence comparison before building to catch range errors early.
*   **Post-Build Sanity Check**: Verifies every residue in the final PDB against the UniProt reference and reports `✅ SUCCESS` or `❌ FAILED`.

```bash
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate modeller_env
cd /path/to/your/project

# Example: Build Wild-Type KRAS (Starting at Residue 1)
colabmda modeller build --pdb-id 4ldj --uniprot-id P01116 --chain A --range 1 169 --uniprot-numbering --outdir structures/4ldj/wt

# Example: Create G12D Mutant (Preserves Numbering)
colabmda modeller mutate --pdb-in structures/4ldj/wt/target.B99990001_with_cryst.pdb --chain A --mut G12D --outdir-mut structures/4ldj/mutants/4ldj_G12D
```

### 3.2. Setup and Run MD
**Environment:** `openmm_env`

```bash
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate openmm_env
cd /path/to/your/project

# 1. Initialize the simulation folder
# For Wild-Type:
colabmda openmm stage --pdb-file structures/4ldj/wt/target.B99990001_with_cryst.pdb --name 4ldj_wt --replica r1
# For Mutant (G12D):
colabmda openmm stage --pdb-file structures/4ldj/mutants/4ldj_G12D/target.B99990001_G12D.pdb --name 4ldj_G12D --replica r1

# 2. Start the Production Run (Example: 10ns)
# For Wild-Type:
colabmda openmm run --name 4ldj_wt --replica r1 --total-ns 10.0 --traj-interval 10 --equil-time 1000 --checkpoint-ps 1000
# For Mutant (G12D):
colabmda openmm run --name 4ldj_G12D --replica r1 --total-ns 10.0 --traj-interval 10 --equil-time 1000 --checkpoint-ps 1000
```

> 💡 **Storage Tip:** For a typical system (e.g., KRAS in water, ~30,000 atoms), a 100ns run at high resolution (1ps) can produce over **36GB** of data. On a free 15GB Google Drive, we recommend using **`--traj-interval 10`** to reduce this to ~3.6GB. Always calculate your storage needs based on your specific system size before starting long runs.

> **Note:** The `run` command includes an **Automated Stability Gate**. It automatically analyzes equilibration logs and aborts if the system hasn't stabilized, saving GPU time.

### 3.3. Merge and Center
Combine trajectory chunks into a single DCD, apply periodic boundary condition (PBC) correction, and center the protein using the robust MDAnalysis-based engine (`--mda`).

```bash
# Standard Merge (Center + Wrap, All Atoms)
# For Wild-Type:
colabmda openmm merge --pdb-dir simulations/4ldj_wt/r1 --center --wrap

# For Mutant (G12D):
colabmda openmm merge --pdb-dir simulations/4ldj_G12D/r1 --center --wrap
```

#### Output Files:
After a standard merge, the following files are created in the simulation replica folder:
*   `prod_full.dcd`: The final, concatenated, and centered/wrapped trajectory file.
*   `prod_full.log`: The consolidated log file containing energy and temperature statistics.

> 💡 **Pro-Tip for Long Runs:**
> Merging processes trajectories frame-by-frame, so it won't crash your RAM. You can merge without striding (`--stride 1`) for full resolution, or use `--stride 10` to create a lightweight file for local viewing.
> 
> When you merge with `--stride 10`, the logs are automatically strided to match, and the analysis command (`colabmda openmm analysis`) will read the logs to correctly infer the time scale.

> 💡 **Protein-Only Trajectory (Optional):**  
> If you wish to save significant disk space (saving >85% storage), you can add the `--protein-only` flag to the merge command. This will extract only the protein atoms, discarding water and ions:  
> `colabmda openmm merge --pdb-dir simulations/4ldj_wt/r1 --center --wrap --protein-only`
> 
> **Outputs with `--protein-only`:**
> * `prod_full.dcd`: Thinned trajectory containing protein atoms only.
> * `prod_full.pdb`: Matching protein-only topology PDB file (crucial for subsequent analysis/visualization to avoid atom mismatch).
> * `prod_full.log`: Consolidated log file.

### 3.4. Post-Merge Guidelines & Verification (FAQ)

#### 3.4.1. 🔍 Verifying Trajectory Frames
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

Alternatively, you can query the frame count using a python one-liner with MDAnalysis:
```bash
python3 -c "import MDAnalysis as mda; u = mda.Universe('prod_full.pdb', 'prod_full.dcd'); print('Frames:', len(u.trajectory))"
```

#### 3.4.2. ⏱️ Understanding Simulation Time vs. Frame Counts
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
    The `status` command will display this discrepancy clearly:
    ```text
    Frames (from logs): 10000        # Original simulation frames
    Merged DCD        : YES (1000 frames) # Thinned frames after striding
    ```

#### 3.4.3. 🧬 Reference Topology Guidelines
When merging trajectories, the reference topology file used and the resulting output files depend on your merge options:

| Merge Mode | Command Flags | Reference Topology | Generated Topology File | Subsequent Command Usage |
| :--- | :--- | :--- | :--- | :--- |
| **Standard Merge** | `colabmda openmm merge --center --wrap` | `solvated.pdb` | *None* (only `prod_full.dcd` is written) | Use `solvated.pdb` for analysis/visualization. |
| **MDAnalysis Merge** | `colabmda openmm merge --mda --center --wrap` | `solvated.pdb` | `prod_full.pdb` (all atoms) | Use `prod_full.pdb` (or `solvated.pdb`). |
| **Protein-Only Merge** | `colabmda openmm merge --protein-only` | `solvated.pdb` | `prod_full.pdb` (protein atoms only) | **Must** use `prod_full.pdb` (since `prod_full.dcd` contains only protein coordinates). |

##### Why does this matter?
If you merge using the `--protein-only` flag, your `prod_full.dcd` will only contain protein coordinates (~2,600 atoms). Attempting to load this trajectory along with the original `solvated.pdb` (~27,000 atoms) in PyMOL or MDAnalysis will result in a **fatal atom mismatch error**. Always match the trajectory file with its corresponding topology file as shown in the table above.


---

## 4. Analysis & Comparison

### 4.1. Single System Analysis
```bash
# For Wild-Type:
colabmda openmm analysis --pdb-id 4ldj_wt
# For Mutant (G12D):
colabmda openmm analysis --pdb-id 4ldj_G12D
```

> ⚠️ **Analysis Tip:** If your plots show the wrong time scale (e.g., 10ns instead of 100ns), provide the frame interval manually. For example, if you ran with `--traj-interval 10`:
> `colabmda openmm analysis --pdb-id 4ldj_wt --interval 10`

### 4.2. WT vs Mutant Comparison
```bash
colabmda openmm compare \
  --series "WT=analysis/single/4ldj_wt/r1,analysis/single/4ldj_wt/r2" \
  --series "G12D=analysis/single/4ldj_G12D/r1,analysis/single/4ldj_G12D/r2" \
  --outdir analysis/compare/wt_vs_g12d_avg
```

---

## 5. Trajectory Visualization in PyMOL

To visualize molecular dynamics trajectories in PyMOL with full secondary structure (ribbon/cartoon) and detailed sidechain representations, you can use the built-in `colabmda openmm view` command.

### 5.1. Local Setup and Installation
Since molecular dynamics simulations are performed on Google Colab, trajectory visualization is run locally on your workstation/laptop. 

> [!NOTE]
> **For Windows Users:** It is highly recommended to run the visualization locally either using **WSL (Windows Subsystem for Linux)** or using native **PyMOL on Windows** (by running the CLI to generate `visualize.pml`, then opening it in the Windows PyMOL GUI).

It is highly recommended to create a dedicated conda environment (e.g. `colabmda_env`) containing both PyMOL and the ColabMDA package:

```bash
# 1. Create a clean environment and install PyMOL
conda create -y -n colabmda_env -c conda-forge python=3.11 pymol-open-source

# 2. Activate the environment
conda activate colabmda_env

# 3. Install ColabMDA package (Lightweight local installation)
python3 -m pip install --upgrade "git+https://github.com/paulshamrat/ColabMDA.git@main"

# 4. Verify PyMOL installation and version
pymol --version
```
*(Note: This local installation is extremely lightweight and does not require heavy simulation engines like OpenMM or MDAnalysis just to view files).*

### 5.2. Running the Visualization in PyMOL
Simply run the view command from your simulation folder (it will automatically look for `prod_full.pdb` and `prod_full.dcd`, generate a PyMOL script, and launch PyMOL):
```bash
# 1. Run view command from within your replica directory:
cd simulations/4ldj_wt/r1
colabmda openmm view

# 2. Or run it by specifying the directory path:
# For Wild-Type:
colabmda openmm view --pdb-dir simulations/4ldj_wt/r1

# For Mutant (e.g., G12D):
colabmda openmm view --pdb-dir simulations/4ldj_G12D/r1
```

> 💡 **Custom Residue/Topology/Trajectory:**  
> By default, the tool highlights residue index 12 (mutant site). You can customize the highlighted residue and load custom trajectory files using flags:  
> `colabmda openmm view --pdb-dir simulations/4ldj_wt/r1 --resi 12 -t prod_full.pdb -x prod_full.dcd`

The command automatically generates a `visualize.pml` script in the folder. 

<details>
<summary><b>Click to view the generated PyMOL script configuration (<code>visualize.pml</code>)</b></summary>

```pymol
# 1. Clear out all the old overlapping objects from memory
reinitialize
bg_color white

# 2. Load the merged trajectory files (prod_full.pdb and prod_full.dcd)
load prod_full.pdb, kras
load_traj prod_full.dcd, kras

# 2b. Align trajectory to the crystal reference to ensure identical viewing orientation (e.g., 4ldj_wt.pdb)
load 4ldj_wt.pdb, crystal_ref
align kras and name CA, crystal_ref and name CA, mobile_state=1, target_state=1
delete crystal_ref

# 3. Hide solvent water and ions
hide everything, all
hide nonbonded, all
hide nb_spheres, all

# 4. Freeze the backbone tumbling rotation frame-by-frame
intra_fit kras and name CA

# 5. Generate your crisp secondary structure ribbon
dss kras
cartoon automatic, kras
show cartoon, kras

# 6. Apply your smooth cyan color scheme
color cyan, kras and name C*
util.cnc("kras")

# 7. Highlight your mutated Cysteine 12 side chain sticks flawlessly
show sticks, kras and resi 12 and not name N+C+O+H
color yellow, kras and resi 12 and name SG
set stick_radius, 0.25

# 8. Focus camera right onto the protein
zoom kras and polymer, buffer=4
```

</details>

### 5.3. Generating Comparative Snapshot Grids
To render publication-quality structural snapshot grids comparing transitions over simulation trajectories (such as WT vs Mutants across specific frames), you can use the `colabmda openmm snapshots` command.

> ⚠️ **Prerequisite:** This command requires the Python `pymol` and `Pillow` (PIL) libraries to be installed in the active environment (e.g., `pymol-viz`).

```bash
# Generate the default 3x8 transition snapshot grid (WT vs G12C vs G12D)
colabmda openmm snapshots
```

#### Customization:
By default, the command uses a built-in template designed for the KRAS WT/G12C/G12D trajectory transition grid. You can customize the behavior by supplying a custom JSON configuration file:

```bash
colabmda openmm snapshots --config my_config.json --output figures/comparison_grid.png
```

#### JSON Configuration Schema:
```json
{
  "align_ref_pdb": "structures/4ldj/wt/target.B99990001_with_cryst.pdb",
  "stable_core_sel": "resi 1-10 or resi 40-55 or resi 80-169",
  "camera_view": [
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0,
    1.0, 0.0, 0.0,
    0.0, 0.0, -50.0
  ],
  "systems": {
    "WT": {
      "pdb": "simulations/4ldj_wt/r1/prod_full.pdb",
      "dcd": "simulations/4ldj_wt/r1/prod_full.dcd",
      "states": [1, 200],
      "times": ["0.00 ns", "2.00 ns"],
      "mut_residue": 12
    }
  }
}
```

---

## 6. Project Strategy
Organize work in three phases:
1. **Preparation**: Build WT first in `structures/<pdbid>/wt/`, then generate mutants.
2. **Simulation**: Run WT and mutants in separate folders under `simulations/`.
3. **Analysis**: Store per-system analysis in `analysis/single/`, then generate overlays in `analysis/compare/`.

---

## 7. Project Structure
```text
/content/drive/MyDrive/openmm/
  structures/
    4ldj/
      wt/          # Wild-type modeled PDBs
      mutants/     # G12D/G12C modeled PDBs
  simulations/
    4ldj_wt/
      r1/          # Replica 1 (em.chk, npt.chk, prod.dcd)
      r2/          # Replica 2
    4ldj_G12D/
      r1/
      r2/
  analysis/
    single/
      4ldj_wt/     # [r1, r2, aggregate] reports
      4ldj_G12D/
    compare/       # Final WT vs Mutant overlays
```

## Acknowledgements
- [OpenMM](https://openmm.org) & [PDBFixer](https://github.com/openmm/pdbfixer)
- [Modeller](https://salilab.org/modeller/)
- [MDAnalysis](https://www.mdanalysis.org) & [MDTraj](https://www.mdtraj.org)
- [NumPy](https://numpy.org), [Matplotlib](https://matplotlib.org), [Biopython](https://biopython.org)
- [Google Colab](https://colab.research.google.com) & [Miniforge/Conda](https://github.com/conda-forge/miniforge)

## Changelog

### v0.1.0 (Initial Beta)
*   **Modular Pipeline:** New modular CLI for EM, NVT, NPT, and Production MD.
*   **Resume-Safe Engine:** Integrated checkpointing logic for fail-safe simulations on Google Colab.
*   **Modeling:** Automated Wild-Type building and mutation support via Modeller.
*   **Analysis:** Robust trajectory merging and comparative RMSD/Rg/RMSF analysis tools.
*   **Professional Standards:** Added CI/CD workflows, Black formatting, and Ruff linting.

---

## Citation

Initially, this repository was based on a GROMACS-on-Colab workflow, but we transformed it to utilize OpenMM to improve feasibility, performance, and stability for automated simulation recovery on cloud platforms.

If you use this tool, please consider citing the underlying studies:

> Paul SK, Saddam M, Rahaman KA, Choi JG, Lee SS, Hasan M. **Molecular modeling, molecular dynamics simulation, and essential dynamics analysis of grancalcin: An upregulated biomarker in experimental autoimmune encephalomyelitis mice.** *Heliyon*. 2022 Oct 23;8(10):e11232. doi: [10.1016/j.heliyon.2022.e11232](https://doi.org/10.1016/j.heliyon.2022.e11232). PMID: 36340004; PMCID: PMC9626934.
>
> [![Cited by](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fapi.openalex.org%2Fworks%2Fhttps%3A%2F%2Fdoi.org%2F10.1016%2Fj.heliyon.2022.e11232&query=%24.cited_by_count&label=cited%20by&color=blue)](https://doi.org/10.1016/j.heliyon.2022.e11232)
>
> Paul SK, Saddam M, Tabassum N, Hasan M. **Molecular dynamics simulation of wild and mutant proteasome subunit beta type 8 (PSMB8) protein: Implications for restoration of inflammation in experimental autoimmune encephalomyelitis pathogenesis.** *Heliyon*. 2024 Dec 15;11(1):e41166. doi: [10.1016/j.heliyon.2024.e41166](https://doi.org/10.1016/j.heliyon.2024.e41166). PMID: 39802026; PMCID: PMC11719297.
>
> [![Cited by](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fapi.openalex.org%2Fworks%2Fhttps%3A%2F%2Fdoi.org%2F10.1016%2Fj.heliyon.2024.e41166&query=%24.cited_by_count&label=cited%20by&color=blue)](https://doi.org/10.1016/j.heliyon.2024.e41166)

---

<p align="center">
  <img src="https://komarev.com/ghpvc/?username=paulshamrat&repo=ColabMDA&label=Visitors&color=blue&style=flat" alt="Visitors" />
</p>
