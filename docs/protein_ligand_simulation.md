# Protein-Ligand Simulation Guide

This guide describes how to model, setup, run, and analyze Molecular Dynamics (MD) simulations for protein-ligand complexes in **ColabMDA**.

> 💡 **Where to Run the Commands:**
> *   **Staging, Running, and Merging (Sections 1-4.4):** These commands are typically run on **Google Colab** (either natively in the **Colab Terminal** [⋮ → Terminal] or inside notebook cells by prefixing the commands with an exclamation mark `!`).
> *   **Trajectory Visualization (Section 4.5):** PyMOL trajectory viewing and snapshot generation are run **locally** on your workstation/laptop (which has a graphical display).
>
> For details on setting up these environments, see the [Installation & Setup Guide](installation.md).

---

## 1. Physical Model & Force Fields

ColabMDA uses a dual-forcefield setup to parameterize protein-ligand systems:
*   **Protein & Ions:** AMBER14SB (`amber14-all.xml` + `amber14/tip3p.xml`)
*   **Small Molecule Ligand:** General AMBER Force Field (GAFF2, version `gaff-2.11`)
*   **Ligand Charges:** AM1-BCC charge model generated on-the-fly using the OpenFF Toolkit
*   **Water Model:** TIP3P
*   **Cofactors/Metal Ions:** Standard AMBER parameters (e.g. for $Mg^{2+}$)

---

## 2. Installation & Setup

### 2.1. Google Colab (Cloud Setup)
To run a protein-ligand simulation on Google Colab, you must install the additional small-molecule parameterization packages (`openff-toolkit` and `openmmforcefields`). Enable this by setting `WITH_LIGAND=1` when running the bootstrap setup script:

```bash
# 1. Install scientific environment with ligand support
cd /content
curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/pwpl/scripts/bootstrap_colab_openmm_gpu.sh -o bootstrap_colab_openmm_gpu.sh
WITH_LIGAND=1 WITH_MODELLER=1 bash bootstrap_colab_openmm_gpu.sh latest

# 2. Install ColabMDA package
python3 -m pip install --upgrade "git+https://github.com/paulshamrat/ColabMDA.git@pwpl"
```

### 2.2. Local Workstation / HPC
For local workstations or cluster environments, activate your standard OpenMM Conda environment and install the required dependencies:

```bash
conda activate openmm_env
mamba install -y -c conda-forge openff-toolkit openmmforcefields
```

---

## 3. Directory Layout & Modeling Strategy

### 3.1. Recommended Structure Folder Layout
For protein-ligand complexes, we recommend organizing your files under the standard `structures/` folder:

```text
structures/
  4ldj/
    wt/
      target.B99990001_with_cryst.pdb   # Modeled WT protein structure (no ligand/MG)
    4ldj_orig.pdb                       # Raw crystal structure (used to extract MG cofactors)
    gdp.sdf                             # Ligand file in SDF format (with correct bond orders)
```

### 3.2. Why use SDF for the ligand?
Standard PDB files do not store bond order or formal charge information. Parameterizing small molecules with AM1-BCC requires explicit chemical identities, which are preserved in **SDF (Structure Data File)** formats. Always convert your ligand to SDF (using tools like Open Babel, RDKit, or PyMOL) before running simulations.

### 3.3. Ensuring Native Coordinate Alignment
To guarantee that the ligand is placed exactly in its native binding pocket:
1. **Coordinate Preservation:** When extracting the ligand (e.g., GDP) from the raw crystal PDB (`4ldj_orig.pdb`) to save as an SDF file, you must **preserve its exact 3D coordinates** without moving or translation.
2. **Automated Coordinate Alignment:** During `colabmda modeller build`, the engine automatically aligns the finished homology model (and any mutants) back to the exact coordinate space of the template crystal structure PDB (`4ldj_orig.pdb`).
3. **Perfect Merging:** Because the modeller engine aligns the protein coordinates back to the original crystal frame, the protein pocket and the SDF ligand share the identical spatial coordinate system. When OpenMM merges them during the setup phase, they align perfectly at the native coordinate binding site.

---

## 4. Step-by-Step Simulation Workflow

### 4.1. Step 1: Build Protein Structures (WT and Mutants)
**Environment:** `modeller_env`

Build your protein coordinates with **Biological Numbering** and automated QC:

```bash
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate modeller_env

# 1. Build Wild-Type KRAS (Starting at Residue 1 matching UniProt indexing)
colabmda modeller build --pdb-id 4ldj --uniprot-id P01116 --chain A --range 1 169 --uniprot-numbering --outdir structures/4ldj/wt

# 2. Create Mutant (G12D) structure from Wild-Type template (if simulating mutants)
colabmda modeller mutate --pdb-in structures/4ldj/wt/target.B99990001_with_cryst.pdb --chain A --mut G12D --outdir-mut structures/4ldj/mutants/4ldj_G12D
```

### 4.2. Step 2: Stage the Protein Topology
**Environment:** `openmm_env`

First, stage the modeled protein structure (which contains no ligand or crystal water). Name the simulation folder using the naming convention `<pdbid>_<cofactor>` (e.g. `4ldj_gdp`):

```bash
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate openmm_env

colabmda openmm stage --pdb-file structures/4ldj/wt/target.B99990001_with_cryst.pdb --name 4ldj_gdp --replica r1
```
*Note: This creates the simulation directory `simulations/4ldj_gdp/r1/` and copies the protein PDB into it.*

### 4.3. Step 3: Run the Simulation (EM, NVT, NPT, MD)
**Environment:** `openmm_env`

Run the simulation by supplying the path to the ligand SDF file and adding the `--keep-mg` flag (if your system contains active Magnesium cofactors to be preserved from the raw crystal structure):

```bash
# Production MD with GDP ligand and Magnesium preservation
colabmda openmm run --name 4ldj_gdp --replica r1 --ligand structures/4ldj/gdp.sdf --keep-mg --total-ns 100.0 --traj-interval 10 --equil-time 1000 --checkpoint-ps 1000
```

> 💡 **Magnesium Preservation (`--keep-mg`):**  
> Default PDBFixer runs strip out all heterogens (including structural metal ions like $Mg^{2+}$). The `--keep-mg` flag bypasses this by extracting the Magnesium atoms from the raw crystal structure PDB (`structures/4ldj/4ldj_orig.pdb` or the raw file inside the simulation workspace) and re-inserting them into the cleaned protein structure prior to parameterization.

### 4.4. Step 4: Merge and Center Trajectories
**Environment:** `openmm_env`

After the simulation finishes, concatenate the trajectory chunks and center the protein-ligand complex:

```bash
colabmda openmm merge --pdb-dir simulations/4ldj_gdp/r1 --center --wrap
```

### 4.5. Step 5: Trajectory Visualization in PyMOL
**Environment:** Local workstation / conda environment with PyMOL

To visualize the trajectory of your protein-ligand system, run the view command from your local workstation:

```bash
colabmda openmm view --pdb-dir simulations/4ldj_gdp/r1
```

You can view the ligand and its interactions inside PyMOL using a custom script setup. Below is the generated `visualize.pml` configuration for protein-ligand systems:

```python
# 1. Clear out old objects
reinitialize
bg_color white

# 2. Load the complex trajectory
load prod_full.pdb, complex
load_traj prod_full.dcd, complex

# 3. Represent protein as cartoon/ribbon
hide everything, all
show cartoon, complex and polymer
color cyan, complex and polymer and name C*
util.cnc("complex")

# 4. Represent Ligand (GDP) as sticks
show sticks, complex and organic
color yellow, complex and organic and name C*
util.cnc("complex and organic")

# 5. Represent Magnesium (MG) cofactors as spheres
show spheres, complex and resn MG
color green, complex and resn MG
set sphere_scale, 0.35, complex and resn MG

# 6. Fit camera and zoom onto binding pocket
zoom complex and (organic or resn MG), buffer=6
```

---

## 5. How It Works (Under the Hood)

During the energy minimization (`em`) stage:

### 5.1. Magnesium Extraction
If `--keep-mg` is active, ColabMDA reads the raw starting structure PDB and filters out coordinates corresponding to `MG` residues. It appends these coordinate records back to the cleaned protein structure before passing it to `PDBFixer`.

### 5.2. OpenFF Loading
The SDF ligand file is parsed using the OpenFF Toolkit's `Molecule` class to retrieve correct connectivity, aromaticity, and formal charges.

### 5.3. System Generator & Solvation
An OpenMM `SystemGenerator` is instantiated with the AMBER14SB forcefield and the GAFF2 small-molecule template. The ligand's topology and positions are merged with the protein-magnesium complex using `Modeller`. The composite system is then solvated with TIP3P water and $0.15 \text{ M}$ neutralizing $\text{NaCl}$ ions.

### 5.4. Parameterization
`SystemGenerator` parameterizes the combined structure, registering GAFF2 parameters for the ligand and AMBER14SB parameters for the protein and ions on-the-fly. The system is then serialized to `system.xml` and minimized.
