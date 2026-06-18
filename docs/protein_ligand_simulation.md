# Protein-Ligand Simulation Guide

This guide describes how to model, setup, run, and analyze Molecular Dynamics (MD) simulations for protein-ligand complexes in water using **ColabMDA**.

> 💡 **Where to Run the Commands:**
> *   **Staging, Running, Merging, and Analysis (Sections 1-4):** These commands are run on **Google Colab** (either natively in the **Colab Terminal** [⋮ → Terminal] or inside a remote shell using the Colab CLI `colab console -s kras-sim`).
> *   **Trajectory Visualization (Section 5):** PyMOL trajectory viewing and snapshot generation are run **locally** on your workstation/laptop (which has a graphical display).
>
> For details on setting up these environments, see the [Installation & Setup Guide](installation.md).

---

## 1. Simulation Workflow Overview

The standard simulation workflow consists of:
1. **Model:** Build structures (Wild-Type & Mutants) in `structures/` using Modeller.
2. **Stage:** Initialize the simulation directory structure under `simulations/`.
3. **Run:** Execute the resume-safe production MD (EM, NVT, NPT, MD stages) with GAFF2 ligand parameterization and $Mg^{2+}$ cofactor preservation.
4. **Merge:** Concatenate raw trajectory chunks, apply periodic boundary condition (PBC) wrap/alignment, and center the complex.
5. **Analyze:** Compute and plot RMSD, Radius of Gyration ($R_g$), and RMSF for the complex.

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

### 2.2. Stage Simulation Workspace
**Environment:** `openmm_env`

Stage the WT simulation workspace to create the directory structure and copy the cleaned starting structure:
```bash
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate openmm_env

# Stage simulation directory
colabmda openmm stage --pdb-file structures/4ldj/wt/target.B99990001_with_cryst.pdb --name 4ldj_gdp --replica r1
```
*Note: This creates the simulation directory `simulations/4ldj_gdp/r1/` and copies the protein PDB into it.*

> 💡 **Tip: Unified One-Command Execution (Alternative)**
> If you prefer to run the entire simulation sequentially (EM, NVT, NPT, QC, and MD) using a single, resume-safe command instead of executing each step manually, you can run:
> ```bash
> colabmda openmm run --name 4ldj_gdp --replica r1 --ligand structures/4ldj/gdp.sdf --keep-mg --total-ns 10.0 --traj-interval 10 --equil-time 1000 --checkpoint-ps 1000
> ```

### 2.3. Energy Minimization (EM)
**Environment:** `openmm_env`

Stages the ligand and cofactor, parameterizes them, solvates, and minimizes potential energy to resolve steric clashes:
```bash
colabmda openmm em --name 4ldj_gdp --workdir simulations/4ldj_gdp/r1 --ligand structures/4ldj/gdp.sdf --keep-mg
```
*Under the Hood:* Parameterizes the protein with **AMBER14SB**, the ligand with **GAFF2** (using **AM1-BCC** charges via the OpenFF Toolkit), and water with the **TIP3P** model. Solvates with a 1.0 nm cubic box and adds $0.15 \text{ M}$ NaCl. Outputs `solvated.pdb`, `system.xml`, and the minimized state `em.chk`.

### 2.4. NVT Equilibration
**Environment:** `openmm_env`

Relaxes solvent molecules and adjusts temperature to 300 K while keeping the protein-ligand complex restrained:
```bash
colabmda openmm nvt --name 4ldj_gdp --workdir simulations/4ldj_gdp/r1 --equil-time 100 --seed 42
```
*Parameters:* `--equil-time` (length in ps), `--seed` (random seed for velocity assignment). Temperature is maintained at $300 \text{ K}$ using a Langevin Middle Integrator (collision frequency $1.0 \text{ ps}^{-1}$).

### 2.5. NPT Equilibration
**Environment:** `openmm_env`

Relaxes the simulation box volume and adjusts system density under pressure:
```bash
colabmda openmm npt --name 4ldj_gdp --workdir simulations/4ldj_gdp/r1 --equil-time 100
```
*Parameters:* Pressure is maintained at $1.0 \text{ atm}$ using an OpenMM Monte Carlo Barostat (volume adjustment attempts every 25 steps).

### 2.6. Equilibration Validation & QC Check
**Environment:** `openmm_env`

Plots and validates the stability of thermodynamic parameters (temperature, density, potential energy) before starting production MD:
```bash
colabmda openmm check-equil --name 4ldj_gdp --workdir simulations/4ldj_gdp/r1
```

### 2.7. Production Molecular Dynamics (MD)
**Environment:** `openmm_env`

Runs the production trajectory in chunked, resume-safe segments:
```bash
colabmda openmm md --name 4ldj_gdp --workdir simulations/4ldj_gdp/r1 --total-ns 10.0 --traj-interval 10 --checkpoint-ps 1000
```
*Parameters:* `--total-ns` (total duration in ns), `--traj-interval` (coordinate saving frequency in ps), `--checkpoint-ps` (duration of each chunk before writing checkpoint to enable safe restarts). All bonds involving hydrogen atoms are constrained using the SHAKE-like algorithm (`HBonds`), allowing a stable 2 fs timestep. Long-range electrostatics are handled via Particle Mesh Ewald (PME) with a $1.0\text{ nm}$ cutoff.

> 💡 **Multi-Replica Acceleration (r2, r3, ...):**
> To run multiple independent replicas without repeating the CPU-intensive solvation, minimization, and equilibration steps, you can directly spawn them from `r1`'s equilibrated state. Simply run the unified `run` command for the next replica:
> ```bash
> colabmda openmm run --name 4ldj_gdp --replica r2 --total-ns 10.0
> ```
> This automatically skips the EM/NVT/NPT preparation steps, inherits `system.xml`, `solvated.pdb`, and the initial NPT checkpoint from `r1`, and initializes the new simulation with randomized velocities (using a new random seed) to ensure independent trajectories.

### 2.8. Merge and Center
**Environment:** `openmm_env`

Combine trajectory chunks into a single DCD, apply periodic boundary condition (PBC) correction, and center the protein-ligand complex:
```bash
colabmda openmm merge --pdb-dir simulations/4ldj_gdp/r1 --center --wrap
```

#### Output Files:
After a standard merge, the following files are created in the simulation replica folder:
*   `prod_full.dcd`: The final, concatenated, and centered/wrapped trajectory file.
*   `prod_full.log`: The consolidated log file containing energy and temperature statistics.

---

## 3. Post-Merge Guidelines & Verification (FAQ)

### 3.1. Directory Layout & Modeling Strategy

#### Recommended Structure Folder Layout
For protein-ligand complexes, we recommend organizing your files under the standard `structures/` folder:

```text
structures/
  4ldj/
    wt/
      target.B99990001_with_cryst.pdb   # Modeled WT protein structure (no ligand/MG)
    4ldj_orig.pdb                       # Raw crystal structure (used to extract MG cofactors)
    gdp.sdf                             # Ligand file in SDF format (with correct bond orders)
```

#### Why use SDF for the ligand?
Standard PDB files do not store bond order or formal charge information. Parameterizing small molecules with AM1-BCC requires explicit chemical identities, which are preserved in **SDF (Structure Data File)** formats. Always convert your ligand to SDF (using tools like Open Babel, RDKit, or PyMOL) before running simulations.

#### Ensuring Native Coordinate Alignment
To guarantee that the ligand is placed exactly in its native binding pocket:
1. **Coordinate Preservation:** When extracting the ligand (e.g., GDP) from the raw crystal PDB (`4ldj_orig.pdb`) to save as an SDF file, you must **preserve its exact 3D coordinates** without moving or translation.
2. **Automated Coordinate Alignment:** During `colabmda modeller build`, the engine automatically aligns the finished homology model (and any mutants) back to the exact coordinate space of the template crystal structure PDB (`4ldj_orig.pdb`).
3. **Perfect Merging:** Because the modeller engine aligns the protein coordinates back to the original crystal frame, the protein pocket and the SDF ligand share the identical spatial coordinate system. When OpenMM merges them during the setup phase, they align perfectly at the native coordinate binding site.

### 3.2. 🔍 Verifying Trajectory Frames
To quickly check the status, simulation time, and number of frames in your merged files, run:
```bash
colabmda openmm status --pdb-dir simulations/4ldj_gdp/r1
```

### 3.3. ⏱️ Understanding Simulation Time vs. Frame Counts
To calculate how much simulation time (in nanoseconds) your trajectory represents or how many frames you should expect, use this quick guide:

*   **Total Simulation Time**: Controlled by `--total-ns` (e.g., `100.0` ns = `100,000` ps).
*   **Frame Saving Frequency**: Controlled by `--traj-interval` in picoseconds (default is `10.0` ps = `0.01` ns).
*   **Calculating Expected Frames**:
    $$\text{Expected Frames} = \frac{\text{Total Time (ps)}}{\text{Trajectory Saving Interval (ps)}}$$

### 3.4. 🧬 Reference Topology Guidelines
When merging trajectories, the reference topology file used and the resulting output files depend on your merge options:

| Merge Mode | Command Flags | Reference Topology | Generated Topology File | Subsequent Command Usage |
| :--- | :--- | :--- | :--- | :--- |
| **Standard Merge** | `colabmda openmm merge --center --wrap` | `solvated.pdb` | *None* (only `prod_full.dcd` is written) | Use `solvated.pdb` for analysis/visualization. |
| **MDAnalysis Merge** | `colabmda openmm merge --mda --center --wrap` | `solvated.pdb` | `prod_full.pdb` (all atoms) | Use `prod_full.pdb` (or `solvated.pdb`). |

---

## 4. Analysis & Comparison

### 4.1. Single System Analysis
```bash
colabmda openmm analysis --pdb-id 4ldj_gdp
```

### 4.2. WT vs Mutant Comparison
```bash
colabmda openmm compare --series "WT=analysis/single/4ldj_gdp/r1" --outdir analysis/compare/wt_vs_mutant_avg
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
colabmda openmm view --pdb-dir simulations/4ldj_gdp/r1
```

Inside PyMOL, the ligand and magnesium ions are automatically rendered with distinct representations (sticks and spheres, respectively) for easy binding-pocket inspection:

```python
# Generated visualize.pml for protein-ligand systems:
reinitialize
bg_color white
load prod_full.pdb, complex
load_traj prod_full.dcd, complex
hide everything, all
show cartoon, complex and polymer
color cyan, complex and polymer and name C*
util.cnc("complex")
show sticks, complex and organic
color yellow, complex and organic and name C*
util.cnc("complex and organic")
show spheres, complex and resn MG
color green, complex and resn MG
set sphere_scale, 0.35, complex and resn MG
zoom complex and (organic or resn MG), buffer=6
```

### 5.3. Superimposition Check & Alignment Verification (PyMOL)
To run the panel superimposition grid analysis locally on your workstation, execute:
```bash
colabmda openmm snapshots
```

#### Alignment Verification Script
Below is the python script utilized by the superimposition check:

```python
import os
from PIL import Image, ImageOps, ImageDraw, ImageFont
import pymol
from pymol import cmd

# Initialize PyMOL
pymol.finish_launching(['pymol', '-qc'])
cmd.reinitialize()
cmd.bg_color('white')

# Set standard manuscript camera matrix (preserves identical scale and framing)
camera_view = [
    -0.597478449,   -0.617595792,   -0.511463583,
    -0.502994895,    0.785388291,   -0.360776752,
     0.624511778,    0.041705273,   -0.779902637,
     0.000050262,   -0.000109358, -182.850265503,
    31.930021286,   32.076828003,   34.52369679,
   141.553222656,  224.141143799,  -20.000000000
]

stable_core_sel = 'name CA and (resi 1-24 or resi 41-54 or resi 81-166)'

# Helper function to align and render a system with exact camera view
def render_panel(name, is_crystal=False, is_model=False, is_superimposed=False, is_ligand_pocket=False):
    cmd.reinitialize()
    cmd.bg_color('white')
    
    # Load reference structure (for alignment only)
    cmd.load("structures/g12c_protein.pdb", "ref_g12c")
    
    # Unified Load & Align:
    cmd.load("structures/4ldj/wt/4ldj_orig.pdb", "crystal")
    cmd.align(f"crystal and {stable_core_sel}", f"ref_g12c and {stable_core_sel}")
    
    cmd.load("structures/4ldj/wt/target.B99990001_with_cryst.pdb", "modeled_wt")
    cmd.matrix_copy("crystal", "modeled_wt")
    
    cmd.load("structures/4ldj/gdp.sdf", "transferred_gdp")
    cmd.matrix_copy("crystal", "transferred_gdp")
    
    cmd.create("mg_cofactor", "crystal and resn MG")
    
    # Selective Clean-up
    cmd.delete("ref_g12c")
    if is_crystal:
        cmd.delete("modeled_wt")
        cmd.delete("transferred_gdp")
        cmd.delete("mg_cofactor")
    elif is_model:
        cmd.delete("crystal")
    elif is_ligand_pocket:
        cmd.delete("crystal")
    elif is_superimposed:
        pass

    # Hide all helper structures
    cmd.hide("everything", "all")
    
    # Show primary structures
    if is_crystal or is_superimposed:
        cmd.show("cartoon", "crystal and chain A")
        cmd.dss("crystal")
    if is_model or is_superimposed or is_ligand_pocket:
        cmd.show("cartoon", "modeled_wt and chain A")
        cmd.dss("modeled_wt")
    
    # Stylings
    if is_crystal:
        cmd.color("gray90", "crystal and chain A")
        cmd.color("orange", "crystal and resid 57-75")
        cmd.color("yellow", "crystal and resid 12")
        cmd.show("sticks", "crystal and resid 12 and not name N+C+O+H")
        
        cmd.show("sticks", "crystal and resn GDP")
        cmd.color("wheat", "crystal and resn GDP and name C*")
        cmd.show("spheres", "crystal and resn MG")
        cmd.color("forest", "crystal and resn MG")
        cmd.set("sphere_scale", 0.35, "crystal and resn MG")
        
    elif is_model:
        cmd.color("cyan", "modeled_wt and chain A")
        cmd.color("marine", "modeled_wt and resid 57-75")
        cmd.color("yellow", "modeled_wt and resid 12")
        cmd.show("sticks", "modeled_wt and resid 12 and not name N+C+O+H")
        
        cmd.show("sticks", "transferred_gdp")
        cmd.color("magenta", "transferred_gdp and name C*")
        cmd.show("spheres", "mg_cofactor")
        cmd.color("forest", "mg_cofactor")
        cmd.set("sphere_scale", 0.35, "mg_cofactor")
        
    elif is_ligand_pocket:
        cmd.color("cyan", "modeled_wt and chain A")
        cmd.color("marine", "modeled_wt and resid 57-75")
        cmd.set("cartoon_transparency", 0.5, "modeled_wt")
        
        pocket_resi_sel = "modeled_wt and (resid 16 or resid 17 or resid 35 or resid 57 or resid 119 or resid 146)"
        cmd.show("sticks", pocket_resi_sel)
        cmd.color("cyan", pocket_resi_sel + " and name C*")
        cmd.show("labels", pocket_resi_sel + " and name CA")
        cmd.set("label_color", "black")
        cmd.set("label_size", 14)
        
        cmd.show("sticks", "transferred_gdp")
        cmd.color("magenta", "transferred_gdp and name C*")
        cmd.set("stick_radius", 0.3, "transferred_gdp")
        
        cmd.show("spheres", "mg_cofactor")
        cmd.color("forest", "mg_cofactor")
        cmd.set("sphere_scale", 0.45, "mg_cofactor")
        
    elif is_superimposed:
        cmd.color("gray90", "crystal and chain A")
        cmd.color("orange", "crystal and resid 57-75")
        cmd.color("yellow", "crystal and resid 12")
        cmd.set("cartoon_transparency", 0.6, "crystal")
        
        cmd.color("cyan", "modeled_wt and chain A")
        cmd.color("marine", "modeled_wt and resid 57-75")
        cmd.color("yellow", "modeled_wt and resid 12")
        cmd.set("cartoon_transparency", 0.0, "modeled_wt")
        
        cmd.show("sticks", "(crystal or modeled_wt) and resid 12 and not name N+C+O+H")
        
        cmd.show("sticks", "crystal and resn GDP")
        cmd.color("wheat", "crystal and resn GDP and name C*")
        cmd.set("stick_radius", 0.18, "crystal and resn GDP")
        
        cmd.show("sticks", "transferred_gdp")
        cmd.color("magenta", "transferred_gdp and name C*")
        cmd.set("stick_radius", 0.28, "transferred_gdp")
        
        cmd.show("spheres", "mg_cofactor")
        cmd.color("forest", "mg_cofactor")
        cmd.set("sphere_scale", 0.35, "mg_cofactor")

    # Common stick settings
    cmd.set("stick_radius", 0.25, "resid 12")
    if is_crystal or is_superimposed:
        cmd.set("stick_radius", 0.18, "crystal and resn GDP")
    if is_model or is_superimposed or is_ligand_pocket:
        cmd.set("stick_radius", 0.28, "transferred_gdp")
    cmd.color("blue", "name N*")
    cmd.color("red", "name O*")
    cmd.color("orange", "name P*")
    cmd.color("white", "name H*")

    cmd.set('ray_shadows', 1)
    cmd.set('ray_trace_mode', 1)
    cmd.set('ray_trace_gain', 0.5)
    cmd.set('depth_cue', 1)
    cmd.set('specular', 0.3)
    cmd.set('ray_opaque_background', 1)

    if is_ligand_pocket:
        cmd.set_view(camera_view)
        cmd.zoom("transferred_gdp or mg_cofactor", buffer=4.0)
    else:
        cmd.set_view(camera_view)
        
    cmd.png(f"panel_{name}.png", width=1200, height=1200, ray=1)

# Render 4 panels
render_panel("A", is_crystal=True)
render_panel("B", is_model=True)
render_panel("C", is_ligand_pocket=True)
render_panel("D", is_superimposed=True)
cmd.quit()

# Crop and combine
cropped_images = {}
for name in ["A", "B", "C", "D"]:
    img = Image.open(f"panel_{name}.png")
    gray = img.convert('L')
    inverted = ImageOps.invert(gray)
    bbox = inverted.getbbox()
    if bbox:
        padding = 40
        w, h = img.size
        left = max(0, bbox[0] - padding)
        top = max(0, bbox[1] - padding)
        right = min(w, bbox[2] + padding)
        bottom = min(h, bbox[3] + padding)
        
        side = max(right - left, bottom - top)
        cx, cy = (left + right) // 2, (top + bottom) // 2
        left = max(0, cx - side // 2)
        top = max(0, cy - side // 2)
        right = min(w, left + side)
        bottom = min(h, top + side)
        
        cropped = img.crop((left, top, right, bottom))
        cropped_images[name] = cropped.resize((800, 800), Image.Resampling.LANCZOS)

composite = Image.new("RGB", (1600, 1600), "white")
composite.paste(cropped_images["A"], (0, 0))
composite.paste(cropped_images["B"], (800, 0))
composite.paste(cropped_images["C"], (0, 800))
composite.paste(cropped_images["D"], (800, 800))

draw = ImageDraw.Draw(composite)

try:
    font = ImageFont.truetype("arial.ttf", 32)
except IOError:
    font = ImageFont.load_default()

labels = {
    "A": (30, 30, "A. Crystal Template (PDB: 4LDJ)"),
    "B": (830, 30, "B. Homology Model (WT)"),
    "C": (30, 830, "C. Active Site & Transferred Ligand"),
    "D": (830, 830, "D. Superimposed Alignment (RMSD: 0.132 A)")
}

for panel, (x, y, text) in labels.items():
    if hasattr(draw, "textsize"):
        tw, th = draw.textsize(text, font=font)
    else:
        left, top, right, bottom = font.getbbox(text)
        tw = right - left
        th = bottom - top
    draw.rectangle([x-10, y-5, x+tw+10, y+th+5], fill="white")
    draw.text((x, y), text, fill="black", font=font)

composite.save("kras_superimposed.png")
```

#### Quantitative Metrics
You can calculate the $C_\alpha$ Root Mean Square Deviation (RMSD) between the crystal template and the homology model in python to quantitatively verify coordinate preservation:

```python
from Bio.PDB import PDBParser
import numpy as np

parser = PDBParser(QUIET=True)
c_struct = parser.get_structure('c', 'structures/4ldj/wt/4ldj_orig.pdb')[0]['A']
m_struct = parser.get_structure('m', 'structures/4ldj/wt/target.B99990001_with_cryst.pdb')[0]['A']

# Extract C-alpha coordinates for matching residues (1 to 169)
c_coords = np.array([r['CA'].get_coord() for r in c_struct if r.id[0] == ' ' and r.id[1] >= 1 and r.id[1] <= 169])
m_coords = np.array([r['CA'].get_coord() for r in m_struct if r.id[0] == ' ' and r.id[1] >= 1 and r.id[1] <= 169])

diff = c_coords - m_coords
rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
print(f"C-alpha RMSD: {rmsd:.4f} Angstroms")
```

For a correct homology model built using `colabmda modeller build` against template `4LDJ`:
* **All-atom $C_\alpha$ RMSD**: **`0.1322 Å`** (Confirms near-perfect structural backbone conservation).
* **Ligand/Mg²⁺ Pocket Fit**: **`0.0000 Å`** deviation (The transferred GDP ligand and $Mg^{2+}$ sphere coordinates overlap perfectly with the crystallographic positions inside the binding pocket).
