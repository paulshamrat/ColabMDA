# OpenMM Protein-Water MD Workflow

This folder provides a streamlined workflow for running molecular dynamics (MD) simulations of proteins in explicit water using OpenMM and PDBFixer. The workflow is designed for ease of use and reproducibility, from PDB download and cleaning to full MD production runs.

## Installation

### 1. Notebook Setup (Colab)
1. Open a new Colab notebook in Google Drive.
2. In the first cell, mount Google Drive and check GPU allocation:
   ```python
   from google.colab import drive
   drive.mount('/content/drive')
   !nvidia-smi
   ```
   *(This ensures your session has access to Google Drive and a GPU is available.)*

> **Note:** All environment setup and package installation (including Conda, Mamba, OpenMM, and analysis libraries) should be performed in the Colab Terminal, not in notebook cells.

### 2. Installation on Terminal

In the Colab Terminal (⋮ → Terminal), run each step one at a time:

```bash
## 00 installation
#Step 1: Download & install Miniforge (Conda)
wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh && bash /tmp/miniforge.sh -b -p "$HOME/miniforge3"

#Step 2: Initialize Conda in this shell
export PATH="$HOME/miniforge3/bin:$PATH" && source "$HOME/miniforge3/etc/profile.d/conda.sh"

#Step 3: Install Mamba into the base environment
conda install -y -n base -c conda-forge mamba

#Step 4: Install CUDA-enabled OpenMM and OpenMMTools
mamba install -y -c conda-forge cudatoolkit=11.8 openmm openmmtools

#Step 5: Install PDBFixer (conda, fallback to pip)
conda install -y -c conda-forge pdbfixer || pip install pdbfixer

#Step 6: Install MDAnalysis, MDTraj, NumPy, Matplotlib, and Biopython
mamba install -y -c conda-forge mdanalysis mdtraj numpy matplotlib biopython

#Step 7: Verify installations
python3 - << 'EOF'
from openmm import Platform; print("OpenMM platforms:", [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())])
import MDAnalysis, mdtraj, Bio; print("MDAnalysis:", MDAnalysis.__version__, "MDTraj:", mdtraj.__version__, "Biopython:", Bio.__version__)
EOF
```




# Download Required Scripts

After setting up your working directory, download all necessary workflow scripts using wget:

```bash
wget https://raw.githubusercontent.com/paulshamrat/ColabMDA/refs/heads/main/openmm/openmm_proteinwater/openmm_proteinwater.py
wget https://raw.githubusercontent.com/paulshamrat/ColabMDA/refs/heads/main/openmm/openmm_proteinwater/openmm_trajanalysis.py
wget https://raw.githubusercontent.com/paulshamrat/ColabMDA/refs/heads/main/openmm/openmm_proteinwater/openmm_trajmerge.py
wget https://raw.githubusercontent.com/paulshamrat/ColabMDA/refs/heads/main/openmm/openmm_proteinwater/pdbfixer_cleaning.py
wget https://raw.githubusercontent.com/paulshamrat/ColabMDA/refs/heads/main/openmm/openmm_proteinwater/openmm_rmsd.py
wget https://raw.githubusercontent.com/paulshamrat/ColabMDA/refs/heads/main/openmm/openmm_proteinwater/openmm_rg.py
wget https://raw.githubusercontent.com/paulshamrat/ColabMDA/refs/heads/main/openmm/openmm_proteinwater/openmm_rmsf.py
wget https://raw.githubusercontent.com/paulshamrat/ColabMDA/refs/heads/main/openmm/openmm_proteinwater/openmm_compare_plots.py
```

These scripts enable the full workflow: preprocessing, simulation, merging, and analysis.

## Workflow Overview

1. **Preprocess the PDB Structure**
   - Download and clean a PDB structure, removing heterogens (except water), building missing residues/atoms, and adding hydrogens at pH 7.0.
   - Command:
     ```bash
     python3 pdbfixer_cleaning.py <pdbid>
     # Example:
     python3 pdbfixer_cleaning.py 4ldj
     ```
   - Output: Creates a folder `<pdbid>/` containing `<pdbid>.pdb` (raw) and `<pdbid>_cleaned.pdb` (processed).

2. **Run Chunked MD Simulation**
   - Perform minimization, equilibration (NVT & NPT), and chunked production MD in the same folder.
   - Command:
     ```bash
     python3 openmm_proteinwater.py <pdbid> \
       --total-ns 1 --traj-interval 10.0 --equil-time 10.0 --checkpoint-ps 10.0
     # Example:
     python3 openmm_proteinwater.py 4ldj \
       --total-ns 1 --traj-interval 10.0 --equil-time 10.0 --checkpoint-ps 10.0
     ```
   - Output: Multiple trajectory chunks (`prod_*.dcd`), logs, checkpoints, and final system files in `<pdbid>/`.

3. **Merge Trajectory Chunks**
   - After simulation, merge all trajectory chunks into a single trajectory and log file.
   - Command:
     ```bash
     python3 openmm_trajmerge.py <pdbid> \
       --topology <pdbid>/solvated.pdb \
       --stride 1 \
       --out-traj prod_full.dcd \
       --out-log prod_full.log
     # Example:
     python3 openmm_trajmerge.py 4ldj_wt \
       --topology 4ldj_wt/solvated.pdb \
       --stride 1 \
       --out-traj prod_full.dcd \
       --out-log prod_full.log
     ```
   - Output: Merged trajectory (`prod_full.dcd`) and log (`prod_full.log`).
   - **Optimization Tip:** For long simulations (e.g., 100 ns), use **`--stride 10`** to keep the file size manageable (~1,000 frames total).
   - **Note:** If you use a stride during merging, you **must** multiply your analysis interval (see below).

4. **Analyze Trajectory**
   - Perform analysis on the merged trajectory. The script automatically detects the simulation length and frame interval.
   - All results are standardized to **Ångströms (Å)** with publication-quality visuals (600 DPI PNG & vector PDF).
   - **The Golden Rule for --interval:**
     `Analysis Interval (ps) = Simulation traj-interval * Merging stride`
     *Example: 10ps recording interval with stride 10 = **--interval 100*** 

   - Command:
     ```bash
     python3 openmm_trajanalysis.py <pdbid> \
       --topology <pdbid>/solvated.pdb \
       --trajectory <pdbid>/prod_full.dcd \
       --interval <calculated_value> \
       --outdir <dir_name>
     ```

   - Examples for your Project:
     ```bash
     # Wild Type (WT) - 50.3 ns trajectory (100ps recording + Stride 1)
     python3 openmm_trajanalysis.py 4ldj_wt \
       --topology 4ldj_wt/solvated.pdb \
       --trajectory 4ldj_wt/prod_full.dcd \
       --interval 100.0 \
       --outdir analysis_4ldj_wt
     
     # G12C Mutant - 100 ns target (Assuming 100ps recording + Stride 10)
     python3 openmm_trajanalysis.py 4ldj_g12c \
       --topology 4ldj_g12c/solvated.pdb \
       --trajectory 4ldj_g12c/prod_full.dcd \
       --interval 1000.0 \
       --outdir analysis_4ldj_g12c
     ```
   - Features:
     - **Trajectory-First Visuals:** Crisp, solid lines to emphasize simulation dynamics.
     - **Scientific Fixes:** Proper protein-only Rg selection and C-alpha RMSF superposition.

5. **Optional: Separate RMSD / Rg / RMSF Scripts**
   - All standalone scripts now follow the same high-fidelity Ångström standards.
   - RMSD:
     ```bash
     python3 openmm_rmsd.py <pdbid> --outdir analysis_<pdbid>_rmsd
     ```
   - Rg:
     ```bash
     python3 openmm_rg.py <pdbid> --outdir analysis_<pdbid>_rg
     ```
   - RMSF (with optional zoom controls):
     ```bash
     python3 openmm_rmsf.py <pdbid> --resid-min 2 --ylim 0 5.0 --outdir analysis_<pdbid>_rmsf
     ```

6. **Optional: Compare WT vs Mutants (Overlay Plots)**
   - After generating `rmsd.csv`, `rg.csv`, and `rmsf.csv` for each system:
     ```bash
     python3 openmm_compare_plots.py \
       --series WT=analysis_4ldj_wt \
       --series G12C=analysis_4ldj_G12C \
       --series G12D=analysis_4ldj_G12D \
       --outdir analysis_compare_wt_vs_mutants
     ```


## Acknowledgements

This workflow was made possible by the following open-source tools and libraries:

- **OpenMM**: High-performance molecular dynamics engine for biomolecular simulations ([openmm.org](https://openmm.org/)).
- **PDBFixer**: Tool for fixing and preparing PDB files for simulation ([github.com/openmm/pdbfixer](https://github.com/openmm/pdbfixer)).
- **MDAnalysis**: Python library for analysis of molecular dynamics trajectories ([www.mdanalysis.org](https://www.mdanalysis.org/)).
- **MDTraj**: Modern, open library for analyzing molecular dynamics trajectories ([mdtraj.org](https://mdtraj.org/)).
- **NumPy**: Fundamental package for scientific computing with Python ([numpy.org](https://numpy.org/)).
- **Matplotlib**: Comprehensive library for creating static, animated, and interactive visualizations in Python ([matplotlib.org](https://matplotlib.org/)).
- **Biopython**: Collection of Python tools for computational biology ([biopython.org](https://biopython.org/)).
- **Google Colab**: Free cloud-based Jupyter notebook environment with GPU support ([colab.research.google.com](https://colab.research.google.com/)).
- **Miniforge/Conda/Mamba**: Package and environment management for reproducible scientific workflows.

We gratefully acknowledge the developers and maintainers of these projects for their contributions to the scientific and open-source communities.

---
For questions or issues, please contact the repository maintainer.
