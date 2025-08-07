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
     python3 openmm_proteinwater_chunked_auto.py <pdbid> \
       --total-ns 1 --traj-interval 10.0 --equil-time 10.0 --checkpoint-ps 10.0
     # Example:
     python3 openmm_proteinwater_chunked_auto.py 4ldj \
       --total-ns 1 --traj-interval 10.0 --equil-time 10.0 --checkpoint-ps 10.0
     ```
   - Output: Multiple trajectory chunks (`prod_*.dcd`), logs, checkpoints, and final system files in `<pdbid>/`.

3. **Merge Trajectory Chunks**
   - After simulation, merge all trajectory chunks into a single trajectory and log file.
   - Command:
     ```bash
     python3 openmm_trajmerge.py <pdbid> \
       --topology <pdbid>/solvated.pdb \
       --out-traj prod_full.dcd \
       --out-log prod_full.log
     # Example:
     python3 openmm_trajmerge.py 4ldj \
       --topology 4ldj/solvated.pdb \
       --out-traj prod_full.dcd \
       --out-log prod_full.log
     ```
   - Output: Merged trajectory (`prod_full.dcd`) and log (`prod_full.log`).

4. **Analyze Trajectory**
   - Perform analysis on the merged trajectory using the analysis script.
   - Command:
     ```bash
     python3 openmm_trajanalysis.py <pdbid> \
       --interval 1.0 \
       --topology <pdbid>/solvated.pdb \
       --trajectory <pdbid>/prod_full.dcd \
       --outdir analysis_<pdbid>
     # Example:
     python3 openmm_trajanalysis.py 4ldj \
       --interval 10.0 \
       --topology 4ldj/solvated.pdb \
       --trajectory 4ldj/prod_full.dcd \
       --outdir analysis_4ldj
     ```
   - Output: Analysis results in `analysis_<pdbid>/`.


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
