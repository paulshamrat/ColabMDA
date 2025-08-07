# ColabMDA: Google Colaboratory–Based Molecular Dynamics Simulation & Analysis

![Flowchart](https://github.com/paulshamrat/ColabMDA/blob/main/images/flowchart.png)

---


## 🚩 Update July 2025

> ⚠️ **CMake Issue (2025) Resolved:**
>
> The standard GROMACS installation is currently unavailable due to CMake version issues. However, you can seamlessly run molecular dynamics simulations using the robust OpenMM Colab workflow as an alternative protocol.
>
> 👉 **See the detailed OpenMM Colab workflow below or [`openmm/README.md`](openmm/README.md) for a complete guide.**

### Notebooks for PSMB8 Publication

These notebooks were developed for the published PSMB8 study (see below):

| Status         | Notebook                          | Description                                                | Colab Link |
|---------------|------------------------------------|------------------------------------------------------------|------------|
| ✅ Simulation  | `05-colabmd-simulation-2024.ipynb` | Execute MD runs for PSMB8 wild-type & G210V mutant         | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/paulshamrat/ColabMDA/blob/main/notebooks/05-colabmd-simulation-2024.ipynb) |
| ✅ Analysis    | `03-colabmd-analysis.ipynb`        | Process & visualize trajectories using MDAnalysis & MDTraj  | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/paulshamrat/ColabMDA/blob/main/notebooks/03-colabmd-analysis.ipynb) |
| ❌ Retired     | `04-colab-gmx-install.ipynb`       | (Retired) GROMACS installation with GPU support on Colab; replaced by the OpenMM protocol below | — |










# OpenMM Protein-Water MD Workflow (Summary)

This workflow enables GPU-accelerated molecular dynamics simulations of proteins in explicit water using OpenMM and PDBFixer on Google Colab. The protocol covers:

- Environment setup and package installation in Colab Terminal
- Downloading all required workflow scripts
- Preprocessing PDB structures
- Running chunked MD simulations with checkpoint/restart support
- Merging trajectory chunks
- Trajectory analysis and visualization

For full details, step-by-step instructions, and example commands, please see [`openmm/openmm_proteinwater/readme_proteinwater.md`](openmm/openmm_proteinwater/readme_proteinwater.md).

---

## Want to Cite Us?

If you use this workflow or data in your research, please cite:

Paul, S. K., Saddam, M., Tabassum, N., & Hasan, M. (2024). Molecular dynamics simulation of wild and mutant proteasome subunit beta type 8 (PSMB8) protein: Implications for restoration of inflammation in experimental autoimmune encephalomyelitis pathogenesis. Heliyon, 11(1), e41166. https://doi.org/10.1016/j.heliyon.2024.e41166

For full details and scripts, see [`openmm/README.md`](openmm/README.md).

----


## Overview

This repository provides end-to-end Jupyter notebooks and scripts for running and analyzing molecular dynamics (MD) simulations in Google Colaboratory, supporting both GROMACS (original protocol) and OpenMM (alternative protocol for GPU-accelerated MD on Colab). It includes:

- **GROMACS workflow** (for PSMB8 publication):
  - Installation of GROMACS 2023.x on Colab
  - MD simulation of PSMB8 wild-type and G210V mutant
  - Trajectory analysis using MDAnalysis and MDTraj
- **OpenMM workflow** (recommended for current Colab GPU usage):
  - GPU-accelerated MD simulation and analysis with OpenMM
  - Checkpoint/restart, RMSD, RMSF, and Rg analysis

All notebooks accompany the published study:

> **Molecular dynamics simulation of wild and mutant proteasome subunit beta type 8 (PSMB8) protein: Implications for restoration of inflammation in experimental autoimmune encephalomyelitis pathogenesis**  
> _Heliyon_, 11 (2025) e41166 • [https://doi.org/10.1016/j.heliyon.2024.e41166](https://www.sciencedirect.com/science/article/pii/S2405844024171976)

---

### Repository

🔗 https://github.com/paulshamrat/ColabMDA

---

### Data & Code Archive

- **Simulation Dataset**  
  Dataset of MD simulations for PSMB8 (3UNF) and its G210V mutant in EAE  
  [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8070983.svg)](https://zenodo.org/records/8157201)

- **Supplementary Files**  
  All input files and analysis scripts used in the paper are included under the `notebooks/` and `data/` directories.

---


### Note & Acknowledgements

We would like to thank the authors who developed the Jupyter notebook framework for molecular dynamics simulation on Google Colab, as well as the OpenMM development team for their open-source molecular simulation toolkit. Please always refer to the original GROMACS and OpenMM manuals for simulation guidance. We are grateful to the authors of the following articles and software, which made it possible to adapt this MD simulation protocol:

- OpenMM: Eastman, P., et al. (2017). OpenMM 7: Rapid development of high performance algorithms for molecular dynamics. PLOS Computational Biology, 13(7), e1005659.
- Engelberger F. et al., J. Chem. Educ. 98(5):1801–1807 (2021)
- Lemkul J. A., Living J. Comput. Mol. Sci. 1(1):5068 (2019)
- Arantes P. R. et al., J. Chem. Inf. Model. 61(10):4852–4856 (2021)
- Gowers R. J. et al., Proc. 15th Python in Science Conf., 98–105 (2016)
- Abraham M. J. et al., SoftwareX 1:19–25 (2015)

---

_Last tested on: 2025-07-14_

![Visitor Badge](https://visitor-badge.laobi.icu/badge?page_id=paulshamrat.ColabMDA)


