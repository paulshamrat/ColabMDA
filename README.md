# ColabMDA: Google Colaboratory–Based Molecular Dynamics Simulation & Analysis

![Flowchart](https://github.com/paulshamrat/ColabMDA/blob/main/images/flowchart.png)

---

## Update July 2025 [Ongoing]

⚠️ **Notice:** CMake issues are under investigation; GROMACS installation is a work in progress.

| Status                             | Notebook                          | Description                                                | Colab Link                                                                                                                                           |
|------------------------------------|-----------------------------------|------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------|
| 🚧 **New (WIP): Updated CMake**     | `04.1_colab_gmx_install.ipynb`      | Install & configure GROMACS 2023.x with patched CMake      | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/paulshamrat/ColabMDA/blob/main/notebooks/04.1_colab_gmx_install.ipynb) |
| ❌ **Retired: CMake Issue v072025** | `04-colab-gmx-install.ipynb`      | Legacy install failing due to CMake version mismatch       | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/paulshamrat/ColabMDA/blob/main/notebooks/04-colab-gmx-install.ipynb)    |
| ✅ **Simulation (2024)**            | `05-colabmd-simulation-2024.ipynb`| Execute MD runs for PSMB8 wild-type & G210V mutant         | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/paulshamrat/ColabMDA/blob/main/notebooks/05-colabmd-simulation-2024.ipynb) |
| ✅ **Analysis**                     | `03-colabmd-analysis.ipynb`       | Process & visualize trajectories using MDAnalysis & MDTraj | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/paulshamrat/ColabMDA/blob/main/notebooks/03-colabmd-analysis.ipynb)     |


## Other Notebooks

| Notebook                  | Description                                | Colab Link                                                                                                                                       |
|---------------------------|--------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------|
| `250417_modeller.ipynb`   | Rebuild missing residues using Modeller     | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/paulshamrat/ColabMDA/blob/main/notebooks/250417_modeller.ipynb) |


## Overview

This repository provides end-to-end Jupyter notebooks for running and analyzing GROMACS-based molecular dynamics (MD) simulations entirely in Google Colaboratory. It includes:

- **Installation** of GROMACS 2023.x on Colab  
- **MD Simulation** of PSMB8 wild-type and G210V mutant  
- **Trajectory Analysis** using MDAnalysis and MDTraj  

All notebooks accompany the published study:

> **Molecular dynamics simulation of wild and mutant proteasome subunit beta type 8 (PSMB8) protein: Implications for restoration of inflammation in experimental autoimmune encephalomyelitis pathogenesis**  
> _Heliyon_, 11 (2025) e41166 • [https://doi.org/10.1016/j.heliyon.2024.e41166](https://www-sciencedirect-com.libproxy.clemson.edu/science/article/pii/S2405844024171976)  

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

We would like to thank the authors who developed the Jupyter notebook framework for molecular dynamics simulation on Google Colab. Please always refer to the original GROMACS manual for simulation guidance. We are grateful to the authors of the following articles, which made it possible to adapt this MD simulation protocol:

- Engelberger F. et al., J. Chem. Educ. 98(5):1801–1807 (2021)  
- Lemkul J. A., Living J. Comput. Mol. Sci. 1(1):5068 (2019)  
- Arantes P. R. et al., J. Chem. Inf. Model. 61(10):4852–4856 (2021)  
- Gowers R. J. et al., Proc. 15th Python in Science Conf., 98–105 (2016)  
- Abraham M. J. et al., SoftwareX 1:19–25 (2015)

---

_Last tested on: 2025-07-14_

![Visitor Badge](https://visitor-badge.laobi.icu/badge?page_id=paulshamrat.ColabMDA)


