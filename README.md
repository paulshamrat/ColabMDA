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
| **Dependencies** | `OpenMM`, `Modeller`, `MDAnalysis`, `MDTraj`, `OpenFF (Optional)`, `OpenMMForceFields (Optional)` |
| **Platform** | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/paulshamrat/ColabMDA/blob/main/notebooks/05-colabmd-simulation-2024.ipynb) `Linux` `HPC (SLURM)` |
| **Structure** | [Source Layout](#-project-structure) |

### 📂 Project Structure
ColabMDA is organized into clear, functional modules:
*   `src/colabmda/modeller/`: Homology modeling engine (Biological numbering supported).
*   `src/colabmda/openmm/`: OpenMM simulation engines (Modular EM/NVT/NPT/MD).
*   `envs/`: Automated installation scripts for scientific environments.
*   `scripts/`: Quick-start bootstrap scripts for Google Colab.
*   `notebooks/`: Ready-to-use Colab notebooks for simulation and analysis.

---

## 📖 User Manual & Simulation Guides

To keep our workspace clean and maintainable, the detailed manuals have been separated into dedicated guides:

### 🛠️ [Installation & Setup Guide](docs/installation.md)
Detailed walkthroughs for Google Colab environment setup, license keys, local workstations setup, and SLURM configurations.

### 💧 [Protein-in-Water Simulation Guide](docs/protein_water_simulation.md)
Step-by-step instructions for staging protein systems in water, running resume-safe production MD, trajectory merging, stability gating, and PyMOL visualization.

### 🧬 [Protein-Ligand Simulation Guide](docs/protein_ligand_simulation.md)
Instructions for modeling small molecule ligands using dynamic GAFF2/AM1-BCC parameterization and retaining structural Magnesium ($Mg^{2+}$) cofactors.

---

## 📋 Recommended Directory Strategy

We recommend organizing your research project directories as follows:

```text
/content/drive/MyDrive/ColabMDA/
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

---

## Acknowledgements
- [OpenMM](https://openmm.org) & [PDBFixer](https://github.com/openmm/pdbfixer)
- [Modeller](https://salilab.org/modeller/)
- [MDAnalysis](https://www.mdanalysis.org) & [MDTraj](https://www.mdtraj.org)
- [NumPy](https://numpy.org), [Matplotlib](https://matplotlib.org), [Biopython](https://biopython.org)
- [Google Colab](https://colab.research.google.com/assets/colab-badge.svg) & [Miniforge/Conda](https://github.com/conda-forge/miniforge)

---

## Citation

Initially, this repository was based on a GROMACS-on-Colab workflow, but we transformed it to utilize OpenMM to improve feasibility, performance, and stability for automated simulation recovery on cloud platforms.

If you use this tool, please consider citing the underlying studies:

> Paul SK, Saddam M, Rahaman KA, Choi JG, Lee SS, Hasan M. **Molecular modeling, molecular dynamics simulation, and essential dynamics analysis of grancalcin: An upregulated biomarker in experimental autoimmune encephalomyelitis mice.** *Heliyon*. 2022 Oct 23;8(10):e11232. doi: [10.1016/j.heliyon.2022.e11232](https://doi.org/10.1016/j.heliyon.2022.e11232). PMID: 36340004; PMCID: PMC9626934.
>
> [![Cited by](https://img.shields.io/endpoint?url=https%3A%2F%2Fapi.juleskreuer.eu%2Fcitation-badge.php%3Fshield%26doi%3D10.1016%2Fj.heliyon.2022.e11232&label=cited%20by)](https://juleskreuer.eu/projects/citation-badge)
>
> Paul SK, Saddam M, Tabassum N, Hasan M. **Molecular dynamics simulation of wild and mutant proteasome subunit beta type 8 (PSMB8) protein: Implications for restoration of inflammation in experimental autoimmune encephalomyelitis pathogenesis.** *Heliyon*. 2024 Dec 15;11(1):e41166. doi: [10.1016/j.heliyon.2024.e41166](https://doi.org/10.1016/j.heliyon.2024.e41166). PMID: 39802026; PMCID: PMC11719297.
>
> [![Cited by](https://img.shields.io/endpoint?url=https%3A%2F%2Fapi.juleskreuer.eu%2Fcitation-badge.php%3Fshield%26doi%3D10.1016%2Fj.heliyon.2024.e41166&label=cited%20by)](https://juleskreuer.eu/projects/citation-badge)

---

<p align="center">
  <img src="https://komarev.com/ghpvc/?username=paulshamrat&repo=ColabMDA&label=Visitors&color=blue&style=flat" alt="Visitors" />
</p>
