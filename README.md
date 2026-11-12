# ColabMDA

[![Release](https://img.shields.io/github/v/tag/paulshamrat/ColabMDA)](https://github.com/paulshamrat/ColabMDA/tags) [![Python CI](https://github.com/paulshamrat/ColabMDA/actions/workflows/python-test.yml/badge.svg?branch=main)](https://github.com/paulshamrat/ColabMDA/actions/workflows/python-test.yml) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![GitHub](https://img.shields.io/badge/GitHub-Repository-blue?logo=github)](https://github.com/paulshamrat/ColabMDA) [![Issues](https://img.shields.io/github/issues/paulshamrat/ColabMDA)](https://github.com/paulshamrat/ColabMDA/issues) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black) [![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff) [![GitHub Stars](https://img.shields.io/github/stars/paulshamrat/ColabMDA)](https://github.com/paulshamrat/ColabMDA/stargazers)

**ColabMDA** is a specialized framework for running high-quality Molecular Dynamics simulations on Google Colab without the fear of losing your progress. The biggest challenge with Google Colab is session timeouts and disconnections, which frequently destroy hours of simulation work. **ColabMDA fixes this** with a "resume-safe" architecture that continuously saves checkpoints directly to Google Drive. If your session expires, you can resume exactly where you left off with one simple command. From homology modeling and point mutations to 24-stage equilibration and publication-ready comparative figures, ColabMDA handles the entire pipeline using free cloud GPUs.

---

## 1. Getting started

These tutorials guide you from protein structure preparation to molecular dynamics
simulation and analysis. Start with Trp-cage for a quick test, or use KRAS to compare
the WT and G12D variants. Each tutorial includes installation and Google login
instructions. Use Linux or
Windows through [WSL2](https://learn.microsoft.com/en-us/windows/wsl/install) for
the local commands. Developer setup and support are in [CONTRIBUTING.md](CONTRIBUTING.md).

| # | Tutorial | Preparation | Simulation | Purpose |
| :---: | :--- | :--- | :--- | :--- |
| 01 | [Trp-cage apo](docs/tutorials/colab/trpcage_apo.md) | Colab | Colab | Quick test |
| 02 | [KRAS apo](docs/tutorials/colab/kras_apo_colab.md) | Colab | Colab | Cloud workflow |
| 03 | [KRAS apo](docs/tutorials/colab/kras_apo.md) | Local | Colab | Local preparation, cloud simulation |

## 2. Workflow overview

Fetch experimental coordinates and the target sequence, build the required models,
prepare the solvated system, equilibrate, run checkpointed production, then merge
and analyze the trajectories. Each tutorial shows the files produced after each
stage and ends with a collapsed directory tree.

Local preparation can be transferred to Colab for GPU simulation. Results and
checkpoints stay on Google Drive; analyzed deliverables can be pulled back for
local inspection. The simulation engine also supports workstation and HPC use.

## 3. Project directory structure

<details>
<summary><b>📁 View Standard Project Directory Layout</b></summary>

```text
research/projects/kras-example/          # Working project, separate from the tool repository
└── data/
    └── apo/                            # Protein in water
        ├── structure/P01116_4LDJ/       # Experimental inputs and WT/G12D models
        │   ├── raw/                    # Archive, template and sequence
        │   ├── variants/               # WT and G12D model folders
        │   └── STRUCTURE_DONE          # Completed fetch and mapping
        ├── simulation/                 # One folder per variant
        │   ├── WT/                     # prep/, eq/, md/r1/, md/r2/, md/r3/
        │   └── G12D/                   # Same simulation stages
        └── analysis/                   # QC, replica metrics, summaries and figures
            ├── WT/                     # eq/, r1–r3/, aggregate/
            ├── G12D/                   # Same analysis stages
            └── compare/                # WT/G12D comparisons
```
</details>

<details>
<summary><b>📜 Changelog History (v0.3.0, v0.2.0, v0.1.0)</b></summary>

### v0.3.0 (Native Biophysical Architecture & 24-Stage Equilibration)
*   **100% Native Domain Outputs:** Zero-JSON and Zero-TSV scientific output tree (`data/`) using strictly native biophysical formats (`.pdb`, `.prmtop`, `.inpcrd`, `.dcd`, `.xml`, `.log`, `.csv`, `.png`).
*   **Unified CLI (`cmda`):** Single-keyword scientific verbs (`fetch`, `model`, `prep`, `stage`, `sim`, `status`, `merge`, `analyze`, `compare`, `push`, `pull`, `view`, `snapshots`).
*   **24-Stage Staged Equilibration Ladder:** Progressive NVT heating ($50\text{ K} \to 300\text{ K}$) and NPT restraint release ladder with automated headless QC plotting (`qc.png`).
*   **Automated Solvation & Titration:** Amber `ff19SB` parameters, OPC 4-point water box, and physiological PROPKA 3 pKa titration at pH 7.4.
*   **Core Workflows:** Verified 5-minute Trp-cage cloud quickstart and KRAS Apo (WT vs. G12D) workflows for both Google Colab and local workstations.

### v0.2.0 (Kinase Architectures, Amber ff19SB & Resilience)
*   **Universal Setup Script (`setup.sh`):** Single bootstrap script automatically configuring both local Linux workstations and Google Colab cloud instances.
*   **Local Preparation Shift:** Moving `cmda prep` (PROPKA pH 7.4 + Amber `tleap` solvation) locally into Phase 1, enabling local PyMOL visual inspection of the water box and keeping Google Colab 100% unattended.
*   **Separation of Concerns:** Clean separation of `cmda fetch` (pure structural retrieval and mapping) from `cmda model` (coordinate building and mutations).
*   **One-Command Sync (`cmda push` & `cmda pull`):** Direct streaming between local workstation and remote Google Drive via Google Colab CLI, eliminating manual archiving and rclone setup.
*   **Cofactor Donor Isolation:** Dedicated `donors/<template>/` directory hierarchy preserving donor crystal integrity.
*   **Target Hierarchy (`<uniprot>_<pdb>`):** Clean segregation for multi-template and multi-protein workspaces.
*   **Resume-Safe Checkpointing:** Automated recovery from Google Colab session timeouts.

### v0.1.0 (Initial Beta)
*   Modular pipeline for EM, NVT, NPT, and Production MD.
*   Checkpointing logic for fail-safe simulations on Google Colab.
*   Wild-Type building and mutation support via Modeller.
</details>

---

## Acknowledgements
- [OpenMM](https://openmm.org) & [PDBFixer](https://github.com/openmm/pdbfixer)
- [Modeller](https://salilab.org/modeller/)
- [AmberTools](https://ambermd.org/AmberTools.php)
- [MDAnalysis](https://www.mdanalysis.org) & [MDTraj](https://www.mdtraj.org)
- [PDB2PQR & PROPKA](https://server.poissonboltzmann.org/)
- [NumPy](https://numpy.org), [Matplotlib](https://matplotlib.org), [Biopython](https://biopython.org)
- [Google Colab](https://colab.research.google.com) & [Miniforge/Conda](https://github.com/conda-forge/miniforge)

---

## Citation

If you use ColabMDA in your research, please cite the underlying studies:

> Paul SK, Saddam M, Rahaman KA, Choi JG, Lee SS, Hasan M. **Molecular modeling, molecular dynamics simulation, and essential dynamics analysis of grancalcin: An upregulated biomarker in experimental autoimmune encephalomyelitis mice.** *Heliyon*. 2022 Oct 23;8(10):e11232. doi: [10.1016/j.heliyon.2022.e11232](https://doi.org/10.1016/j.heliyon.2022.e11232). PMID: 36340004; PMCID: PMC9626934.
>
> [![Cited by](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fapi.openalex.org%2Fworks%2Fhttps%3A%2F%2Fdoi.org%2F10.1016%2Fj.heliyon.2022.e11232&query=%24.cited_by_count&label=cited%20by&color=blue)](https://doi.org/10.1016/j.heliyon.2022.e11232)
>
> Paul SK, Saddam M, Tabassum N, Hasan M. **Molecular dynamics simulation of wild and mutant proteasome subunit beta type 8 (PSMB8) protein: Implications for restoration of inflammation in experimental autoimmune encephalomyelitis pathogenesis.** *Heliyon*. 2024 Dec 15;11(1):e41166. doi: [10.1016/j.heliyon.2024.e41166](https://doi.org/10.1016/j.heliyon.2024.e41166). PMID: 39802026; PMCID: PMC11719297.
>
> [![Cited by](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fapi.openalex.org%2Fworks%2Fhttps%3A%2F%2Fdoi.org%2F10.1016%2Fj.heliyon.2024.e41166&query=%24.cited_by_count&label=cited%20by&color=blue)](https://doi.org/10.1016/j.heliyon.2024.e41166)

---

<p align="center">
  <img src="https://komarev.com/ghpvc/?username=paulshamrat&repo=ColabMDA&label=Visitors&color=blue&style=flat-square" alt="Visitors" />
</p>
