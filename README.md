# ColabMDA

[![Release](https://img.shields.io/github/v/tag/paulshamrat/ColabMDA)](https://github.com/paulshamrat/ColabMDA/tags) [![Python CI](https://github.com/paulshamrat/ColabMDA/actions/workflows/python-test.yml/badge.svg?branch=main)](https://github.com/paulshamrat/ColabMDA/actions/workflows/python-test.yml) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![GitHub](https://img.shields.io/badge/GitHub-Repository-blue?logo=github)](https://github.com/paulshamrat/ColabMDA) [![Issues](https://img.shields.io/github/issues/paulshamrat/ColabMDA)](https://github.com/paulshamrat/ColabMDA/issues) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black) [![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff) [![GitHub Stars](https://img.shields.io/github/stars/paulshamrat/ColabMDA)](https://github.com/paulshamrat/ColabMDA/stargazers)

**ColabMDA** is a specialized framework for running high-quality Molecular Dynamics simulations on Google Colab without the fear of losing your progress. The biggest challenge with Google Colab is session timeouts and disconnections, which frequently destroy hours of simulation work. **ColabMDA fixes this** with a "resume-safe" architecture that continuously saves checkpoints directly to Google Drive. If your session expires, you can resume exactly where you left off with one simple command. From homology modeling and point mutations to 24-stage equilibration and publication-ready comparative figures, ColabMDA handles the entire pipeline using free cloud GPUs.

---

## 1. Quick Installation

Run this one-time setup on your local machine (Linux, macOS, WSL2):

```bash
# 1. Install cmda, OpenMM, AmberTools, PyMOL, and Modeller:
bash -c "$(curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/dev1/scripts/setup.sh)"

# 2. Activate environment:
source ~/miniforge3/etc/profile.d/conda.sh
conda activate cmda

# 3. Log in to Google Colab (one-time browser login):
colab sessions
```

> [!TIP]
> **Modeller Key:** Free for academic research ([register here](https://salilab.org/modeller/registration.html)). Enter it when prompted, or run:  
> `KEY_MODELLER="YOUR_KEY" bash -c "$(curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/dev1/scripts/setup.sh)"`

> [!NOTE]
> **Cloud Setup:** Google Colab GPU instances are configured during Phase 3. See [`docs/installation.md`](docs/installation.md).

<details>
<summary><b>🛠️ Install from Source</b></summary>

```bash
git clone --branch dev1 https://github.com/paulshamrat/ColabMDA.git ~/ColabMDA
cd ~/ColabMDA
bash scripts/setup.sh
```
</details>

---

## 2. Documentation & Standalone Guides

All in-depth scientific tutorials, protonation tables, parameters, and platform deployment guides are modularized in `docs/`:

| # | Guide | Description | Link |
| :-: | :--- | :--- | :--- |
| 01 | **Installation** | Local workstation & Google Colab GPU setup | [`docs/installation.md`](docs/installation.md) |
| 02 | **KRAS Apo** | Protein-only GTPase tutorial (WT vs G12D) | [`docs/tutorials/kras_colab.md`](docs/tutorials/kras_colab.md) |
| 03 | **CDKL5 Apo** | Full kinase in water (WT vs L119R) | [`docs/tutorials/apo.md`](docs/tutorials/apo.md) |
| 04 | **CDKL5 Holo** | Kinase with ATP + $\text{Mg}^{2+}$ cofactors | [`docs/tutorials/atpmg.md`](docs/tutorials/atpmg.md) |
| 05 | **Google Colab** | Cloud GPU provisioning & Drive sync | [`docs/platforms/colab.md`](docs/platforms/colab.md) |
| 06 | **HPC & SLURM** | Cluster submission & batch arrays | [`docs/platforms/hpc.md`](docs/platforms/hpc.md) |

---

<details>
<summary><b>📁 View Standard Project Directory Layout</b></summary>

```text
research/projects/cmda-example/           # Project root in Google Drive / local workspace
└── data/
    ├── apo/                              # Track 1: Apo state (protein in water)
    │   ├── structure/
    │   │   └── O76039_4BGQ/              # Biological Target (UniProt) + Structural Template (PDB)
    │   │       ├── raw/                  # Untouched experimental inputs (.cif, template.pdb, uniprot.fasta)
    │   │       ├── variants/             # Modeled coordinates (WT/WT.pdb, L119R/L119R.pdb)
    │   │       └── STRUCTURE_DONE        # Verification completion marker
    │   ├── simulation/
    │   │   ├── WT/
    │   │   │   ├── 01_prep/              # Solvated OpenMM system (model.pdb, system.xml, integrator.xml)
    │   │   │   ├── 02_eq/                # 24-stage equilibration hub (eq.pdb, state.xml, stages/)
    │   │   │   └── 03_md/                # Production MD replicas (md.dcd, md.log, state.xml, chunks/)
    │   │   └── L119R/
    │   │       ├── 01_prep/
    │   │       ├── 02_eq/
    │   │       └── 03_md/{r1, r2, r3}/
    │   └── analysis/                     # Mirrored reports (<variant>/02_eq/qc.png, qc.txt) & flat compare/
    │
    └── atpmg/                            # Track 2: Holo state (Kinase + ATP + Mg²⁺)
        ├── structure/
        │   ├── O76039_4BGQ/              # Inherited from apo (identical starting coordinates)
        │   └── cofactors/atpmg/          # Homologous donor (CDK2 8FP5: atp.pdb, mg.pdb)
        ├── simulation/
        │   └── WT/01_prep/               # tleap combine { REC LIG MG } + Hu 2024 parameters
        └── analysis/
```
</details>

<details>
<summary><b>📜 Changelog History (v0.2.0, v0.1.0)</b></summary>

### v0.2.0 (Kinase Architectures, Amber ff19SB & Resilience)
*   **Universal Setup Script (`setup.sh`):** Single bootstrap script automatically configuring both local Linux workstations and Google Colab cloud instances.
*   **Local Preparation Shift:** Moving `cmda prep` (PROPKA pH 7.4 + Amber `tleap` solvation) locally into Phase 1, enabling local PyMOL visual inspection of the water box and keeping Google Colab 100% unattended.
*   **Separation of Concerns:** Clean separation of `cmda fetch` (pure structural retrieval and mapping) from `cmda model` (coordinate building and mutations).
*   **One-Command Sync (`cmda push` & `cmda pull`):** Direct streaming between local workstation and remote Google Drive via Google Colab CLI, eliminating manual archiving and rclone setup.
*   **Cofactor Donor Isolation:** Dedicated `donors/<template>/` directory hierarchy preserving donor crystal integrity.
*   **Target Hierarchy (`<uniprot>_<pdb>`):** Clean segregation for multi-template and multi-protein workspaces.
*   **Unified CLI (`cmda`):** Intuitive top-level commands (`fetch`, `model`, `assemble`, `prep`, `push`, `stage`, `sim`, `status`, `merge`, `analyze`, `compare`, `pull`, `view`, `snapshots`).
*   **Dual State Architecture:** Clean separation between unliganded **`apo`** dynamics and catalytic **`atpmg`** holo dynamics.
*   **Automated Triplicate Comparison:** Single-call `cmda compare --state <state> --variants WT,MUT` generating publication-grade Mean ± SD overlay plots.
*   **24-Stage Equilibration Hub:** Robust heating and restraint release ladder guaranteeing stable density and temperature.
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
