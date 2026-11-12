# ColabMDA

[![Release](https://img.shields.io/github/v/tag/paulshamrat/ColabMDA)](https://github.com/paulshamrat/ColabMDA/tags) [![Python CI](https://github.com/paulshamrat/ColabMDA/actions/workflows/python-test.yml/badge.svg?branch=main)](https://github.com/paulshamrat/ColabMDA/actions/workflows/python-test.yml) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![GitHub](https://img.shields.io/badge/GitHub-Repository-blue?logo=github)](https://github.com/paulshamrat/ColabMDA) [![Issues](https://img.shields.io/github/issues/paulshamrat/ColabMDA)](https://github.com/paulshamrat/ColabMDA/issues) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black) [![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff) [![GitHub Stars](https://img.shields.io/github/stars/paulshamrat/ColabMDA)](https://github.com/paulshamrat/ColabMDA/stargazers)

**ColabMDA** is a specialized framework for running high-quality Molecular Dynamics simulations on Google Colab without the fear of losing your progress. The biggest challenge with Google Colab is session timeouts and disconnections, which frequently destroy hours of simulation work. **ColabMDA fixes this** with a "resume-safe" architecture that continuously saves checkpoints directly to Google Drive. If your session expires, you can resume exactly where you left off with one simple command. From homology modeling and point mutations to 24-stage equilibration and publication-ready comparative figures, ColabMDA handles the entire pipeline using free cloud GPUs.

---

## 1. Installation

ColabMDA divides work between your local computer and Google Colab:

### Prerequisites

On a fresh Ubuntu or WSL2 installation, install these system tools first:

```bash
sudo apt-get update
sudo apt-get install -y git curl wget ca-certificates
```

Git is needed when pip installs ColabMDA from its GitHub repository. The installer
checks that `git`, `curl`, and `wget` are on `PATH` before downloading Miniforge,
changing environments, or prompting for a MODELLER key. These checks also apply
when running the installer on Colab.

### Local Workstation
Run this one-time setup in your terminal on **Linux (Ubuntu)** or Windows via **[WSL2](https://learn.microsoft.com/en-us/windows/wsl/install)**:

```bash
# 1. Install cmda, OpenMM, AmberTools, PyMOL, and Modeller:
bash -c "$(curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/main/scripts/setup.sh)"

# 2. Activate environment:
source ~/miniforge3/etc/profile.d/conda.sh
conda activate cmda

# 3. Log in to Google Colab (one-time browser login):
colab sessions
```

> [!TIP]
> **Modeller Key:** Free for academic research ([register here](https://salilab.org/modeller/registration.html)). Enter it when prompted, or run:  
> `KEY_MODELLER="YOUR_KEY" bash -c "$(curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/main/scripts/setup.sh)"`

### Google Colab (Cloud)
Run this once inside your Colab terminal (GPU or CPU runtime) in ~90 seconds:

```bash
# 1. Install cmda, OpenMM (CUDA/CPU), and AmberTools:
curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/main/scripts/setup.sh | bash

# 2. Activate environment:
source /root/miniforge3/etc/profile.d/conda.sh
conda activate cmda
```

> [!TIP]
> **GPU vs. CPU Runtime:** Use a GPU runtime for running simulations (`cmda sim`). For pulling deliverables (`cmda pull`) or running quick checks, you can launch an instant, free CPU session without waiting for GPU quota.

<details>
<summary><b>🛠️ Install from Source & Run Tests</b></summary>

```bash
git clone https://github.com/paulshamrat/ColabMDA.git ~/ColabMDA
cd ~/ColabMDA
bash scripts/setup.sh

# Run test suite:
pytest tests/ -v
```
</details>

---

## 2. Tutorials

All in-depth scientific tutorials, protonation tables, and simulation parameters are modularized in `docs/tutorials/colab/`:

| # | Level | Tutorial | Platform | Description | Link |
| :-: | :--- | :--- | :--- | :--- | :--- |
| 01 | Quickstart | **Trp-Cage: Protein in Water** | Colab | 5-minute GPU simulation of a 20-residue miniprotein | [`trpcage_apo.md`](docs/tutorials/colab/trpcage_apo.md) |
| 02 | Beginner | **KRAS: Protein in Water** | Colab | Direct cloud simulation and analysis on Google Drive | [`kras_apo_colab.md`](docs/tutorials/colab/kras_apo_colab.md) |
| 03 | Beginner | **KRAS: Protein in Water** | Local + Colab | Local modeling and preparation with cloud GPU simulation | [`kras_apo.md`](docs/tutorials/colab/kras_apo.md) |

---

<details>
<summary><b>📁 View Standard Project Directory Layout</b></summary>

```text
research/projects/cmda-example/           # Project root in Google Drive / local workspace
└── data/
    ├── apo/                              # Track 1: Apo state (protein in water)
    │   ├── structure/
    │   │   └── O76039_4BGQ/                      # Biological Target (UniProt) + Structural Template (PDB)
    │   │       ├── raw/                          # Untouched experimental inputs (cmda fetch)
    │   │       │   ├── 4bgq.cif                  # Raw crystallographic archive record
    │   │       │   ├── template.pdb              # Cleaned monomer (waters/buffers stripped)
    │   │       │   └── uniprot.fasta             # Canonical sequence for biological numbering
    │   │       │
    │   │       ├── variants/                     # Modeled coordinates (cmda model)
    │   │       │   ├── WT/
    │   │       │   │   ├── WT.pdb                # Canonical modeled Wild-Type coordinate
    │   │       │   │   ├── target.B99990001.pdb  # Modeller build
    │   │       │   │   ├── target.B99990001_with_cryst.pdb # Model with crystal headers
    │   │       │   │   ├── target.ini            # Starting coordinate mesh
    │   │       │   │   ├── target.rsr            # Spatial restraints
    │   │       │   │   ├── target.sch            # Optimization schedule
    │   │       │   │   ├── target.D00000001      # Restraint progress log
    │   │       │   │   ├── target.V99990001      # Restraint violation profile
    │   │       │   │   ├── target.ali            # PIR target sequence alignment
    │   │       │   │   ├── alignment.ali         # PIR template pairwise alignment
    │   │       │   │   ├── O76039.fasta          # Sequence reference
    │   │       │   │   ├── 4bgq_clean.pdb        # Sanitized template
    │   │       │   │   ├── 4bgq_orig.pdb         # Raw template backup
    │   │       │   │   ├── 4bgq.pdb              # Isolated chain coordinate
    │   │       │   │   ├── pipeline.log          # Modeller build log
    │   │       │   │   └── MODEL_DONE            # Modeling completion gate
    │   │       │   │
    │   │       │   └── L119R/                    # Kinase-inactivating point mutant
    │   │       │       ├── L119R.pdb             # Mutated coordinate
    │   │       │       ├── WT_L119R.pdb          # Direct mutation coordinate
    │   │       │       ├── WT_L119R_with_cryst.pdb # Mutant with crystal headers
    │   │       │       ├── mutate.log            # Rotamer sampling & energy log
    │   │       │       └── MODEL_DONE            # Mutant completion gate
    │   │       │
    │   │       └── STRUCTURE_DONE                # Verification completion marker
    │   ├── simulation/
    │   │   ├── WT/
    │   │   │   ├── prep/                 # Solvated Amber system
    │   │   │   │   ├── model.pdb         # Input receptor coordinates
    │   │   │   │   ├── solvated.pdb      # Solvated visualization PDB
    │   │   │   │   ├── system.prmtop     # Amber ff19SB parameters
    │   │   │   │   ├── system.inpcrd     # Initial orthogonal coordinates
    │   │   │   │   └── system.xml        # Serialized OpenMM system
    │   │   │   │
    │   │   │   ├── eq/                   # 24-stage equilibration
    │   │   │   │   ├── eq.pdb            # Final relaxed coordinates
    │   │   │   │   ├── system.xml        # OpenMM system topology & parameters
    │   │   │   │   ├── state.xml         # OpenMM thermodynamic state snapshot
    │   │   │   │   └── stages/           # Granular restart checkpoints (state_01.xml ... state_24.xml)
    │   │   │   │
    │   │   │   └── md/                   # Production MD replicas
    │   │   │       ├── r1/               # Production replica 1
    │   │   │       │   ├── md.dcd        # Continuous merged trajectory
    │   │   │       │   ├── md.log        # Continuous merged thermodynamic log
    │   │   │       │   ├── state.xml     # Latest OpenMM simulation state snapshot
    │   │   │       │   ├── prod.chk      # OpenMM checkpoint binary for resume
    │   │   │       │   ├── stripped.prmtop      # Solvent-stripped topology for PyMOL & analysis
    │   │   │       │   ├── prod_keep_ligand.dcd # Continuous PBC-autoimaged stripped trajectory
    │   │   │       │   ├── ca.prmtop     # Cα backbone topology
    │   │   │       │   ├── prod_ca.dcd   # Cα backbone trajectory
    │   │   │       │   ├── MERGE_DONE    # Merge completion marker
    │   │   │       │   └── chunks/       # Raw crash-recovery chunks (0001.dcd, 0001.log, ...)
    │   │   │       ├── r2/               # Independent velocity seed
    │   │   │       └── r3/               # Triplicate simulation
    │   │   └── L119R/
    │   │       ├── prep/
    │   │       ├── eq/
    │   │       └── md/{r1, r2, r3}/
    │   └── analysis/                     # Mirrored deliverables (Zero JSON policy)
    │       ├── WT/
    │       │   ├── eq/                   # qc.png, qc.txt
    │       │   └── md/                   # r1/ (rmsd.png, rg.png, rmsf.png, summary.txt), r2/, r3/
    │       ├── L119R/
    │       │   ├── eq/
    │       │   └── md/{r1, r2, r3}/
    │       └── compare/                  # Flat cross-variant comparison (rmsd.png, rg.png, rmsf.png, summary.txt)
    │
    └── atpmg/                            # Track 2: Holo state (Kinase + ATP + Mg²⁺)
        ├── structure/
        │   ├── O76039_4BGQ/              # Inherited from apo (identical starting coordinates)
        │   └── cofactors/atpmg/          # Homologous donor (CDK2 8FP5: atp.pdb, mg.pdb)
        ├── simulation/
        │   └── WT/prep/                  # tleap combine { REC LIG MG } + Hu 2024 parameters
        └── analysis/
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
