# ColabMDA

**ColabMDA** is a specialized tool that lets you run high-quality Molecular Dynamics simulations on Google Colab without the fear of losing your work. The biggest problem with Colab is that it disconnects, often destroying hours of simulation data. **ColabMDA fixes this** with a "resume-safe" system that automatically saves your progress to Google Drive. If your session expires, you can resume exactly where you left off with one simple command. From modeling protein mutations to generating publication-ready analysis, ColabMDA handles the complex setup for you, making it the easiest way to get high-quality MD results using free cloud GPUs.

📖 **Full Documentation:** Visit our official manual at [colabmda.readthedocs.io](https://colabmda.readthedocs.io/)

---

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
| **Dependencies** | `OpenMM`, `AmberTools (tleap, cpptraj)`, `Modeller`, `MDAnalysis`, `MDTraj`, `PDB2PQR / PROPKA`, `OpenFF (Optional)` |
| **Platform** | `Linux` · `Google Colab (CLI / Terminal)` · `HPC (SLURM)` |
| **Structure** | [Source Layout](#-project-structure) |

### 📂 Project Structure
ColabMDA is organized into clear, functional modules:
*   `src/colabmda/cli.py`: Central unified CLI (`colabmda` and shortcut `cmda`).
*   `src/colabmda/core/`: System taxonomy (`SystemSpec`), project layout (`ProjectPaths`), and divalent ion parameters.
*   `src/colabmda/amber/`: AmberTools parameterization (`tleap`, Hu 2024 quantum nucleotide parameters, `cpptraj` autoimaging & solvent stripping).
*   `src/colabmda/modeller/`: Homology modeling engine (UniProt biological numbering, mutagenesis, geometric substrate docking).
*   `src/colabmda/openmm/`: OpenMM simulation engines (Modular prep/EM/NVT/NPT/QC/MD/analysis).
*   `src/colabmda/remote/`: Google Colab GPU remote session management via `colab` CLI.
*   `scripts/`: Single standalone bootstrap script (`scripts/setup.sh`).
*   `tests/`: Unit and regression test suite (`pytest`).
*   `docs/`: Sphinx and ReadTheDocs scientific documentation source.

---

## 1. Quick Start (Google Colab Remote Execution)

### 1.1. Local Terminal Runner (Recommended — Zero Browser)
Drive remote Colab GPU simulations directly from your local terminal using the `colab` CLI:

```bash
# 1. Provision a remote GPU instance and mount Google Drive
colab new --gpu T4
colab drivemount /content/drive

# 2. Bootstrap OpenMM, AmberTools, and ColabMDA on the VM (one-time)
cmda colab setup

# 3. Connect to your remote Colab GPU shell
colab console
```

> [!TIP]
> **Headless Execution:** If you prefer running simulations directly from your local terminal without entering `colab console`, you can dispatch any command with `cmda colab exec "<command>"`.

> [!IMPORTANT]
> **Modeller License (Required for Homology Modeling & Mutagenesis):**
> A Modeller license key is **required** to run homology modeling and point-mutation workflows (`cmda modeller`).
> 🔑 **Get a Free Academic License:** Register at [salilab.org/modeller/registration.html](https://salilab.org/modeller/registration.html) (Licenses are free for academic use and sent instantly via email).
> To provision MODELLER during bootstrap:
> ```bash
> KEY_MODELLER="YOUR_LICENSE_KEY" WITH_MODELLER=1 cmda colab setup
> ```

<details>
<summary><b>🌐 Alternative: Running inside the Colab Browser UI</b></summary>

If you prefer using the Google Colab web browser interface:
1. **Enable GPU:** Go to `Runtime` -> `Change runtime type` -> select **T4 GPU**.
2. **Mount Drive:** Run in a notebook cell:
   ```python
   from google.colab import drive
   drive.mount('/content/drive')
   ```
3. **Install ColabMDA:** Open Colab Terminal (**⋮** menu -> **Terminal**) and run:
   ```bash
   cd /content
   curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/main/scripts/setup.sh -o setup.sh
   bash setup.sh latest
   ```
</details>

---

### 💡 Tip: Timeout Immunity & Seamless Resumption
* **Timeout-Proof:** Running via `cmda colab exec` keeps the remote kernel in a `BUSY` execution state, permanently preventing Colab's 30-minute idle VM reaper from shutting down your run.
* **Granular Recovery:** If a Colab session ever expires at its 12-hour platform limit, simply restart the command. ColabMDA automatically detects saved `stage_*.state.xml` checkpoints on Google Drive and fast-forwards directly to the unfinished step.

---

## 2. Workstation & HPC Installation

<details>
<summary><b>💻 A. Local Workstation Setup (Conda / Mamba)</b></summary>

ColabMDA runs natively on Linux workstations with NVIDIA GPUs. Use the provided `environment.yml` to provision the production environment:
```bash
git clone https://github.com/paulshamrat/ColabMDA.git
cd ColabMDA
mamba env create -f environment.yml
conda activate colabmda
```
</details>

<details>
<summary><b>🧪 B. Running the Test Suite</b></summary>

To verify force fields, peptide docking, Amber leap parameterization, and CLI execution:
```bash
pytest tests/ -v
```
</details>

<details>
<summary><b>🏢 C. HPC Usage (SLURM)</b></summary>

ColabMDA generates decoupled SLURM batch scripts automatically under `<project>/scripts/submit.slurm`. Since trajectories are chunked with automatic checkpointing, jobs run reliably within cluster walltime limits.
</details>

---

## 3. Recommended Directory Strategy & Architectural Principles

ColabMDA enforces strict architectural separation to ensure reproducibility across local workstations, cloud GPUs, and HPC clusters:

1. **Zero Data in `src/`:** Structural inputs (`*.pdb`), topologies (`*.prmtop`), and trajectories belong strictly in the user's project workspace. Non-standard residues like phosphotyrosine (`pY171` -> `PTR171`) are constructed mathematically with zero hardcoded PDBs in source code.
2. **Decoupled Execution Scripts:** Execution shell scripts and SLURM batch scripts reside exclusively in `<project>/scripts/`, never inside simulation data folders.
3. **Single Bootstrap Script:** Standalone environment provisioning is handled by a single script: `scripts/setup.sh`.

```text
<project_root>/                               # E.g. test-project (local) or /content/drive/MyDrive/projects/test-project (Colab)
├── scripts/                                  # Decoupled execution & batch submission scripts
│   ├── run.sh                                # Standalone runner for local workstation or Colab
│   ├── submit_eq.slurm                       # HPC SLURM equilibration batch script
│   ├── submit_prod.slurm                     # HPC SLURM production array job script
│   └── submit_chain.sh                       # Chained submission helper (EQ -> Array Prod)
│
└── data/                                     # Single root for all scientific domain data
    ├── structure/                            # Raw and modeled structural inputs (Zero raw PDBs in src/)
    │   ├── 4bgq/                             # Target kinase structure (PDB: 4BGQ)
    │   │   ├── wt/                           # Wild-Type modeled receptor (4bgq_wt.pdb, prep/4bgq_wt_cleaned.pdb)
    │   │   └── mutants/4bgq_l119r/           # Pathogenic hinge mutant (4bgq_l119r.pdb, prep/4bgq_l119r_cleaned.pdb)
    │   ├── template/                         # Structural templates (4BGQ.pdb, 8FP5.pdb)
    │   ├── peptide/                          # Substrate peptides (EB2.pdb)
    │   └── ligand/                           # Small molecules & nucleotides (ATP.sdf, ADP.sdf)
    ├── simulation/                           # Staged & executed simulation data
    │   ├── apo/                              # Track 1: Unliganded kinase state
    │   │   ├── WT/                           # Wild-Type CDKL5
    │   │   │   ├── prep/                     # Parameterization & Amber topology (model.pdb, solvated.pdb, system.prmtop, prep.json)
    │   │   │   ├── eq/                       # Shared 24-stage equilibration hub (eq.pdb, qc.png, qc.json)
    │   │   │   ├── r1/                       # Production replica 1 (prod.dcd, prod.log, restart.chk)
    │   │   │   ├── r2/                       # Production replica 2 (independent velocity seed from eq/)
    │   │   │   └── r3/                       # Production replica 3
    │   │   └── L119R/                        # Pathogenic hinge mutant (same modular structure)
    │   └── atp_mg/                           # Track 2: Pre-catalytic holo state (CDKL5 + ATP + Mg2+)
    │       ├── WT/                           # Holo WT (prep/, eq/, r1/, r2/, r3/)
    │       └── L119R/                        # Holo L119R (prep/, eq/, r1/, r2/, r3/)
    └── analysis/                             # Trajectory metrics, CSVs, and publication plots
        ├── apo/                              # Unliganded trajectory analysis
        │   ├── WT/                           # Per-variant metrics (rmsd.csv, rmsf.csv, rg.csv)
        │   ├── L119R/                        # Per-variant metrics (rmsd.csv, rmsf.csv, rg.csv)
        │   └── compare/                      # Cross-variant comparative plots (wt_vs_l119r.png)
        ├── atp_mg/                           # Holo trajectory analysis
        │   ├── WT/                           # Per-variant metrics (rmsd.csv, rmsf.csv, rg.csv)
        │   ├── L119R/                        # Per-variant metrics (rmsd.csv, rmsf.csv, rg.csv)
        │   └── compare/                      # Cross-variant comparative plots (wt_vs_l119r.png)
        └── compare_apo_vs_holo/              # Conformational response to nucleotide binding (apo vs holo)
```

> 💡 **Why Shared `eq/`?**
> Equilibration (EM, NVT, NPT, QC) only needs to be run once per system. Replicas (`r1`, `r2`, `r3`) fork directly from `eq/npt.state.xml` with distinct random velocity seeds, saving substantial GPU time while maintaining statistical independence.

---

## 4. Unified Biomolecular Simulation Engine (`cmda`)

ColabMDA provides a streamlined, single-keyword CLI accessible via `cmda` (or `colabmda`) that unifies structure preparation, parameterization, OpenMM simulation dispatch, trajectory post-processing, and multi-replica analysis.

### 4.1. Core Capabilities: Advanced Biomolecular Modeling & Dynamics
* **Universal Protein & Mutant Support:** Prepares and simulates arbitrary wild-type and mutant protein structures with automated biological residue renumbering and Amber ff19SB/ff14SB parameterization.
* **Post-Translational Phosphorylation (`pY171`):** Mathematically converts target tyrosine residues to Amber `PTR` with sp3-hybridized oxygens and $-2.0e$ net charge, dynamically routing Amber `phosaa19SB` (or `phosaa14SB`) force fields.
* **Quantum Nucleotide Parameters:** Integrates Hu et al. 2024 quantum-calculated parameters for ATP, ADP, and GDP.
* **Divalent Metal Ion Catalysis:** Supports catalytic dual $Mg^{2+}$ or $Mn^{2+}$ coordination using Amber 12-6-4 / standard LJ parameters.
* **Complex & Substrate Peptide Docking:** Docks substrate peptides (e.g. EB2 `KRPSSAPSG` at Ser222) via continuous geometric backbone vector propagation ($d_{\text{C-N}} = 1.33\text{ Å}$), eliminating atom overlap and NaN minimization errors.
* **Ultra-Fast Solvation & Neutralization:** Uses Amber `tleap` with randomized ion placement (`addIonsRand`), solvating systems in under 2 seconds.
* **Dual-Stream Trajectory Stripping:** Automates `cpptraj` PBC autoimaging and strips bulk water, producing lightweight dry (`keep_ligand`) and hydration-shell trajectories (>90% storage savings).

---

### 4.2. Single-Keyword Command Reference (`cmda`)

| Command | Syntax | Description |
| :--- | :--- | :--- |
| **`modeller`** | `cmda modeller <build\|mutate>` | Builds homology models (UniProt biological numbering) and introduces point mutations. |
| **`stage`** | `cmda stage --state <tag> --variant <var>` | Generates decoupled project trees and runner scripts (`run.sh`, `submit.slurm`). |
| **`prep`** | `cmda prep --state <tag> --variant <var>` | Assembles coordinates and parameterizes Amber topology (`system.prmtop`, `system.inpcrd`). |
| **`sim`** | `cmda sim --state <tag> --variant <var> --stage <stages>` | Dispatches OpenMM simulations (`em`, `nvt`, `npt`, `qc`, `md`, or `eq`). |
| **`merge`** | `cmda merge --state <tag> --variant <var>` | Wraps PBC and strips bulk solvent via `cpptraj` (`stripped.system.prmtop`). |
| **`analyze`** | `cmda analyze --state <tag> --variant <var>` | Computes RMSD, $R_g$, RMSF, PCA, and contact maps. |
| **`compare`** | `cmda compare --state <tag> --variants <v1,v2>` | Renders multi-replica or cross-variant comparative overlay plots. |
| **`colab`** | `cmda colab <new\|setup\|exec\|stop>` | Provisions remote GPU VMs, bootstraps environments, and dispatches timeout-proof simulations. |

---

### 4.3. Foundation: Homology Modeling & Biological Mutagenesis (`cmda modeller`)

ColabMDA structures simulation campaigns into clear, progressive biophysical tracks—from unliganded baseline dynamics to fully active signaling complexes. Here we demonstrate the full pipeline using the human **CDKL5 catalytic kinase domain (UniProt `O76039`, residues 1–303, template PDB `4BGQ`)** comparing the **Wild-Type (`WT`)** against the prioritized pathogenic hinge-region variant **`L119R`**.

```bash
# 1. Build WT CDKL5 catalytic kinase domain (residues 1-303 matching UniProt O76039)
cmda modeller build --pdb-id 4bgq --uniprot-id O76039 --chain A --range 1 303 --uniprot-numbering --outdir test-project/data/structure/protein

# 2. Mutate to create prioritized pathogenic variant (L119R) preserving biological numbering
cmda modeller mutate --pdb-in test-project/data/structure/protein/target.B99990001_with_cryst.pdb --chain A --mut L119R --outdir-mut test-project/data/structure/protein
```

---

### 4.4. Track 1: Unliganded Apo Dynamics Pipeline (`--state apo`)

#### 4.4.1. Biophysical Scope & Input Requirements
The **Apo state** represents the unliganded kinase domain in explicit OPC water with neutralizing counterions ($0.15\text{ M NaCl}$). It establishes the baseline conformational plasticity and quantifies intrinsic mutation-induced destabilization in the N-lobe (residues 13–56) and activation loop (residues 151–191) without cofactor stabilization.

#### 4.4.2. Staging Project Directory Trees
Scaffold independent directory hierarchies and batch scripts for `WT` and `L119R`:
```bash
cmda stage --state apo --variant WT --root test-project
cmda stage --state apo --variant L119R --root test-project
```

#### 4.4.3. Apo Receptor Preparation & Parameterization
`cmda prep` automatically detects the modeled receptors in `data/structure/protein/`, evaluates 3D electrostatic microenvironments with PROPKA 3 at physiological pH 7.4, strips unneeded crystallographic waters, and executes Amber `tleap` in <2 seconds:
```bash
# Amber ff19SB + OPC 4-point water + randomized 0.15 M NaCl neutralization
cmda prep --state apo --variant WT --root test-project
cmda prep --state apo --variant L119R --root test-project
```
*Outputs in `data/simulation/apo/<variant>/prep/`:* `model.pdb`, `solvated.pdb`, `system.prmtop`, `system.inpcrd`, and `prep.json` (audit manifest).

#### 4.4.4. 24-Stage Multi-Step Equilibration
Runs the rigorous 24-stage VarMDyn relaxation protocol (minimization $\rightarrow$ gradual thermalization to 300 K $\rightarrow$ staged restraint release on backbone and sidechains down to unconstrained NPT):
```bash
cmda sim --state apo --variant WT --stage eq --root test-project
cmda sim --state apo --variant L119R --stage eq --root test-project
```
*Outputs in `data/simulation/apo/<variant>/eq/`:* `solvated.pdb`, `eq.pdb` (equilibrated seed for production), `qc.png` (4-panel stability plot), and `qc.json`.

#### 4.4.5. Production MD Triplicates (`r1`, `r2`, `r3`)
Execute resume-safe production replicas, forking independently from `eq/`:
```bash
# Triplicate production MD (100 ns each)
cmda sim --state apo --variant WT --stage md --replica r1,r2,r3 --total-ns 100.0 --root test-project
cmda sim --state apo --variant L119R --stage md --replica r1,r2,r3 --total-ns 100.0 --root test-project
```

---

### 4.5. Track 2: Pre-Catalytic Holo Dynamics Pipeline (`--state atp_mg`)

#### 4.5.1. Biophysical Scope & Active-Site Chemistry
The **Holo state** models CDKL5 with bound ATP and catalytic divalent magnesium ions ($\text{Mg}^{2+}$). This state investigates nucleotide-induced active site closure, the positioning of the catalytic Lys42–Glu60 salt bridge, and how mutations like `L119R` perturb cofactor coordination.

#### 4.5.2. ATP & $\text{Mg}^{2+}$ Coordinate Transfer
Rather than blind docking, ColabMDA replicates the validated `cdkl5-activation` / `varmdyn` protocol:
1. Performs core $\text{C}\alpha$ alignment (residues 30–220) between CDKL5 and the active CDK2 crystal structure (**PDB: `8FP5`**).
2. Transfers the exact coordinates of ATP (`atponly.pdb`) and catalytic $\text{Mg}^{2+}$ ions (`mgonly.pdb`) into the CDKL5 binding cleft.
3. Automatically writes the complex coordinate into `data/simulation/atp_mg/<variant>/prep/model.pdb`.

#### 4.5.3. Staging Holo Directory Trees
```bash
cmda stage --state atp_mg --variant WT --root test-project
cmda stage --state atp_mg --variant L119R --root test-project
```

#### 4.5.4. Holo Parameterization (Hu et al. 2024 Quantum Parameters)
`cmda prep` routes specialized force field parameters for the nucleotide-cofactor complex:
* **Protein:** Amber `ff19SB`
* **Nucleotide:** Hu et al. 2024 quantum-calculated ATP parameters (`hu2024_atp.frcmod` + unit library)
* **Divalent Cations:** Amber 12-6-4 Lennard-Jones parameters for $\text{Mg}^{2+}$
* **Neutralization:** Automatically accounts for the $-4.0e$ net charge of ATP and $+2.0e$ / $+4.0e$ of $\text{Mg}^{2+}$:
```bash
cmda prep --state atp_mg --variant WT --root test-project
cmda prep --state atp_mg --variant L119R --root test-project
```

#### 4.5.5. Cofactor-Aware Equilibration & Production MD
Equilibrates the complex using 24-stage multi-step relaxation with dedicated positional restraints on ATP heavy atoms and $\text{Mg}^{2+}$ coordination geometry during heating:
```bash
# 1. Equilibrate complex in shared eq/ hub
cmda sim --state atp_mg --variant WT --stage eq --root test-project
cmda sim --state atp_mg --variant L119R --stage eq --root test-project

# 2. Run triplicate production MD (100 ns each)
cmda sim --state atp_mg --variant WT --stage md --replica r1,r2,r3 --total-ns 100.0 --root test-project
cmda sim --state atp_mg --variant L119R --stage md --replica r1,r2,r3 --total-ns 100.0 --root test-project
```

---

### 4.6. Progressive State Expansion Roadmap

ColabMDA's declarative taxonomy natively supports the next stages in the CDKL5 activation cascade without changing CLI syntax:

1. **State 3: Post-Catalytic / Autophosphorylated Kinase (`--state adp_2mg_py171`)**
   - Incorporates autophosphorylation of the activation loop at Tyr171 (`pY171` $\rightarrow$ Amber `PTR` parameterized via `phosaa19SB`) alongside Hu 2024 ADP and $2\text{Mg}^{2+}$.
2. **State 4: Active Quaternary Signaling Complex (`--state atp_2mg_py171_peb2`)**
   - Combines phosphorylated kinase `pY171`, quantum ATP, $2\text{Mg}^{2+}$, and continuous geometric docking of the human EB2 Ser222 substrate peptide (`KRPSSAPSG` motif).

---

### 4.7. Trajectory Post-Processing & Comparative Analytics

#### 4.7.1. PBC Autoimaging & Bulk Solvent Stripping (`cmda merge`)
Uses `cpptraj` to center the solute, wrap across periodic boundary conditions, and strip bulk water while preserving ligand and active-site ions (>90% storage savings):
```bash
# Apo trajectories
cmda merge --state apo --variant WT --replicas r1,r2,r3 --root test-project
cmda merge --state apo --variant L119R --replicas r1,r2,r3 --root test-project

# Holo trajectories
cmda merge --state atp_mg --variant WT --replicas r1,r2,r3 --root test-project
cmda merge --state atp_mg --variant L119R --replicas r1,r2,r3 --root test-project
```

#### 4.7.2. Trajectory Metrics & Convergence (`cmda analyze`)
Extracts time-series metrics (Backbone RMSD, Radius of Gyration $R_g$, per-residue RMSF) and publication summary plots:
```bash
cmda analyze --state apo --variant WT --root test-project
cmda analyze --state apo --variant L119R --root test-project
cmda analyze --state atp_mg --variant WT --root test-project
cmda analyze --state atp_mg --variant L119R --root test-project
```

#### 4.7.3. Cross-Variant & Cross-State Comparative Overlays (`cmda compare`)
```bash
# 1. Compare Wild-Type vs Hinge Mutant in Apo state
cmda compare --state apo --variants WT,L119R --root test-project

# 2. Compare Wild-Type vs Hinge Mutant in Holo state
cmda compare --state atp_mg --variants WT,L119R --root test-project
```

> 💡 **Specifying Project Root (`--root`):**
> You can pass either a relative folder in your workspace or a complete path:
> * **Relative path (Workstation / CWD):** `--root test-project`
> * **Complete path (Google Colab):** `--root /content/drive/MyDrive/projects/test-project` (guarantees persistent storage on Google Drive across VM restarts)
> * **Complete path (HPC Cluster):** `--root /scratch/$USER/projects/test-project` (high-speed NVMe scratch storage)

---

### 4.8. Pilot Run Protocol: 1 ns Equilibration + 1 ns Production (WT & L119R)

Before initiating full 100 ns production runs, execute a standardized **1 ns Equilibration + 1 ns Production Pilot Run** to benchmark GPU performance, verify active-site coordination, and validate the automated autoimaging and analysis pipelines.

#### 4.8.1. Pilot Run Benchmark Parameters
* **Equilibration Duration:** 1 ns total (`--equil-time 1000` = 1,000 ps across NVT thermalization and NPT density equilibration).
* **Production Duration:** 1 ns total (`--total-ns 1.0` = 500,000 steps at 2.0 fs timestep).
* **Trajectory Frame Interval:** `--traj-interval 10` (10 ps = 100 frames) or `--traj-interval 1.0` (1,000 frames).
* **Expected Walltime:** ~1.5–2.5 minutes per replica on modern NVIDIA GPUs (T4 / V100 / A100 / RTX).

#### 4.8.2. Environment 1: Google Colab (`--root /content/drive/MyDrive/projects/test-project`)
On Google Colab, mount Google Drive to guarantee persistent storage across VM disconnections:

```bash
# ==============================================================================
# Google Colab Pilot Execution (WT & L119R across Apo and Holo atp_mg)
# ==============================================================================
PROJECT_DIR="/content/drive/MyDrive/projects/test-project"

# 1. Stage project directories and scripts
cmda stage --state apo --variant WT --root "${PROJECT_DIR}"
cmda stage --state apo --variant L119R --root "${PROJECT_DIR}"
cmda stage --state atp_mg --variant WT --root "${PROJECT_DIR}"
cmda stage --state atp_mg --variant L119R --root "${PROJECT_DIR}"

# 2. Build parameter topologies (Amber ff19SB / Hu 2024 / OPC 4-point water)
cmda prep --state apo --variant WT --root "${PROJECT_DIR}"
cmda prep --state apo --variant L119R --root "${PROJECT_DIR}"
cmda prep --state atp_mg --variant WT --root "${PROJECT_DIR}"
cmda prep --state atp_mg --variant L119R --root "${PROJECT_DIR}"

# 3. 1 ns Equilibration Hub (EM -> NVT -> NPT -> QC)
cmda sim --state apo --variant WT --stage eq --equil-time 1000 --root "${PROJECT_DIR}"
cmda sim --state apo --variant L119R --stage eq --equil-time 1000 --root "${PROJECT_DIR}"
cmda sim --state atp_mg --variant WT --stage eq --equil-time 1000 --root "${PROJECT_DIR}"
cmda sim --state atp_mg --variant L119R --stage eq --equil-time 1000 --root "${PROJECT_DIR}"

# 4. 1 ns Production MD Pilot Replicas (r1)
cmda sim --state apo --variant WT --stage md --replica r1 --total-ns 1.0 --traj-interval 10 --root "${PROJECT_DIR}"
cmda sim --state apo --variant L119R --stage md --replica r1 --total-ns 1.0 --traj-interval 10 --root "${PROJECT_DIR}"
cmda sim --state atp_mg --variant WT --stage md --replica r1 --total-ns 1.0 --traj-interval 10 --root "${PROJECT_DIR}"
cmda sim --state atp_mg --variant L119R --stage md --replica r1 --total-ns 1.0 --traj-interval 10 --root "${PROJECT_DIR}"

# 5. Automated Autoimaging & Trajectory Post-Processing
cmda merge --state apo --variant WT --replicas r1 --root "${PROJECT_DIR}"
cmda merge --state apo --variant L119R --replicas r1 --root "${PROJECT_DIR}"
cmda merge --state atp_mg --variant WT --replicas r1 --root "${PROJECT_DIR}"
cmda merge --state atp_mg --variant L119R --replicas r1 --root "${PROJECT_DIR}"

# 6. Biophysical Analysis & Comparative Overlays
cmda analyze --state apo --variant WT --root "${PROJECT_DIR}"
cmda analyze --state apo --variant L119R --root "${PROJECT_DIR}"
cmda analyze --state atp_mg --variant WT --root "${PROJECT_DIR}"
cmda analyze --state atp_mg --variant L119R --root "${PROJECT_DIR}"
cmda compare --state apo --variants WT,L119R --root "${PROJECT_DIR}"
cmda compare --state atp_mg --variants WT,L119R --root "${PROJECT_DIR}"
```

#### 4.8.3. Environment 2: Palmetto HPC Scratch (`--root /scratch/$USER/projects/test-project`)
On the Clemson University Palmetto cluster, execute high-throughput simulations on scratch NVMe/Lustre storage:

```bash
# ==============================================================================
# Palmetto HPC Scratch Pilot Execution
# ==============================================================================
export SCRATCH_PROJECT="/scratch/$USER/projects/test-project"

# 1. Stage project directories and generate SLURM runner scripts
cmda stage --state apo --variant WT --root "${SCRATCH_PROJECT}"
cmda stage --state apo --variant L119R --root "${SCRATCH_PROJECT}"
cmda stage --state atp_mg --variant WT --root "${SCRATCH_PROJECT}"
cmda stage --state atp_mg --variant L119R --root "${SCRATCH_PROJECT}"

# 2. Build parameter topologies
cmda prep --state apo --variant WT --root "${SCRATCH_PROJECT}"
cmda prep --state apo --variant L119R --root "${SCRATCH_PROJECT}"
cmda prep --state atp_mg --variant WT --root "${SCRATCH_PROJECT}"
cmda prep --state atp_mg --variant L119R --root "${SCRATCH_PROJECT}"

# 3. Interactive Pilot Execution on GPU Node (salloc / qsub)
# Request GPU node: salloc --nodes=1 --ntasks=1 --cpus-per-task=8 --gpus=1 --mem=32gb --time=01:00:00
cmda sim --state apo --variant WT --stage eq --equil-time 1000 --root "${SCRATCH_PROJECT}"
cmda sim --state apo --variant WT --stage md --replica r1 --total-ns 1.0 --traj-interval 10 --root "${SCRATCH_PROJECT}"
cmda sim --state atp_mg --variant WT --stage eq --equil-time 1000 --root "${SCRATCH_PROJECT}"
cmda sim --state atp_mg --variant WT --stage md --replica r1 --total-ns 1.0 --traj-interval 10 --root "${SCRATCH_PROJECT}"

# 4. Or Submit via Chained SLURM Batch Helpers (EQ -> Array Prod)
cd "${SCRATCH_PROJECT}"
bash scripts/submit_chain.sh
```

#### 4.8.4. Environment 3: Modular OpenMM CLI Execution (`<pdbid>_<variant>`)
For users employing the modular OpenMM engine with `<pdbid>_<variant>` convention:

```bash
# 1. 1 ns Equilibration in shared eq/ hub (4bgq_wt)
colabmda openmm stage --pdb-file test-project/data/structure/4bgq/wt/prep/4bgq_wt_cleaned.pdb --name 4bgq_wt --replica r1
colabmda openmm em --name 4bgq_wt --protein-ff ff19SB --water-model opc
colabmda openmm nvt --name 4bgq_wt --equil-time 500
colabmda openmm npt --name 4bgq_wt --equil-time 500
colabmda openmm qc --name 4bgq_wt

# 2. 1 ns Production MD replica 1
colabmda openmm md --name 4bgq_wt --replica r1 --total-ns 1.0 --traj-interval 10

# 3. Repeat for Hinge Mutant (4bgq_l119r)
colabmda openmm stage --pdb-file test-project/data/structure/4bgq/mutants/4bgq_l119r/prep/4bgq_l119r_cleaned.pdb --name 4bgq_l119r --replica r1
colabmda openmm em --name 4bgq_l119r --protein-ff ff19SB --water-model opc
colabmda openmm nvt --name 4bgq_l119r --equil-time 500
colabmda openmm npt --name 4bgq_l119r --equil-time 500
colabmda openmm qc --name 4bgq_l119r
colabmda openmm md --name 4bgq_l119r --replica r1 --total-ns 1.0 --traj-interval 10
```

---

## 5. Granular & Legacy Workflows (`colabmda openmm` / `modeller`)

For users requiring fine-grained control over individual modeling steps or legacy PDB workflows:

### 5.1. Homology Modeling & Point Mutations (Modeller)
**Environment:** `modeller_env`

The Modeller workflow includes **Biological Residue Numbering** and **Automated Quality Control**:
*   **`--uniprot-numbering`**: Physically renumbers residues to match UniProt biological indexing.
*   **Automatic Alignment Summary**: Confirms sequence alignment before building.
*   **Post-Build QC**: Verifies modeled coordinates and stereochemical energies.

```bash
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate modeller_env

# 1. Build Wild-Type CDKL5 (Residues 1–303 matching UniProt O76039)
colabmda modeller build --pdb-id 4bgq --uniprot-id O76039 --chain A --range 1 303 --uniprot-numbering --outdir test-project/data/structure/4bgq/wt

# 2. Create Point Mutant (e.g. L119R) preserving biological numbering
colabmda modeller mutate --pdb-in test-project/data/structure/4bgq/wt/target.B99990001_with_cryst.pdb --chain A --mut L119R --outdir-mut test-project/data/structure/4bgq/mutants/4bgq_l119r
```

---

### 5.2. Structure Preparation & 3D pKa Titration (PROPKA at pH 7.4)
**Environment:** `openmm_env`

Accurate protonation states are critical for biophysical simulation fidelity. ColabMDA natively integrates **PDBFixer**, **PDB2PQR**, and **PROPKA 3** to calculate local 3D electrostatic microenvironments and assign explicit AMBER titratable residue forms (`HID`, `HIE`, `HIP`, `ASH`, `GLH`, `LYN`, `CYX`):

```bash
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate openmm_env

# Run 3D pKa titration at physiological pH 7.4
colabmda openmm prep --pdb-file test-project/data/structure/4bgq/wt/4bgq_wt.pdb --outdir test-project/data/structure/4bgq/wt/prep --name 4bgq_wt --ph 7.4
colabmda openmm prep --pdb-file test-project/data/structure/4bgq/mutants/4bgq_l119r/4bgq_l119r.pdb --outdir test-project/data/structure/4bgq/mutants/4bgq_l119r/prep --name 4bgq_l119r --ph 7.4
```

*Key features:*
* **Automated Coordinate Repair:** PDBFixer resolves missing heavy atoms, non-standard residues, and chain breaks.
* **3D Electrostatic Titration:** Evaluates local hydrogen-bonding networks using PROPKA 3 at target pH (default: `7.4`).
* **Explicit AMBER Forms:** Differentiates neutral Histidine tautomers (`HID`/`HIE`), charged Histidine (`HIP`), protonated Asp/Glu (`ASH`/`GLH`), neutral Lysine (`LYN`), and disulfide bridges (`CYX`).
* **Audit Manifest:** Produces `protonation_summary.json` alongside `<name>_cleaned.pdb` recording every assigned titration state.

---

### 5.3. Simulation Execution: Modular Stages vs. All-in-One

ColabMDA provides both fine-grained single-keyword commands and an automated all-in-one runner:

#### Single-Keyword Modular Stages

| Stage | Command | Description |
| :--- | :--- | :--- |
| **Stage** | `colabmda openmm stage --pdb-file <file> --name <system> --replica r1` | Initialize folder structure in `data/simulation/` and `eq/`. |
| **EM** | `colabmda openmm em --name <system> --protein-ff ff19SB --water-model opc` | Energy minimization & box solvation in `eq/`. |
| **NVT** | `colabmda openmm nvt --name <system> --equil-time 100` | NVT thermalization to 300 K. |
| **NPT** | `colabmda openmm npt --name <system> --equil-time 100` | NPT density equilibration with Monte Carlo barostat. |
| **QC** | `colabmda openmm qc --name <system>` | Stability verification & publication QC plots (`qc.png`). |
| **MD** | `colabmda openmm md --name <system> --replica r1 --total-ns 100 --traj-interval 10` | Resume-safe chunked production MD. |
| **Merge** | `colabmda openmm merge --name <system> --replica r1 --center --wrap` | Concatenate chunks into centered, PBC-wrapped trajectory. |
| **Status** | `colabmda openmm status --name <system> --replica r1` | Verify chunks, topology, frame counts, and resume readiness. |
| **Analysis** | `colabmda openmm analysis --name <system> --replica r1` | Single-system trajectory metrics (RMSD, $R_g$, RMSF). |

#### Step-by-Step Modular Execution:

```bash
# 1. Initialize simulation folder and eq hub (WT)
colabmda openmm stage --pdb-file test-project/data/structure/4bgq/wt/prep/4bgq_wt_cleaned.pdb --name 4bgq_wt --replica r1

# 2. Energy Minimization in shared eq/ hub (ff19SB / OPC 4-point water)
colabmda openmm em --name 4bgq_wt --protein-ff ff19SB --water-model opc

# 3. NVT Thermalization (100 ps heating to 300 K)
colabmda openmm nvt --name 4bgq_wt --equil-time 100

# 4. NPT Density Equilibration (100 ps with Monte Carlo barostat)
colabmda openmm npt --name 4bgq_wt --equil-time 100

# 5. Equilibration Quality Control Gating
colabmda openmm qc --name 4bgq_wt

# 6. Production MD across Triplicates (r1, r2, r3 fork automatically from eq/)
colabmda openmm md --name 4bgq_wt --replica r1 --total-ns 100.0 --traj-interval 10
colabmda openmm md --name 4bgq_wt --replica r2 --total-ns 100.0 --traj-interval 10
colabmda openmm md --name 4bgq_wt --replica r3 --total-ns 100.0 --traj-interval 10

# 7. Merge Trajectory Chunks for all replicas (Center + Wrap)
colabmda openmm merge --name 4bgq_wt --replica r1 --center --wrap
colabmda openmm merge --name 4bgq_wt --replica r2 --center --wrap
colabmda openmm merge --name 4bgq_wt --replica r3 --center --wrap

# 8. Check Status & Frame Counts
colabmda openmm status --name 4bgq_wt --replica r1

# 9. Biophysical Trajectory Analysis (RMSD, Rg, RMSF)
colabmda openmm analysis --name 4bgq_wt --replica r1
colabmda openmm analysis --name 4bgq_wt --replica r2
colabmda openmm analysis --name 4bgq_wt --replica r3
```

#### All-in-One Alternative (`run`):
For automated, hands-off staging and execution from EM through Production MD in one call:
```bash
colabmda openmm run --name 4bgq_wt --replica r1 --total-ns 100.0 --traj-interval 10 --equil-time 1000 --checkpoint-ps 1000
```

> 💡 **Storage Tip:** For a typical system (e.g., CDKL5 kinase domain in water, ~35,000 atoms), a 100ns run at high resolution (1ps) can produce over **40GB** of data. On a free 15GB Google Drive, we recommend using **`--traj-interval 10`** to reduce this to ~4GB. Always calculate your storage needs based on your specific system size before starting long runs.

> 🛡️ **Stability Gate:** The `qc` command (and built-in gating in `run`) automatically analyzes equilibration logs and aborts if the system has not stabilized, preventing wasted cloud GPU time.

---

### 5.4. Merge and Center Trajectories

Combine trajectory chunks into a single DCD, apply periodic boundary condition (PBC) correction, and center the protein using the robust MDAnalysis-based engine:

```bash
# Standard Merge (Center + Wrap, All Atoms)
colabmda openmm merge --name 4bgq_wt --replica r1 --center --wrap
colabmda openmm merge --name 4bgq_wt --replica r2 --center --wrap
colabmda openmm merge --name 4bgq_wt --replica r3 --center --wrap

# For Mutant (L119R):
colabmda openmm merge --name 4bgq_l119r --replica r1 --center --wrap
colabmda openmm merge --name 4bgq_l119r --replica r2 --center --wrap
colabmda openmm merge --name 4bgq_l119r --replica r3 --center --wrap
```

#### Output Files:
*   `prod_full.dcd`: The final, concatenated, and centered/wrapped trajectory file.
*   `prod_full.log`: The consolidated log file containing energy and temperature statistics.

> 💡 **Pro-Tip for Long Runs:**
> Merging processes trajectories frame-by-frame, so it will not crash your RAM. You can merge without striding (`--stride 1`) for full resolution, or use `--stride 10` to create a lightweight file for local viewing.
>
> When you merge with `--stride 10`, the logs are automatically strided to match, and the analysis command (`colabmda openmm analysis`) will read the logs to correctly infer the time scale.

> 💡 **Protein-Only Trajectory (Optional):**
> To save >85% storage, add the `--protein-only` flag to extract only protein atoms, discarding solvent water and ions:
> `colabmda openmm merge --name 4bgq_wt --replica r1 --center --wrap --protein-only`

> **Outputs with `--protein-only`:**
> * `prod_full.dcd`: Thinned trajectory containing protein atoms only.
> * `prod_full.pdb`: Matching protein-only topology PDB file (crucial for subsequent analysis/visualization to avoid atom mismatch).
> * `prod_full.log`: Consolidated log file.

---

### 5.5. Post-Merge Guidelines & Verification (FAQ)

#### 5.5.1. 🔍 Verifying Trajectory Frames
To quickly check the status, simulation time, and number of frames in your merged files, run:
```bash
colabmda openmm status --name 4bgq_wt --replica r1
```
This prints a comprehensive status report including topology stats, chunks, and exact frame counts:
```text
[STATUS]
  Workdir          : data/simulation/4bgq_wt/r1
  Chunks (DCD/log) : 100 / 100
  Topology File    : solvated.pdb
                     └─ 33450 atoms, 10240 residues
                     └─ (303 protein, 10120 water, 62 ions)
  Trajectory File  : prod_full.dcd (10000 frames)
  Log File         : prod_full.log (10000 frames)
  Frames (from logs): 10000
```

Alternatively, query the frame count directly with MDAnalysis:
```bash
python3 -c "import MDAnalysis as mda; u = mda.Universe('data/simulation/4bgq_wt/r1/prod_full.pdb', 'data/simulation/4bgq_wt/r1/prod_full.dcd'); print('Frames:', len(u.trajectory))"
```

#### 5.5.2. ⏱️ Understanding Simulation Time vs. Frame Counts
Use this formula to calculate expected frames:
$$\text{Expected Frames} = \frac{\text{Total Time (ps)}}{\text{Trajectory Saving Interval (ps)}}$$
*Example:* A $100.0\text{ ns}$ ($100,000\text{ ps}$) simulation with a $10.0\text{ ps}$ saving interval:
$$\text{Expected Frames} = \frac{100,000\text{ ps}}{10\text{ ps}} = 10,000\text{ frames}$$

*Striding on Merged Trajectories:*
$$\text{Merged Frames} = \frac{\text{Total Frames}}{\text{Stride}} = \frac{10,000}{10} = 1,000\text{ frames}$$

#### 5.5.3. 🧬 Reference Topology Guidelines
When merging trajectories, the reference topology file used and the resulting output files depend on your merge options:

| Merge Mode | Command Flags | Reference Topology | Generated Topology File | Subsequent Command Usage |
| :--- | :--- | :--- | :--- | :--- |
| **Standard Merge** | `colabmda openmm merge --center --wrap` | `solvated.pdb` | *None* (only `prod_full.dcd` is written) | Use `solvated.pdb` for analysis/visualization. |
| **MDAnalysis Merge** | `colabmda openmm merge --mda --center --wrap` | `solvated.pdb` | `prod_full.pdb` (all atoms) | Use `prod_full.pdb` (or `solvated.pdb`). |
| **Protein-Only Merge** | `colabmda openmm merge --protein-only` | `solvated.pdb` | `prod_full.pdb` (protein atoms only) | **Must** use `prod_full.pdb` (since `prod_full.dcd` contains only protein coordinates). |

> ⚠️ **Why does this matter?**
> If you merge using `--protein-only`, `prod_full.dcd` contains only protein coordinates (~4,800 atoms). Attempting to load this trajectory along with `solvated.pdb` (~33,000 atoms) in PyMOL or MDAnalysis will result in a **fatal atom mismatch error**. Always match the trajectory file with its corresponding topology file as shown in the table above.

---

### 5.6. Analysis & Comparison

#### 5.6.1. Single System Analysis (RMSD, Rg, RMSF)
```bash
# For Wild-Type Triplicates:
colabmda openmm analysis --name 4bgq_wt --replica r1
colabmda openmm analysis --name 4bgq_wt --replica r2
colabmda openmm analysis --name 4bgq_wt --replica r3

# For Mutant Triplicates (L119R):
colabmda openmm analysis --name 4bgq_l119r --replica r1
colabmda openmm analysis --name 4bgq_l119r --replica r2
colabmda openmm analysis --name 4bgq_l119r --replica r3
```

> ⚠️ **Analysis Tip:** If your plots show the wrong time scale (e.g., 10ns instead of 100ns), provide the frame interval manually with `--interval 10`.

#### 5.6.2. WT vs Mutant Comparison (Triplicate Statistics)
```bash
colabmda openmm compare \
  --series "WT=data/analysis/4bgq_wt/r1,data/analysis/4bgq_wt/r2,data/analysis/4bgq_wt/r3" \
  --series "L119R=data/analysis/4bgq_l119r/r1,data/analysis/4bgq_l119r/r2,data/analysis/4bgq_l119r/r3" \
  --outdir data/analysis/compare/wt_vs_l119r
```

---

## 6. Trajectory Visualization in PyMOL

To visualize molecular dynamics trajectories in PyMOL with full secondary structure (ribbon/cartoon) and detailed sidechain representations, use the built-in `colabmda openmm view` command.

### 6.1. Local Setup and Installation
Since molecular dynamics simulations are performed on Google Colab, trajectory visualization is run locally on your workstation/laptop.

> [!NOTE]
> **For Windows Users:** It is highly recommended to run visualization locally either using **WSL (Windows Subsystem for Linux)** or using native **PyMOL on Windows** (by running the CLI to generate `visualize.pml`, then opening it in the Windows PyMOL GUI).

Create a dedicated conda environment containing PyMOL and ColabMDA:

```bash
# 1. Create a clean environment and install PyMOL
conda create -y -n colabmda_env -c conda-forge python=3.11 pymol-open-source

# 2. Activate the environment
conda activate colabmda_env

# 3. Install ColabMDA package (Lightweight local installation)
python3 -m pip install --upgrade "git+https://github.com/paulshamrat/ColabMDA.git@main"

# 4. Verify PyMOL installation
pymol --version
```

### 6.2. Running the Visualization in PyMOL
Run the view command from your simulation folder (it looks for `prod_full.pdb` and `prod_full.dcd`, generates a PyMOL script, and launches PyMOL):

```bash
# For Wild-Type:
colabmda openmm view --pdb-dir data/simulation/4bgq_wt/r1

# For Mutant (L119R, highlighting residue 119):
colabmda openmm view --pdb-dir data/simulation/4bgq_l119r/r1 --resi 119
```

The command automatically generates a `visualize.pml` script in the folder.

<details>
<summary><b>Click to view the generated PyMOL script configuration (<code>visualize.pml</code>)</b></summary>

```text
# 1. Clear out memory
reinitialize
bg_color white

# 2. Load merged trajectory files
load prod_full.pdb, 4bgq
load_traj prod_full.dcd, 4bgq

# 2b. Align trajectory to reference (e.g., 4bgq_wt.pdb)
load 4bgq_wt.pdb, crystal_ref
align 4bgq and name CA, crystal_ref and name CA, mobile_state=1, target_state=1
delete crystal_ref

# 3. Hide solvent water and ions
hide everything, all
hide nonbonded, all
hide nb_spheres, all

# 4. Freeze backbone tumbling frame-by-frame
intra_fit 4bgq and name CA

# 5. Generate crisp secondary structure ribbon
dss 4bgq
cartoon automatic, 4bgq
show cartoon, 4bgq

# 6. Apply smooth cyan color scheme
color cyan, 4bgq and name C*
util.cnc("4bgq")

# 7. Highlight mutated residue 119 side chain sticks
show sticks, 4bgq and resi 119 and not name N+C+O+H
color yellow, 4bgq and resi 119 and name NH1+NH2
set stick_radius, 0.25

# 8. Focus camera onto protein
zoom 4bgq and polymer, buffer=4
```

</details>

### 6.3. Generating Comparative Snapshot Grids
Render publication-quality structural snapshot grids comparing transitions over trajectories (such as WT vs Mutants across specific frames):

```bash
# Generate the default transition snapshot grid (WT vs L119R)
colabmda openmm snapshots
```

#### Customization via JSON Configuration:
```bash
colabmda openmm snapshots --config my_config.json --output figures/comparison_grid.png
```

```json
{
  "align_ref_pdb": "test-project/data/structure/4bgq/wt/target.B99990001_with_cryst.pdb",
  "stable_core_sel": "resi 20-35 or resi 85-110 or resi 180-260",
  "camera_view": [
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0,
    1.0, 0.0, 0.0,
    0.0, 0.0, -50.0
  ],
  "systems": {
    "WT": {
      "pdb": "data/simulation/4bgq_wt/r1/prod_full.pdb",
      "dcd": "data/simulation/4bgq_wt/r1/prod_full.dcd",
      "states": [1, 200],
      "times": ["0.00 ns", "2.00 ns"],
      "mut_residue": 119
    }
  }
}
```

---

## Acknowledgements
- [OpenMM](https://openmm.org) & [PDBFixer](https://github.com/openmm/pdbfixer)
- [Modeller](https://salilab.org/modeller/)
- [MDAnalysis](https://www.mdanalysis.org) & [MDTraj](https://www.mdtraj.org)
- [PDB2PQR & PROPKA](https://server.poissonboltzmann.org/)
- [NumPy](https://numpy.org), [Matplotlib](https://matplotlib.org), [Biopython](https://biopython.org)
- [Google Colab](https://colab.research.google.com) & [Miniforge/Conda](https://github.com/conda-forge/miniforge)

---

## Changelog

### v0.3.0 (Kinase-Substrate MD & Unified Architecture)
*   **Unified `cmda` CLI:** Single-keyword entrypoints (`stage`, `prep`, `sim`, `merge`, `analyze`, `compare`, `colab`) accessible via both `cmda` and `colabmda`.
*   **Kinase-Substrate Dynamics Engine:** Automated CDKL5 phosphotyrosine (`pY171` -> Amber `PTR` with `phosaa19SB`), quantum-parameterized nucleotides (Hu 2024 ATP/ADP/GDP), catalytic dual divalent ions ($Mg^{2+}$, $Mn^{2+}$), and continuous geometric substrate peptide docking (EB2 `KRPSSAPSG` at Ser222).
*   **High-Performance AmberTools Integration:** Automated Amber `tleap` parameterization with randomized ion neutralization (`addIonsRand`, solvating in <2s) and automated `cpptraj` PBC autoimaging and dual-stream solvent stripping (>90% storage savings).
*   **Strict Decoupled Architecture:** Guaranteed zero data files in `src/`, decoupled execution scripts in `<project>/scripts/` (`run.sh`, `submit.slurm`), and a single bootstrap script (`scripts/setup.sh`).
*   **Comprehensive Test Suite:** 40+ unit and regression tests passing across all modeling, parameterization, and simulation modules.

### v0.2.0 (Modular Refactor & Resilience)
*   **Single-Keyword Nomenclature:** Unified modular commands (`prep`, `stage`, `em`, `nvt`, `npt`, `qc`, `md`) and standardized `equil/` directory layout.
*   **PROPKA 3D pKa Titration:** Integrated PDBFixer and PDB2PQR/PROPKA 3 for physiological pH 7.4 protonation state assignment.
*   **Equilibration QC Gating:** Added automated `colabmda openmm qc` module generating publication-grade stability verification (`qc.png`).
*   **Decoupled Modules:** Isolated virtual screening workflows onto dedicated `docking` branch for lightweight performance.
*   **Automated Test Suite:** Added unit and regression test suite under `tests/` (`test_cli.py`, `test_forcefields.py`, `test_qc.py`).
*   **Force Field Resilience:** Integrated Amber14/ff19SB fallbacks and OPC 4-point water support.

### v0.1.0 (Initial Beta)
*   **Modular Pipeline:** New modular CLI for EM, NVT, NPT, and Production MD.
*   **Resume-Safe Engine:** Integrated checkpointing logic for fail-safe simulations on Google Colab.
*   **Modeling:** Automated Wild-Type building and mutation support via Modeller.
*   **Analysis:** Robust trajectory merging and comparative RMSD/Rg/RMSF analysis tools.
*   **Professional Standards:** Added CI/CD workflows, Black formatting, and Ruff linting.

---

## Citation

Initially, this repository was based on a GROMACS-on-Colab workflow, but we transformed it to utilize OpenMM to improve feasibility, performance, and stability for automated simulation recovery on cloud platforms.

If you use this tool, please consider citing the underlying studies:

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
