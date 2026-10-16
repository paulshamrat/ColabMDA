# ColabMDA

**ColabMDA** is a high-performance Molecular Dynamics (MD) pipeline designed for Google Colab and HPC clusters. It eliminates lost simulation progress caused by cloud disconnections through automated Google Drive checkpointing and instant one-command resumption. ColabMDA streamlines the entire workflow—from homology modeling and physiological protonation (PROPKA 3) to AmberTools parameterization, OpenMM simulation, trajectory autoimaging (`cpptraj`), and publication analysis.

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
| **Dependencies** | `OpenMM`, `AmberTools (tleap, cpptraj)`, `Modeller`, `MDAnalysis`, `MDTraj`, `PDB2PQR / PROPKA` |
| **Platform** | `Linux` · `Google Colab (CLI / Terminal)` · `HPC (SLURM)` |
| **Structure** | [Source Layout](#-project-structure) |

### 📂 Project Structure
ColabMDA is organized into clear, functional modules:
*   `src/colabmda/cli.py`: Central unified CLI (`cmda` and `colabmda`).
*   `src/colabmda/core/`: System taxonomy (`SystemSpec`), project layout (`ProjectPaths`), and divalent ion parameters.
*   `src/colabmda/amber/`: AmberTools parameterization (`tleap`, Hu 2024 quantum nucleotide parameters, `cpptraj` autoimaging & solvent stripping).
*   `src/colabmda/modeller/`: Homology modeling engine (UniProt biological numbering, mutagenesis, geometric substrate docking).
*   `src/colabmda/openmm/`: OpenMM simulation engines (Modular prep, EM, NVT, NPT, QC gating, chunked MD, analysis).
*   `src/colabmda/remote/`: Google Colab GPU remote session management via `colab` CLI.
*   `scripts/`: Single standalone bootstrap script (`scripts/setup.sh`).
*   `tests/`: Comprehensive unit and regression test suite (`pytest`).
*   `docs/`: Sphinx and ReadTheDocs scientific documentation source.

---

## 1. Platform Setup & Environment Provisioning (Google Colab & HPC)

ColabMDA is engineered to run seamlessly across cloud instances, high-performance computing (HPC) clusters, and local Linux workstations.

### 1.1. Environment A: Google Colab (Remote Cloud GPU)

#### 1.1.1. Local Terminal Runner (Recommended — Zero Browser)
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

#### 1.1.2. Web Browser Notebook Interface Alternative
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

#### 1.1.3. Timeout Immunity & Crash Recovery Logic
* **Timeout-Proof Execution:** Running simulations via `cmda colab exec` keeps the remote kernel in a `BUSY` execution state, permanently preventing Colab's 30-minute idle VM reaper from shutting down your run.
* **Granular Checkpoint Recovery:** If a Colab session expires at its 12-hour platform limit, simply restart the command. ColabMDA automatically detects saved `stage_*.state.xml` checkpoints on Google Drive and fast-forwards directly to the unfinished step.

#### 1.1.4. Modeller Academic License Provisioning
A Modeller license key is **required** to run homology modeling and point-mutation workflows (`cmda modeller`).
* 🔑 **Get a Free Academic License:** Register at [salilab.org/modeller/registration.html](https://salilab.org/modeller/registration.html) (Licenses are free for academic use and sent instantly via email).
* To provision MODELLER during bootstrap:
  ```bash
  KEY_MODELLER="YOUR_LICENSE_KEY" WITH_MODELLER=1 cmda colab setup
  ```

---

### 1.2. Environment B: High-Performance Computing (HPC / SLURM)

#### 1.2.1. Cluster Conda Activation & Scratch Storage
On HPC clusters running SLURM, run simulations on high-speed scratch storage (`/scratch/$USER/projects/<project>`):

```bash
# 1. Activate conda environment
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate colabmda

# 2. Define scratch project workspace
export SCRATCH_PROJECT="/scratch/$USER/projects/test-project"
```

#### 1.2.2. Interactive GPU Node Allocation (`salloc` / `qsub`)
For testing or pilot runs, request an interactive GPU compute node:
```bash
salloc --nodes=1 --ntasks=1 --cpus-per-task=8 --gpus=1 --mem=32gb --time=02:00:00
```

#### 1.2.3. Decoupled Batch Submission (SLURM Array Jobs)
ColabMDA automatically generates decoupled SLURM batch scripts under `<project>/scripts/`:
* `submit_eq.slurm`: Single GPU job executing equilibration (EM -> NVT -> NPT -> QC).
* `submit_prod.slurm`: SLURM job array (`--array=1-3`) running triplicates `r1`, `r2`, and `r3` in parallel.
* `submit_chain.sh`: Chained submission helper (`sbatch --dependency=afterok:$EQ_ID`).

Submit the full production campaign with a single command:
```bash
cd "${SCRATCH_PROJECT}"
bash scripts/submit_chain.sh
```

---

### 1.3. Environment C: Local Linux Workstation & Test Verification

#### 1.3.1. Conda / Mamba Environment Provisioning
ColabMDA runs natively on Linux workstations with NVIDIA GPUs. Provision the environment using `environment.yml`:
```bash
git clone https://github.com/paulshamrat/ColabMDA.git
cd ColabMDA
mamba env create -f environment.yml
conda activate colabmda
```

#### 1.3.2. Automated Test Suite Verification
Verify force field parameters, peptide docking, Amber leap parameterization, and CLI execution across all modules:
```bash
pytest tests/ -v
```
*(All 51 unit and regression tests pass cleanly).*

---

## 2. Standardized Project & Data Architecture

ColabMDA enforces strict architectural separation to ensure reproducibility across local workstations, cloud GPUs, and HPC clusters:

1. **Zero Data in `src/`:** Structural inputs (`*.pdb`), topologies (`*.prmtop`), and trajectories belong strictly in the user's project workspace. Non-standard residues like phosphotyrosine (`pY171` -> `PTR171`) are constructed mathematically with zero hardcoded PDBs in source code.
2. **Decoupled Execution Scripts:** Execution shell scripts and SLURM batch scripts reside exclusively in `<project>/scripts/`, never inside simulation data folders.
3. **Single Bootstrap Script:** Standalone environment provisioning is handled by a single script: `scripts/setup.sh`.

```text
<project_root>/                               # E.g. test-project or /content/drive/MyDrive/projects/test-project
├── scripts/                                  # Decoupled execution & batch submission scripts
│   ├── run.sh                                # Local / Colab runner script
│   ├── submit_eq.slurm                       # HPC SLURM equilibration batch script
│   ├── submit_prod.slurm                     # HPC SLURM production array job script
│   └── submit_chain.sh                       # Chained submission helper (EQ -> Array Prod)
│
└── data/                                     # Single root for all domain data
    ├── structure/                            # Raw and modeled structural inputs (Zero raw PDBs in src/)
    │   ├── 4bgq/                             # Target kinase structure (PDB: 4BGQ)
    │   │   ├── wt/                           # Wild-Type receptor (4bgq_wt.pdb, prep/4bgq_wt_cleaned.pdb)
    │   │   └── mutants/4bgq_l119r/           # Pathogenic hinge mutant (4bgq_l119r.pdb, prep/4bgq_l119r_cleaned.pdb)
    │   ├── template/                         # Structural templates (4BGQ.pdb, 8FP5.pdb)
    │   ├── peptide/                          # Substrate peptides (EB2.pdb)
    │   └── ligand/                           # Small molecules & nucleotides (ATP.sdf, ADP.sdf)
    ├── simulation/                           # Staged & executed simulation data
    │   ├── apo/                              # Track 1: Unliganded kinase state
    │   │   ├── WT/                           # (prep/, eq/, r1/, r2/, r3/)
    │   │   └── L119R/                        # (prep/, eq/, r1/, r2/, r3/)
    │   └── atp_mg/                           # Track 2: Pre-catalytic holo state (CDKL5 + ATP + Mg2+)
    │       ├── WT/                           # (prep/, eq/, r1/, r2/, r3/)
    │       └── L119R/                        # (prep/, eq/, r1/, r2/, r3/)
    └── analysis/                             # Trajectory metrics, CSVs, and publication plots
        ├── apo/                              # WT/, L119R/, compare/
        ├── atp_mg/                           # WT/, L119R/, compare/
        └── compare_apo_vs_holo/              # Conformational response to nucleotide binding
```

> 💡 **Why Shared `eq/`?**
> Equilibration (EM, NVT, NPT, QC) only needs to be run once per system. Replicas (`r1`, `r2`, `r3`) fork directly from `eq/npt.state.xml` with distinct random velocity seeds, saving substantial GPU time while maintaining statistical independence.

---

## 3. The Unified Molecular Dynamics Protocol (`cmda`)

ColabMDA provides a streamlined, single-keyword CLI accessible via `cmda` (or `colabmda`) that unifies structure preparation, parameterization, OpenMM simulation dispatch, trajectory post-processing, and multi-replica analysis.

### 3.1. Single-Keyword Command Reference Table

| Command | Syntax | Description |
| :--- | :--- | :--- |
| **`modeller`** | `cmda modeller <build\|mutate>` | Builds homology models (UniProt biological numbering) and introduces point mutations. |
| **`stage`** | `cmda stage --state <tag> --variant <var>` | Generates decoupled project trees and runner scripts (`run.sh`, `submit.slurm`). |
| **`prep`** | `cmda prep --state <tag> --variant <var>` | Heals coordinates, runs PROPKA 3D pKa titration, and parameterizes Amber topology. |
| **`sim`** | `cmda sim --state <tag> --variant <var> --stage <stages>` | Dispatches OpenMM simulations (`em`, `nvt`, `npt`, `qc`, `md`, or composite `eq`). |
| **`merge`** | `cmda merge --state <tag> --variant <var>` | Wraps PBC and strips bulk solvent via `cpptraj` (`stripped.system.prmtop`). |
| **`analyze`** | `cmda analyze --state <tag> --variant <var>` | Computes RMSD, $R_g$, RMSF, PCA, and contact maps. |
| **`compare`** | `cmda compare --state <tag> --variants <v1,v2>` | Renders multi-replica or cross-variant comparative overlay plots. |
| **`colab`** | `cmda colab <new\|setup\|exec\|stop>` | Provisions remote GPU VMs, bootstraps environments, and dispatches simulations. |

---

### 3.2. Foundation: Homology Modeling & Point Mutagenesis (`cmda modeller`)

ColabMDA structures simulation campaigns into clear, progressive biophysical tracks—from unliganded baseline dynamics to fully active signaling complexes. Here we demonstrate the full pipeline using the human **CDKL5 catalytic kinase domain (UniProt `O76039`, residues 1–303, template PDB `4BGQ`)** comparing the **Wild-Type (`WT`)** against the prioritized pathogenic hinge-region variant **`L119R`**.

#### Biological Residue Numbering (`--uniprot-numbering`)
Crystal structures frequently contain expression tags, truncated N-termini, or crystallographic offsets that offset residue numbering (e.g. residue 119 in clinical reports appearing as residue 108 in raw PDB files). 
* The **`--uniprot-numbering`** flag aligns the template sequence to the canonical UniProt sequence (`O76039`) and physically renumbers coordinates so that mutation indices directly match clinical genetics and biochemical literature.

#### Automated Quality Control (QC)
Modeller automated building confirms:
* **Sequence Alignment Summary:** Confirms 100% sequence identity across resolved crystal regions before building coordinates.
* **Stereochemical Energy Gating:** Verifies modeled loop geometries, DOPE scores, and stereochemical restraint satisfaction.

```bash
# 1. Build WT CDKL5 catalytic kinase domain (residues 1-303 matching UniProt O76039)
cmda modeller build --pdb-id 4bgq --uniprot-id O76039 --chain A --range 1 303 --uniprot-numbering --outdir test-project/data/structure/4bgq/wt

# 2. Mutate to create prioritized pathogenic variant (L119R) preserving biological numbering
cmda modeller mutate --pdb-in test-project/data/structure/4bgq/wt/target.B99990001_with_cryst.pdb --chain A --mut L119R --outdir-mut test-project/data/structure/4bgq/mutants/4bgq_l119r
```

---

### 3.3. Track 1: Unliganded Apo Dynamics Pipeline (`--state apo`)

#### 3.3.1. Biophysical Scope & Input Requirements
The **Apo state** represents the unliganded kinase domain in explicit OPC water with neutralizing counterions ($0.15\text{ M NaCl}$). It establishes the baseline conformational plasticity and quantifies intrinsic mutation-induced destabilization in the N-lobe (residues 13–56) and activation loop (residues 151–191) without cofactor stabilization.

#### 3.3.2. Staging Project Directory Trees
Scaffold independent directory hierarchies and batch scripts for `WT` and `L119R`:
```bash
cmda stage --state apo --variant WT --root test-project
cmda stage --state apo --variant L119R --root test-project
```

#### 3.3.3. Chemical Preparation, Coordinate Repair & 3D pKa Titration (`cmda prep`)
Accurate protonation states are critical for biophysical simulation fidelity. `cmda prep` automates coordinate repair and electrostatics:
* **Automated Coordinate Repair:** PDBFixer resolves missing heavy atoms, non-standard residues, and chain breaks.
* **3D Electrostatic Titration:** Evaluates local hydrogen-bonding microenvironments using **PROPKA 3** at physiological pH (default: `7.4`).
* **Explicit AMBER Forms:** Differentiates neutral Histidine tautomers (`HID`/`HIE`), charged Histidine (`HIP`), protonated Asp/Glu (`ASH`/`GLH`), neutral Lysine (`LYN`), and disulfide bridges (`CYX`).
* **Audit Manifest:** Produces `protonation_summary.json` alongside `prep.json` recording every assigned titration state and net system charge.
* **Ultra-Fast Solvation & Neutralization:** Executes Amber `tleap` with Amber `ff19SB`, `OPC` 4-point water, and randomized 0.15 M NaCl ion placement (`addIonsRand`), solvating systems in under 2 seconds.

```bash
# Amber ff19SB + OPC 4-point water + randomized 0.15 M NaCl neutralization
cmda prep --state apo --variant WT --root test-project
cmda prep --state apo --variant L119R --root test-project
```
*Outputs in `data/simulation/apo/<variant>/prep/`:* `model.pdb`, `solvated.pdb`, `system.prmtop`, `system.inpcrd`, and `prep.json`.

#### 3.3.4. Simulation Execution: Modular Stages vs. Composite Hub (`cmda sim`)
ColabMDA allows running equilibration either as an automated all-in-one hub or as granular individual steps:

##### A. Composite All-in-One Hub (`--stage eq`)
Executes the full 24-stage VarMDyn relaxation protocol (minimization $\rightarrow$ gradual thermalization to 300 K $\rightarrow$ staged restraint release down to unconstrained NPT $\rightarrow$ automated QC gating) in one call:
```bash
cmda sim --state apo --variant WT --stage eq --root test-project
cmda sim --state apo --variant L119R --stage eq --root test-project
```

##### B. Granular Modular Execution (`--stage <stages>`)
For fine-grained control over individual thermodynamic ensembles:
| Stage | Syntax | Description |
| :--- | :--- | :--- |
| **Energy Minimization** | `--stage em` | Harmonic restraint relaxation down to $\Delta E < 10.0\text{ kJ/mol/nm}$. |
| **NVT Thermalization** | `--stage nvt` | 100 ps gradual heating to 300 K with Langevin integrator. |
| **NPT Density Equilibration** | `--stage npt` | 100 ps pressure relaxation to 1.0 bar via Monte Carlo barostat. |
| **Stability QC Gating** | `--stage qc` | Automated log parsing and stability verification (`qc.png`, `qc.json`). |
| **Production MD** | `--stage md` | Resume-safe chunked production MD across triplicates (`r1`, `r2`, `r3`). |

```bash
# Step-by-step modular progression
cmda sim --state apo --variant WT --stage em --root test-project
cmda sim --state apo --variant WT --stage nvt --equil-time 100 --root test-project
cmda sim --state apo --variant WT --stage npt --equil-time 100 --root test-project
cmda sim --state apo --variant WT --stage qc --root test-project
```

> 🛡️ **Stability Gate:** The `qc` stage automatically parses equilibration logs (temperature drift, potential energy, density convergence) and aborts execution if the system has not stabilized, preventing wasted cloud GPU time.

> 💡 **Storage Optimization Pro-Tip:** For a typical kinase system in water (~35,000 atoms), a 100 ns run at 1 ps saving frequency produces over **40 GB** of data. Adding **`--traj-interval 10`** (10 ps interval) reduces this to ~4 GB, fitting comfortably on Google Drive while preserving full conformational statistics.

#### 3.3.5. Triplicate Production MD (`r1`, `r2`, `r3`)
Execute resume-safe production replicas, forking independently from `eq/`:
```bash
# Triplicate production MD (100 ns each)
cmda sim --state apo --variant WT --stage md --replica r1,r2,r3 --total-ns 100.0 --traj-interval 10 --root test-project
cmda sim --state apo --variant L119R --stage md --replica r1,r2,r3 --total-ns 100.0 --traj-interval 10 --root test-project
```

---

### 3.4. Track 2: Pre-Catalytic Holo Dynamics Pipeline (`--state atp_mg`)

#### 3.4.1. Biophysical Scope & Active-Site Chemistry
The **Holo state** models CDKL5 with bound ATP and catalytic divalent magnesium ions ($\text{Mg}^{2+}$). This state investigates nucleotide-induced active site closure, the positioning of the catalytic Lys42–Glu60 salt bridge, and how mutations like `L119R` perturb cofactor coordination.

#### 3.4.2. ATP & $\text{Mg}^{2+}$ Coordinate Transfer Protocol
Rather than blind docking, ColabMDA replicates the validated `cdkl5-activation` / `varmdyn` protocol:
1. Performs core $\text{C}\alpha$ alignment (residues 30–220) between CDKL5 and the active CDK2 crystal structure (**PDB: `8FP5`**).
2. Transfers the exact coordinates of ATP (`atponly.pdb`) and catalytic $\text{Mg}^{2+}$ ions (`mgonly.pdb`) into the CDKL5 binding cleft.
3. Automatically writes the complex coordinate into `data/simulation/atp_mg/<variant>/prep/model.pdb`.

#### 3.4.3. Staging Holo Directory Trees
```bash
cmda stage --state atp_mg --variant WT --root test-project
cmda stage --state atp_mg --variant L119R --root test-project
```

#### 3.4.4. Holo Parameterization (Hu et al. 2024 Quantum Parameters)
`cmda prep` routes specialized force field parameters for the nucleotide-cofactor complex:
* **Protein:** Amber `ff19SB`
* **Nucleotide:** Hu et al. 2024 quantum-calculated ATP parameters (`hu2024_atp.frcmod` + unit library)
* **Divalent Cations:** Amber 12-6-4 Lennard-Jones parameters for $\text{Mg}^{2+}$
* **Neutralization:** Automatically accounts for the $-4.0e$ net charge of ATP and $+2.0e$ of $\text{Mg}^{2+}$:
```bash
cmda prep --state atp_mg --variant WT --root test-project
cmda prep --state atp_mg --variant L119R --root test-project
```

#### 3.4.5. Cofactor-Aware Equilibration & Production MD
Equilibrates the complex using 24-stage multi-step relaxation with dedicated positional restraints on ATP heavy atoms and $\text{Mg}^{2+}$ coordination geometry during heating:
```bash
# 1. Equilibrate complex in shared eq/ hub
cmda sim --state atp_mg --variant WT --stage eq --root test-project
cmda sim --state atp_mg --variant L119R --stage eq --root test-project

# 2. Run triplicate production MD (100 ns each)
cmda sim --state atp_mg --variant WT --stage md --replica r1,r2,r3 --total-ns 100.0 --traj-interval 10 --root test-project
cmda sim --state atp_mg --variant L119R --stage md --replica r1,r2,r3 --total-ns 100.0 --traj-interval 10 --root test-project
```

---

### 3.5. Progressive Activation Cascade Roadmap

ColabMDA's declarative taxonomy natively supports the subsequent stages in the CDKL5 activation cascade without changing CLI syntax:

1. **State 3: Post-Catalytic / Autophosphorylated Kinase (`--state adp_2mg_py171`)**
   - Incorporates autophosphorylation of the activation loop at Tyr171 (`pY171` $\rightarrow$ Amber `PTR` parameterized via `phosaa19SB`) alongside Hu 2024 ADP and $2\text{Mg}^{2+}$.
2. **State 4: Active Quaternary Signaling Complex (`--state atp_2mg_py171_peb2`)**
   - Combines phosphorylated kinase `pY171`, quantum ATP, $2\text{Mg}^{2+}$, and continuous geometric docking of the human EB2 Ser222 substrate peptide (`KRPSSAPSG` motif).

---

### 3.6. Trajectory Autoimaging, Solvent Stripping & Topology Matching (`cmda merge`)

#### 3.6.1. Automated PBC Autoimaging & Solvent Stripping
Uses `cpptraj` to center the solute, wrap across periodic boundary conditions, and strip bulk water while preserving ligand and active-site ions (>90% storage savings):
```bash
# Apo trajectories
cmda merge --state apo --variant WT --replicas r1,r2,r3 --root test-project
cmda merge --state apo --variant L119R --replicas r1,r2,r3 --root test-project

# Holo trajectories
cmda merge --state atp_mg --variant WT --replicas r1,r2,r3 --root test-project
cmda merge --state atp_mg --variant L119R --replicas r1,r2,r3 --root test-project
```

#### 3.6.2. 🧬 Reference Topology Guidelines (Preventing Fatal Atom Mismatch)
When analyzing or visualizing trajectories, the topology file used **must strictly match** the trajectory coordinates:

| Trajectory File | Atoms | Reference Topology | Used For |
| :--- | :--- | :--- | :--- |
| **`prod_full.dcd`** | ~33,000 (all atoms) | `solvated.pdb` | Full explicit solvent analysis, water residence times. |
| **`prod_keep_ligand.dcd`** | ~4,900 (solute + ions) | `stripped.system.prmtop` / `prod_full.pdb` | Fast visualization, backbone RMSD, RMSF, PCA. |

> ⚠️ **Critical Warning:** Attempting to load a solvent-stripped trajectory (`prod_keep_ligand.dcd`) alongside `solvated.pdb` in PyMOL or MDAnalysis will result in a **fatal atom mismatch error**. Always pair stripped trajectories with matching stripped topologies.

---

### 3.7. Biophysical Analytics & Trajectory Verification (`cmda analyze` & `cmda compare`)

#### 3.7.1. Single-System Trajectory Metrics (`cmda analyze`)
Extracts time-series metrics (Backbone RMSD, Radius of Gyration $R_g$, per-residue RMSF) and publication summary plots:
```bash
cmda analyze --state apo --variant WT --root test-project
cmda analyze --state apo --variant L119R --root test-project
cmda analyze --state atp_mg --variant WT --root test-project
cmda analyze --state atp_mg --variant L119R --root test-project
```

#### 3.7.2. Cross-Variant & Cross-State Comparative Overlays (`cmda compare`)
```bash
# 1. Compare Wild-Type vs Hinge Mutant in Apo state
cmda compare --state apo --variants WT,L119R --root test-project

# 2. Compare Wild-Type vs Hinge Mutant in Holo state
cmda compare --state atp_mg --variants WT,L119R --root test-project
```

#### 3.7.3. ⏱️ Mathematical Formula: Simulation Time vs. Expected Frame Counts
Use this formula to calculate expected frames:
$$\text{Expected Frames} = \frac{\text{Total Time (ps)}}{\text{Trajectory Saving Interval (ps)}}$$
*Example:* A $100.0\text{ ns}$ ($100,000\text{ ps}$) simulation with a $10.0\text{ ps}$ saving interval:
$$\text{Expected Frames} = \frac{100,000\text{ ps}}{10\text{ ps}} = 10,000\text{ frames}$$

*Striding on Merged Trajectories:*
$$\text{Merged Frames} = \frac{\text{Total Frames}}{\text{Stride}} = \frac{10,000}{10} = 1,000\text{ frames}$$

#### 3.7.4. 🔍 Verification: Fast MDAnalysis Frame Check
To verify frame counts and coordinate integrity directly in Python:
```bash
python3 -c "import MDAnalysis as mda; u = mda.Universe('data/simulation/apo/WT/r1/prod_full.pdb', 'data/simulation/apo/WT/r1/prod_full.dcd'); print('Frames:', len(u.trajectory))"
```

---

### 3.8. Publication-Quality Trajectory Visualization in PyMOL

#### 3.8.1. Local Workstation Setup
Visualization runs locally on your workstation/laptop:
```bash
# 1. Provision PyMOL environment
conda create -y -n pymol_env -c conda-forge python=3.11 pymol-open-source
conda activate pymol_env

# 2. Install lightweight ColabMDA package
python3 -m pip install --upgrade "git+https://github.com/paulshamrat/ColabMDA.git@main"
```

#### 3.8.2. Launching Trajectory Visualization
```bash
# Visualize Apo WT replica 1
colabmda openmm view --pdb-dir data/simulation/apo/WT/r1

# Visualize Holo L119R replica 1 highlighting mutated Arg119 sidechain
colabmda openmm view --pdb-dir data/simulation/atp_mg/L119R/r1 --resi 119
```

<details>
<summary><b>Click to view the generated PyMOL script configuration (<code>visualize.pml</code>)</b></summary>

```text
# 1. Clear memory and set clean canvas
reinitialize
bg_color white

# 2. Load merged trajectory coordinates
load prod_full.pdb, kinase
load_traj prod_full.dcd, kinase

# 3. Align trajectory to crystal reference
load 4bgq_wt.pdb, crystal_ref
align kinase and name CA, crystal_ref and name CA, mobile_state=1, target_state=1
delete crystal_ref

# 4. Freeze backbone tumbling frame-by-frame
intra_fit kinase and name CA

# 5. Crisp secondary structure ribbon
dss kinase
cartoon automatic, kinase
show cartoon, kinase
color cyan, kinase and name C*
util.cnc("kinase")

# 6. Highlight mutated residue 119 sidechain sticks
show sticks, kinase and resi 119 and not name N+C+O+H
color yellow, kinase and resi 119 and name NH1+NH2
set stick_radius, 0.25

# 7. Focus camera onto kinase domain
zoom kinase and polymer, buffer=4
```
</details>

#### 3.8.3. Comparative Transition Snapshot Grids
Render publication-quality structural snapshot grids comparing conformational transitions across time points:
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
      "pdb": "data/simulation/apo/WT/r1/prod_full.pdb",
      "dcd": "data/simulation/apo/WT/r1/prod_full.dcd",
      "states": [1, 100],
      "times": ["0.00 ns", "1.00 ns"],
      "mut_residue": 119
    }
  }
}
```

---

## 4. Benchmark Pilot Run: 1 ns Equilibration + 1 ns Production (WT & L119R)

Before launching full 100 ns production campaigns, execute a standardized **1 ns Equilibration + 1 ns Production Pilot Run** to benchmark GPU performance, verify active-site coordination, and validate the automated autoimaging and analysis pipelines.

### 4.1. Pilot Run Benchmark Parameters
* **Equilibration Duration:** 1 ns total (`--equil-time 1000` = 1,000 ps across NVT and NPT).
* **Production Duration:** 1 ns total (`--total-ns 1.0` = 500,000 steps at 2.0 fs timestep).
* **Trajectory Frame Interval:** `--traj-interval 10` (10 ps = 100 frames).
* **Expected Walltime:** ~1.5–2.5 minutes per replica on modern NVIDIA GPUs (T4 / V100 / A100 / RTX).

---

### 4.2. Google Colab End-to-End Pilot Execution

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

---

### 4.3. HPC Cluster Scratch End-to-End Pilot Execution

```bash
# ==============================================================================
# HPC Cluster Scratch Pilot Execution (SLURM)
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

---

## 5. Appendix: Backward-Compatibility Mapping (v0.2.0 to v0.3.0)

<details>
<summary><b>Click to expand translation table for legacy <code>colabmda openmm</code> subcommands</b></summary>

For users maintaining legacy shell scripts from earlier versions of ColabMDA, the table below maps legacy multi-keyword commands to their modern single-keyword `cmda` equivalents:

| Legacy Command (v0.2.0) | Modern Command (v0.3.0) | Key Advantages |
| :--- | :--- | :--- |
| `colabmda modeller build ...` | `cmda modeller build ...` | Direct top-level access, biological numbering by default. |
| `colabmda openmm prep ...` | `cmda prep --state <tag> --variant <var>` | Combines PDBFixer, PROPKA 3 at pH 7.4, and Amber tleap solvation into one step. |
| `colabmda openmm stage ...` | `cmda stage --state <tag> --variant <var>` | Auto-creates standard `data/simulation/<state>/<var>` layout and decoupled scripts. |
| `colabmda openmm em/nvt/npt/qc` | `cmda sim --state <tag> --stage eq` | Executes full 24-stage relaxation hub with automated QC gating in one command. |
| `colabmda openmm md ...` | `cmda sim --state <tag> --stage md --replica r1` | Resume-safe checkpointed production MD. |
| `colabmda openmm merge ...` | `cmda merge --state <tag> --variant <var>` | Automated `cpptraj` PBC autoimaging and dual-stream solvent stripping. |
| `colabmda openmm analysis ...` | `cmda analyze --state <tag> --variant <var>` | Generates time-series metrics (RMSD, Rg, RMSF, PCA). |
| `colabmda openmm compare ...` | `cmda compare --state <tag> --variants WT,L119R` | Multi-replica and cross-variant comparative overlays. |
| `colabmda openmm view ...` | `colabmda openmm view --pdb-dir <dir>` | Generates automated `visualize.pml` and launches PyMOL. |

</details>

---

## Acknowledgements
- [OpenMM](https://openmm.org) & [PDBFixer](https://github.com/openmm/pdbfixer)
- [Modeller](https://salilab.org/modeller/)
- [AmberTools](https://ambermd.org/AmberTools.php) (`tleap`, `cpptraj`)
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
*   **Comprehensive Test Suite:** 50+ unit and regression tests passing across all modeling, parameterization, and simulation modules.

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
