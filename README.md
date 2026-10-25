# ColabMDA

**ColabMDA** supports protein modeling, parameterization, checkpointed OpenMM simulations, trajectory processing and analysis on Google Colab, HPC clusters and Linux workstations. This README is the user guide, from platform setup through apo/holo workflows and visualization.


---

## 🛠 Project Information
| Category | Details |
| :--- | :--- |
| **Release** | [![GitHub tag](https://img.shields.io/github/v/tag/paulshamrat/ColabMDA)](https://github.com/paulshamrat/ColabMDA/tags) |
| **Availability** | [![GitHub](https://img.shields.io/badge/GitHub-Repository-blue?logo=github)](https://github.com/paulshamrat/ColabMDA) |
| **Workflows** | [![Python CI](https://github.com/paulshamrat/ColabMDA/actions/workflows/python-test.yml/badge.svg?branch=main)](https://github.com/paulshamrat/ColabMDA/actions/workflows/python-test.yml) |
| **Issues** | [![GitHub issues](https://img.shields.io/github/issues/paulshamrat/ColabMDA)](https://github.com/paulshamrat/ColabMDA/issues) |
| **License** | [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) |
| **Style / Lint** | [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black) [![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff) |
| **Dependencies** | `OpenMM`, `AmberTools (tleap, cpptraj)`, `Modeller`, `MDAnalysis`, `MDTraj`, `PDB2PQR / PROPKA` |
| **Platform** | `Linux` · `Google Colab (CLI / Terminal)` · `HPC (SLURM)` |
| **Structure** | [Source Layout](#project-structure) |

### 📂 Project Structure <a id="project-structure"></a>
ColabMDA is organized into clear, functional modules:
*   `src/colabmda/cli.py`: Central unified CLI (`cmda` and `colabmda`).
*   `src/colabmda/core/`: Project paths, safe example installation, small-protein workflows, state taxonomy and ion definitions.
*   `src/colabmda/amber/`: AmberTools parameterization (`tleap`, Hu 2024 quantum nucleotide parameters, `cpptraj` autoimaging & solvent stripping).
*   `src/colabmda/modeller/`: Homology modeling engine (UniProt biological numbering, mutagenesis, geometric substrate docking).
*   `src/colabmda/openmm/`: Flat operation modules: `prepare.py`, `em.py`, `nvt.py`, `npt.py`, `md.py`, `qc.py`, `merge.py`, `analysis.py`, `figures.py`, `common.py`, `commands.py` and the compatibility adapter `run.py`.
*   `src/colabmda/remote/`: Google Colab GPU remote session management via `colab` CLI.
*   `scripts/`: Single standalone bootstrap script (`scripts/setup.sh`).
*   `examples/{fkbp,hras,trpcage}/scripts/`: Runnable templates copied into the execution project.
*   `examples/hpc/{simulate,analyze}.slurm`: Thin GPU simulation and CPU analysis adapters for the shared scientific commands.
*   `docs/{apo,atpmg}.md`: CDKL5 workflows, with setup under `docs/platforms/{colab,hpc}.md`.
*   `tests/`: Comprehensive unit and regression test suite (`pytest`).

---


CDKL5 command-based workflows are in [Protein in water](docs/apo.md) and
[ATP with one Mg](docs/atpmg.md). Platform instructions are separate:
[Google Colab](docs/platforms/colab.md) and [HPC](docs/platforms/hpc.md).
These guides use `data/{apo,atpmg}/{structure,simulation}` for new runs;
existing numbered examples below retain their historical paths. The full README
baseline remains here alongside the shared guides.

Current CDKL5 checks cover archive mapping, WT/L119R modeling, reviewed ATP/Mg
transfer and actual Amber topology checks for four systems. Reviewed-backbone WT
pilots passed full 24-stage, 1.96 ns equilibration QC on Colab T4 (apo) and Palmetto
A100 (ATP/one-Mg). Production and analysis cover 200 ps/20 frames and 100 ps/10
frames, respectively. Colab also passed an actual interruption/resume check with
unchanged committed chunks. L119R has CPU modeling/preparation checks only.
These short pilots verify the workflow; they do not establish biological convergence.

## 1. Platform Setup

Use the checkout whose instructions you are reading. This development guide includes
unreleased changes on `dev1`; the public `main` branch may contain an older release.
The tool stays separate from the project where simulations run.

### 1.1. Google Colab

#### 1.1.1. Terminal Runner (Recommended)

Prerequisites are Linux/Bash, Git, Python 3.12 or newer for the current Colab CLI,
internet access and a Colab account. The core ColabMDA package supports Python
3.10+, but the local Colab control environment must also meet the CLI requirement.
Complete CLI authentication when requested. Existing installations can activate
their environment instead of creating another one.

**On your local computer**, enter the ColabMDA checkout. In our workspace this is
`~/research/tools/ColabMDA`. A GPU is not required on the local computer:

```bash
cd ~/research/tools/ColabMDA
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e . google-colab-cli
cmda --help
colab sessions
```

Finish local changes and commit the deployment snapshot before creating the GPU
runtime, so it is ready for use. An idle runtime can expire; a replacement needs
fresh Drive authentication even when the previous session was mounted.

Then create a GPU runtime, authenticate Drive and deploy the committed tool:

```bash
colab new -s colabmda-example --gpu T4
colab drivemount --session colabmda-example
cmda colab setup -s colabmda-example --local-source
colab console -s colabmda-example
```

Local-source deployment currently requires a clean, committed `dev1` checkout.
It sends runtime files and example templates with `SNAPSHOT.json` hashes. Run
these commands serially and wait for setup to finish before entering the console.
Setup asks whether to install MODELLER; the small example does not require it.

**Inside the Colab console**, activate the installed environment:

```bash
source /root/miniforge3/etc/profile.d/conda.sh
conda activate colabmda
export MPLBACKEND=Agg
cmda --help
```

The installed snapshot is `/content/ColabMDA`; persistent project data belongs under
`/content/drive/MyDrive/research/projects/`. Section 4 creates and runs the complete
small example there. You can also dispatch foreground shell commands from the local
terminal using `cmda colab exec -s colabmda-example --timeout 21600 'COMMAND'`.

Each remote exec call starts a new shell. Include environment activation, project
variables and `cd` in that call; shell state does not carry between invocations.
A quoted multiline shell block runs in one Bash process and propagates failure.

#### 1.1.2. Web Notebook

With a GPU runtime selected in Colab, mount Drive in a notebook cell:

```python
from google.colab import drive
drive.mount('/content/drive')
```

For a published release, the notebook terminal can download and run that release's
`scripts/setup.sh`. For this private development checkout, use the terminal-runner
local-source deployment above so the installed code matches this guide. Downloading
public `main` does not deploy unpublished local changes.

#### 1.1.3. Crash Recovery & Checkpointing

Keep simulations in the foreground with live output. Sessions can still terminate.
Keep inputs, completed chunks and checkpoints on Drive. On a replacement runtime,
remount Drive and restore the recorded software before attempting recovery.
Checkpoint compatibility across different GPUs or software versions is not assumed.
Section 3.5 explains restart behavior; section 4 provides a short test.

#### 1.1.4. MODELLER License

Model building and mutation require MODELLER and a valid key. Its registration page
is [salilab.org/modeller/registration.html](https://salilab.org/modeller/registration.html).
The installer preserves the optional question and hidden key prompt.

```bash
cmda colab setup -s colabmda-example --local-source --modeller
```

Automated setup accepts `KEY_MODELLER` or `MODELLER_LICENSE` through the process
environment. Do not place real keys in commands saved to Git, logs or provenance.
`--no-modeller` skips the optional remote question. MODELLER uses `modeller_env`;
section 3.2 explains activating it and installing the same tool checkout there.

---

### 1.2. HPC Clusters (SLURM)

#### 1.2.1. Scratch Storage

Install the tool and its Conda dependencies on the cluster according to its software
policy. Put active simulation outputs in scratch, not the software checkout.
Activate the selected MD environment on the allocated compute node, then use the same example scripts
or state-based commands below with an explicit project root.

For Palmetto access from Divergent, use the workspace's `hpc` CLI for inspection,
code deployment and submission. Do not use direct SSH/SCP for Palmetto. Transfer
computational code and inputs only; manuscripts and private editorial material stay local.

#### 1.2.2. Interactive GPU Allocation

On a generic Slurm cluster, a small interactive request might be:

```bash
salloc --nodes=1 --ntasks=1 --cpus-per-task=8 --gpus=1 --mem=32G --time=01:00:00
```

Partition, account, GPU flags and environment activation are cluster-specific.
For Palmetto, obtain the allocation through the approved `hpc` workflow. Validate
the small example before designing a production campaign.

#### 1.2.3. Batch Submission

State-based `cmda stage --scripts` generates `eq.slurm`, `prod.slurm` and `chain.sh`. These are
starting templates for an existing research system, not a scheduler configuration
for every cluster. Review resources, environment, scratch paths, state, variant,
replica bounds, scientific settings and QC before submission.

Stage one state/variant at a time when using those shared script names: later stage
calls replace the shared runner files. The Trp-cage `--example` mode instead copies
the numbered workflow scripts shown in section 2; it does not generate Slurm jobs.

---

### 1.3. Local Workstation

#### 1.3.1. Installation

**On a Linux workstation with an NVIDIA GPU**, enter the tool checkout and run its
single installer with local paths. This ensures it installs the current checkout
and does not try to create a Colab-only `/content` directory:

```bash
cd ~/research/tools/ColabMDA
export COLABMDA_TOOL="$PWD"
COLABMDA_SOURCE="$COLABMDA_TOOL" \
  INSTALL_DIR="$HOME/.local/share/colabmda" \
  WORK_DIR="$HOME/research" \
  bash scripts/setup.sh
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate colabmda
cmda --help
```

With an existing Conda installation, the core environment can also be created
from `environment.yml`, followed by `python -m pip install -e .` in the checkout.
The installer can optionally configure MODELLER and its hidden license prompt.
Divergent has no GPU; use section 1.1 for its control-plane installation and Colab
execution rather than running a protein pilot locally.

#### 1.3.2. Test Suite

Standard tests run on the developer workstation in an environment containing the
scientific dependencies and `pytest`:

```bash
python -m pip install pytest
python -m pytest tests/ -q
```

These CPU tests check software behavior. They do not establish physical convergence
or replace GPU integration testing. Development rules are in [CONTRIBUTING.md](CONTRIBUTING.md)
and [AGENTS.md](AGENTS.md).

---

## 2. Project Architecture

For CDKL5 WT/L119R, use this shared layout on Colab and HPC:

```text
colabmda-example/
├── data/{apo,atpmg}/
│   ├── structure/{raw,mapping,models/{WT,L119R}}/
│   └── simulation/<run>/<variant>/
│       ├── 01_prep/
│       ├── 02_eq/
│       ├── 03_md/{r1,r2,r3}/
│       └── 04_analysis/
└── logs/{apo,atpmg}/<run>/<variant>/
```

Only processing steps are numbered. Use a fresh `COLABMDA_RUN_ID` for changed
inputs or protocol. `cmda stage --state apo --variant WT --root "$PROJECT_ROOT"`
creates the directories and no runner scripts by default. Colab runs direct
commands; the [HPC adapter](examples/hpc/simulate.slurm) runs those same commands
inside Slurm. The adapter supports a separately named `smoke` software check and
requires a GPU; it can select an existing conda environment. Platform setup is in [Colab](docs/platforms/colab.md) and
[HPC](docs/platforms/hpc.md); the state guides hold the scientific workflow.

The numbered example architecture below remains supported for existing runs.


The tool repository contains the implementation and reusable workflow templates.
Each execution project contains its own scripts, inputs, outputs and logs. No
project configuration YAML is required.

To create the complete Trp-cage example on the machine where it will run:

```bash
cmda stage --example trpcage --root "$HOME/research/projects/colabmda-example"
```

On Colab, use `/content/drive/MyDrive/research/projects/colabmda-example` instead.
Staging copies seven small scripts and records their tool version and hashes. It
does not fetch structures, install software or start a simulation. Identical files
are safe to stage again; different existing files stop staging and are preserved.

```text
colabmda-example/
├── scripts/
│   ├── 01_structure/01_fetch.py
│   └── 02_simulation/
│       ├── common.sh
│       ├── 01_prepare.sh
│       ├── 02_equilibrate.sh
│       ├── 03_production.sh
│       ├── 04_analysis.sh
│       └── 05_check.py
├── data/
│   ├── 00_setup/example.json
│   ├── 01_structure/trpcage/           # Raw archive, selected model, input hashes
│   └── 02_simulation/<run>/
│       ├── provenance.json
│       ├── production.json
│       └── apo/WT/
│           ├── 01_prep/
│           ├── 02_eq/
│           ├── 03_md/r1/
│           └── 04_analysis/figures/
└── logs/02_simulation/<run>/
```

Scientific projects can use their existing state/variant layout. `ProjectPaths`
resolves the existing layout; do not create parallel compatibility trees. The
numbered example layout above is the default for newly staged examples. Keep
runtime data and logs outside Git and preserve native raw coordinates and provenance.

The older FKBP12 and H-Ras templates remain available through `cmda stage
--example fkbp` or `--example hras`. They use the numbered example layout.
The CDKL5 commands in section 3 use the state-based layout above.

Independent replicas may share an equilibration endpoint, while each production
replica receives its own velocities, stochastic streams and checkpoint. They are
independent simulation realizations, not proof of converged or uncorrelated sampling.

---

## 3. Simulation Protocol (cmda)

This workflow uses CDKL5 WT and L119R, with two states:

| Track | Protein input | State | Variants |
| :--- | :--- | :--- | :--- |
| Apo | 4BGQ, chain A; UniProt O76039, residues 1–303 | Protein in water (`apo`) | WT, L119R |
| Holo | Same reviewed protein backbone, with ATP and one Mg | ATP/Mg (`atpmg`) | WT, L119R |

Start with WT. After its short validation passes, repeat with L119R. Use matching
protein models and settings when comparing states. Full 24-stage equilibration
and short WT production have been tested on Colab (apo) and HPC (ATP/Mg).
L119R has modeling and preparation checks only; its GPU validation remains to be
run. Short pilots do not establish mutation effects or biological convergence.

Structure sources: [4BGQ](https://www.rcsb.org/structure/4BGQ) and
[O76039](https://www.uniprot.org/uniprotkb/O76039/entry).
The separate [apo](docs/apo.md) and [ATP/Mg](docs/atpmg.md) guides provide the same
workflow by state. Section 4 retains Trp-cage as an optional quick software test.

### 3.1. Command Reference

| Command | Purpose |
| :--- | :--- |
| `cmda fetch` / `check` / `assemble` | Map an archive to UniProt, verify a model or transfer a reviewed ATP/Mg pose. |
| `cmda modeller build` / `mutate` | Model a selected sequence or mutate a supplied structure using MODELLER. |
| `cmda stage --example NAME --root PROJECT` | Install fkbp, hras or trpcage scripts without overwriting user edits. |
| `cmda stage --state TAG --variant NAME --root PROJECT` | Create state/variant folders; add `--scripts` for compatibility runners. |
| `cmda prep --state TAG --variant NAME --root PROJECT` | Assemble reviewed coordinates and parameterize them with AmberTools. |
| `cmda sim --state TAG --variant NAME --stage STAGES --root PROJECT` | Run EM, NVT, NPT, QC or production MD. |
| `cmda merge --state TAG --variant NAME --replicas r1 --root PROJECT` | Autoimage trajectories and write matching stripped topologies. |
| `cmda analyze --state TAG --variant NAME --root PROJECT` | Calculate RMSD, geometric heavy-atom Rg and C-alpha RMSF, then plot the tables. |
| `cmda openmm compare --series LABEL=DIR1,DIR2 --outdir OUT` | Plot comparisons from existing analysis tables. |
| `cmda view --topology FILE --trajectory FILE --resi NUMBER` | Generate a PyMOL view from a matching topology/trajectory pair. |
| `cmda snapshots --config FILE --output FILE` | Render a custom snapshot description through the existing PyMOL helper. |
| `cmda colab setup` / `exec` / `stop` | Set up and control a remote Colab session. |

Use `--help` for each command's arguments. Standard analysis does not calculate
PCA or contact maps. Top-level `cmda compare` currently reports the request but
produces no plots; the implemented plotting command is `cmda openmm compare`.

### 3.2. Homology Modeling & Mutagenesis

#### Biological Numbering

**Inside Colab**, after setup, define the project and run:

```bash
export COLABMDA_TOOL=/content/ColabMDA
export PROJECT_ROOT=/content/drive/MyDrive/research/projects/colabmda-example
export COLABMDA_RUN_ID=cdkl5_reviewed_v1
mkdir -p "$PROJECT_ROOT"
cd "$PROJECT_ROOT"
cmda fetch --state apo --pdb-id 4BGQ --uniprot-id O76039 --chain A \
  --range 1 303 --mutations L119R --root "$PROJECT_ROOT"
cd "$PROJECT_ROOT/data/apo/structure"
```

On HPC, use the scratch project path and installed tool path from the
[platform guide](docs/platforms/hpc.md). Export these variables again in each new
shell. Use a fresh run ID when inputs or settings change.

Fetch keeps the raw mmCIF archive, with PDB fallback, and the UniProt sequence.
`mapping/residues.json` links author positions and insertion codes to UniProt.
For 4BGQ A, label positions 2–304 map to UniProt 1–303: 276 residues have
coordinates and 27 are missing. The L119 WT identity must agree in UniProt and
the mapped coordinates, with a C-alpha present. Raw coordinates stay unchanged;
derived models use UniProt numbering.

For the reviewed study model, supply `raw/receptor.pdb` from the study owner,
record its source and SHA-256, then verify it:

```bash
mkdir -p models/WT
cp -n raw/receptor.pdb models/WT/WT.pdb
cmda check --state apo --variant WT --root "$PROJECT_ROOT"
```

Our source is the resolved file at
`cdkl5-activation/input/varmodel/01_structure/target.B99990001_with_cryst.pdb`.
It is a study input, not bundled with ColabMDA. Use this same reviewed backbone
for both apo and ATP/Mg.

To build a new WT model instead, use an empty `models/WT` folder and a new run ID:

```bash
conda activate modeller_env
python -m pip install -e "$COLABMDA_TOOL"
cd "$PROJECT_ROOT/data/apo/structure"
cmda modeller build --pdb-id 4BGQ --uniprot-id O76039 --chain A \
  --range 1 303 --uniprot-numbering --template raw/template.pdb \
  --sequence raw/uniprot.fasta --outdir models/WT --seed 13031
conda activate colabmda
cmda check --state apo --variant WT --root "$PROJECT_ROOT"
```

MODELLER is optional during installation and uses a hidden license-key prompt.
Building or mutating a model requires it. A new model needs loop and pocket
review before ATP/Mg placement; the tested fresh rebuild had a severe ATP clash
and was rejected for holo use.

Make L119R from the selected WT model:

```bash
conda activate modeller_env
python -m pip install -e "$COLABMDA_TOOL"
cd "$PROJECT_ROOT/data/apo/structure"
cmda modeller mutate --pdb-in models/WT/WT.pdb --chain A --mut L119R \
  --outdir-mut models/L119R --seed 13031
conda activate colabmda
cmda check --state apo --variant L119R --root "$PROJECT_ROOT"
```

`cmda check` verifies the full sequence, numbering and C-alpha coordinates. It
records hashes in `model_check.json` before writing `MODEL_DONE`.

#### Quality Control (QC)

Review alignment, modeled gaps, stereochemistry and model scores. A sequence
match does not validate loop conformation. Review protonation, termini and
cofactor geometry before preparation.

### 3.3. Track 1: Apo Dynamics

#### 3.3.1. Scope & Inputs

Use the verified CDKL5 protein model without ATP or bound Mg. Preparation adds
water and bulk salt. Start with WT; repeat with L119R after WT validation.

#### 3.3.2. Directory Staging

```bash
conda activate colabmda
export STATE=apo
export VARIANT=WT
cd "$PROJECT_ROOT"
cmda stage --state "$STATE" --variant "$VARIANT" --replicas r1,r2,r3 --root "$PROJECT_ROOT"
```

The run ID groups variants in separate folders. Direct staging creates no runner
scripts unless `--scripts` is requested.

#### 3.3.3. Preparation & Protonation (cmda prep)

```bash
cd "$PROJECT_ROOT/data/$STATE/simulation/$COLABMDA_RUN_ID/$VARIANT/01_prep"
cmda prep --state "$STATE" --variant "$VARIANT" --protein cdkl5 \
  --receptor "$PROJECT_ROOT/data/$STATE/structure/models/$VARIANT/$VARIANT.pdb" \
  --root "$PROJECT_ROOT"
```

Defaults are ff19SB, OPC, 14 Å padding, neutralization and an estimated 0.15 M
NaCl pair count. Inspect the actual box, ion counts, protonation and charge.
This route does not run automatic pKa titration. The lower-level `cmda openmm prep`
supports PDBFixer/PDB2PQR/PROPKA; retain its audit if used.

Inspect `model.pdb`, `solvated.pdb`, `system.prmtop`, `system.inpcrd` and `prep.json`.
`topology_check.json` checks charge, composition, ion parameters and output hashes
before `PREP_DONE`. `preparation_input.json` freezes input hashes. Unchanged
completed preparation is reused; changed inputs require a new run ID.

For L119R, set `VARIANT=L119R` and repeat staging and preparation. Each variant
keeps its own preparation, equilibration and production replicas.

### 3.4. Track 2: Holo Dynamics (ATP + Mg²⁺)

#### 3.4.1. Scope & Chemistry

Use the same reviewed WT and L119R backbones with one ATP and one Mg. The `atpmg`
state adds no phosphorylation or substrate peptide. Fetch its mapping and copy
the verified protein models from the apo structure folder:

```bash
cd "$PROJECT_ROOT"
cmda fetch --state atpmg --pdb-id 4BGQ --uniprot-id O76039 --chain A \
  --range 1 303 --mutations L119R --root "$PROJECT_ROOT"
cd "$PROJECT_ROOT/data/atpmg/structure"
for variant in WT L119R; do
  mkdir -p "models/$variant"
  cp -n "$PROJECT_ROOT/data/apo/structure/models/$variant/$variant.pdb" "models/$variant/$variant.pdb"
  cmda check --state atpmg --variant "$variant" --root "$PROJECT_ROOT"
done
```

#### 3.4.2. ATP & Mg²⁺ Placement

Supply three reviewed files in `data/atpmg/structure/raw/`:

| File | Contents |
| :--- | :--- |
| `reference.pdb` | Protein in the reviewed ATP/Mg coordinate frame |
| `ATP.pdb` | One ATP in that frame |
| `MG.pdb` | One Mg in that frame |

These are study inputs supplied by the owner. Record source paths and SHA-256
in `raw/provenance.json`. The [ATP/Mg guide](docs/atpmg.md#32-supply-the-reviewed-study-inputs)
records the original source locations. The activation study transferred the pose
from 8FP5 using kinase-core C-alpha residues 30–220; this workflow reuses that
reviewed pose.

```bash
export STATE=atpmg
export VARIANT=WT
cd "$PROJECT_ROOT/data/atpmg/structure"
cmda assemble --receptor "models/$VARIANT/$VARIANT.pdb" --reference raw/reference.pdb \
  --ligand raw/ATP.pdb --ions raw/MG.pdb --mapping mapping/residues.json \
  --variant "$VARIANT" --core 30 220 --max-rmsd 1.0 \
  --out "models/$VARIANT/complex.pdb"
```

Assembly fits resolved core C-alpha atoms, excluding modeled gaps. It checks
identities, the declared substitution, component counts and severe protein/ATP
clashes. `complex.json` records fit residues, RMSD, short contacts and hashes.
Inspect the pocket and Mg coordination; passing a clash gate is not a complete
geometry assessment. Repeat with `VARIANT=L119R` after WT validation.

#### 3.4.3. Directory Staging & Holo Parameterization

```bash
export STATE=atpmg
export VARIANT=WT
cd "$PROJECT_ROOT"
cmda stage --state "$STATE" --variant "$VARIANT" --replicas r1,r2,r3 --root "$PROJECT_ROOT"
cd "$PROJECT_ROOT/data/$STATE/simulation/$COLABMDA_RUN_ID/$VARIANT/01_prep"
cmda prep --state "$STATE" --variant "$VARIANT" --protein cdkl5 \
  --receptor "$PROJECT_ROOT/data/$STATE/structure/models/$VARIANT/complex.pdb" \
  --root "$PROJECT_ROOT"
```

Supply the assembled complex directly; do not also pass `--template`, which could
reposition its cofactors. Preparation uses the bundled Hu et al. B3 ATP parameters
with ff19SB/OPC and complete water-matched Li–Merz 12–6 ions. Check the actual
topology: ATP has 43 atoms and charge −4; the single Mg has charge +2.
`topology_check.json` also checks neutral total charge and Na/Cl/Mg LJ values.
Loading 12–6–4 LJ values without the extra C4 term is incomplete; this route
uses the complete 12–6 model. See the
[ParmEd C4 documentation](https://parmed.github.io/ParmEd/html/parmed.html#add12-6-4).

Repeat staging and preparation with `VARIANT=L119R` after its assembly is reviewed.

#### 3.4.4. Additional State Conventions

ATP, ADP and GTP resources remain available for appropriate structures. Choose
state tags and ion counts for the actual experimental/model composition. Existing
phosphorylation and substrate-peptide conventions remain advanced capabilities;
they are not required by the small tutorials or validated by their acceptance tests.
For kinase projects, the previously documented composition tags remain available:

| State tag | Requested composition |
| :--- | :--- |
| `adp_2mg_py171` | ADP, two Mg ions and phosphorylated Tyr171. |
| `atp_2mg_py171_peb2` | ATP, two Mg ions, phosphorylated Tyr171 and the EB2 peptide preset. |

These tags request components; review phosphorylation parameters, peptide
placement and coordination before interpreting an activation state.
Generic GAFF parameters are not equivalent to the specialized nucleotide route.
Preserve a matching Amber topology/coordinate pair for independently parameterized
systems, and review protonation, charge, ligand atom names and coordination.

### 3.5. Equilibration, QC & Production

#### 3.5.1. Staged Equilibration

Choose the prepared state and variant in the current shell. For example, to run
apo WT after section 3.3:

```bash
export STATE=apo
export VARIANT=WT
cd "$PROJECT_ROOT/data/$STATE/simulation/$COLABMDA_RUN_ID/$VARIANT/02_eq"
cmda sim --state "$STATE" --variant "$VARIANT" --stage eq --protocol staged \
  --root "$PROJECT_ROOT"
```

Use `STATE=atpmg` for the prepared ATP/Mg system. Repeat with `VARIANT=L119R`
after WT validation. Keep the same state, variant and run ID for later commands.
Do not set `COLABMDA_TEST_1_24` for research. Quick/reduced runs require
`--software-test`; the optional Trp-cage example in section 4 uses that route.

Restraints act on solute heavy atoms; bulk water and monovalent ions remain
unrestrained. Strengths below are in kcal mol⁻¹ Å⁻²:

| Stage | Operation | Conditions | Restraint |
| :--- | :--- | :--- | :--- |
| EM | Initial minimization | No physical time | None |
| 01–03 | Minimization ladder | Constant box | 100, 50, 25 |
| 04–09 | Heating ladder | Six 10 ps NVT intervals, 50–300 K, 1 fs | 25 |
| 10–14 | Minimization ladder | Constant box | 25, 10, 5, 2, 1 |
| 15–18 | NVT release ladder | Four 100 ps intervals, 300 K, 2 fs | 10, 5, 2, 1 |
| 19–22 | NPT release ladder | Four 100 ps intervals, 300 K, 1 bar, 2 fs | 10, 5, 2, 1 |
| 23 | Unrestrained NPT | 100 ps, 300 K, 1 bar, 2 fs | None |
| 24 | Final unrestrained NPT | 1000 ps, 300 K, 1 bar, 2 fs | None |

Numbered state/JSON records support recovery from completed stages. The quick
profile in section 4 instead uses one short NVT and NPT interval; it is not this
full protocol. The legacy combined `openmm run` adapter retains ff14SB/TIP3P defaults
and delegates to quick equilibration; it should not be confused with the staged route.

#### 3.5.2. Physical Settings & QC

Shared dynamics use PME with a 1.0 nm cutoff, switching from 0.8 nm, hydrogen-bond
constraints and rigid water. LangevinMiddle friction is 1/ps. The isotropic NPT
barostat targets 1 bar and attempts a move every 25 steps. Inspect the saved
parameterization record rather than assuming every preparation route is identical.

The native OpenMM construction route offers ff19SB/ff14SB and installed OPC,
TIP3P, TIP3P-FB or TIP4P-Ew XML files. Auto padding selects 1.0, 1.2 or 1.4 nm based
on solute extent; it does not replace a supplied Amber topology or the high-level
LEaP defaults. That native route neutralizes without extra bulk salt by default.
Its four-site OPC construction uses TIP4P-Ew starting geometry with OPC parameters.
Retain `prep.json`, `parameterization.json`, the topology, `system.xml`,
`npt.state.xml` and the run/input/runtime records with the production outputs.

Read `equilibration_qc.json` and the downstream PNGs under `figures/`. Research `cmda sim` enforces staged timing and physical QC before production.
Explicit `--software-test` runs report physical QC without a research gate. A successful exit or
an output marker does not certify physical equilibration. The separate
`cmda openmm qc` command fails on QC failure unless `--warn-only` is selected.
No QC PDF is produced.

#### 3.5.3. Production MD (r1, r2, r3)

After full equilibration passes QC, run a short WT production check:

```bash
cd "$PROJECT_ROOT/data/$STATE/simulation/$COLABMDA_RUN_ID/$VARIANT/03_md"
cmda sim --state "$STATE" --variant "$VARIANT" --stage md --replica r1 \
  --protocol staged --total-ns 0.1 --checkpoint-ps 100 --traj-interval 10 \
  --root "$PROJECT_ROOT"
```

This produces ten frames over 100 ps. Review analysis and cofactor geometry before
running L119R. For a scientific campaign, choose a justified duration and replica
count. The following illustrates three 100 ns targets:

```bash
cmda sim --state "$STATE" --variant "$VARIANT" --stage md --replica r1,r2,r3 \
  --total-ns 100 --checkpoint-ps 1000 --traj-interval 10 --root "$PROJECT_ROOT"
```

Repeat for each reviewed state/variant and its run ID. Targets are per replica.
The high-level route processes the listed replicas sequentially; parallel batch
execution needs an appropriate scheduler design.

#### 3.5.4. Timestep, Steps & Saved Output

The **integration timestep**, **frame-saving interval** and **checkpoint interval**
are different settings:

| Quantity | Current implementation | Meaning |
| :--- | :--- | :--- |
| Integration timestep | 1 fs during staged heating; 2 fs during subsequent dynamics, quick equilibration and production | Physical time advanced by one integration step. |
| Frame interval | `--traj-interval`, in ps; example uses 10 ps | How often coordinates and production measurements are saved. |
| Checkpoint/chunk interval | `--checkpoint-ps`, in ps; example uses 100 ps | Completed interval published before the production checkpoint advances. |
| Production target | `--total-ns`, in ns, per replica | Total requested production time, excluding equilibration. |

Unit conversions are **1 ps = 1,000 fs** and **1 ns = 1,000 ps = 1,000,000 fs**.
At the production timestep of 2 fs = 0.002 ps:

$$N_\mathrm{steps} = \frac{T_\mathrm{ns}\times 10^6}{\Delta t_\mathrm{fs}}$$

$$N_\mathrm{steps/frame} = \frac{\Delta t_\mathrm{frame,ps}\times 10^3}{\Delta t_\mathrm{fs}},\qquad
N_\mathrm{steps/chunk} = \frac{\Delta t_\mathrm{chunk,ps}\times 10^3}{\Delta t_\mathrm{fs}}$$

| Duration or interval | Calculation at 2 fs | Integration steps |
| :--- | :--- | ---: |
| 1 ps | 1,000 / 2 | 500 |
| 10 ps frame interval | 10,000 / 2 | 5,000 |
| 100 ps equilibration interval or chunk | 100,000 / 2 | 50,000 |
| 1 ns production | 1,000,000 / 2 | 500,000 |
| 100 ns production | 100,000,000 / 2 | 50,000,000 |

A 1 ns run with 10 ps frames and 100 ps chunks therefore has **100 frames**,
**10 completed chunks** and **10 frames per chunk**. Its first saved production
frame is at 10 ps, not at 0 ps. A 10 ps heating interval at 1 fs instead requires
10,000 steps. Minimization iterations do not represent physical simulation time.

The normal staged schedule integrates 60 ps heating + 400 ps NVT + 1,500 ps NPT
= **1,960 ps (1.96 ns)** before production. The quick example uses 100 ps NVT +
100 ps NPT = **200 ps (0.2 ns)** before its 1 ns production. Equilibration is not
part of the production target. Three 100 ns replicas provide 300 ns aggregate
production, not one continuous 300 ns trajectory.

Production currently fixes the timestep at 2 fs; there is no production timestep
CLI option. Changing `--traj-interval` changes saved output, not the integrator.
The high-level `--steps` testing override changes the requested production step
count as well as quick equilibration. Use `--total-ns` for research duration. Use intervals corresponding to whole integration steps;
the implementation converts requested intervals to integer step counts.

#### 3.5.5. Interruption & Checkpoint Recovery

`--total-ns` is an absolute production target, excluding equilibration. Extending
0.1 ns to 1 ns requests 1 ns total, not another 1 ns. New replicas start from the
portable NPT state with fresh velocities and stochastic streams, resetting time
to zero. They do not branch from a sibling's binary checkpoint.

A complete trajectory/log chunk is published before advancing `prod.chk`.
After Ctrl+C or a disconnection, rerun the same command with compatible inputs,
settings, hardware and software. The unfinished chunk is replayed from the last
committed checkpoint. Verify unchanged committed chunks and continuous times.
Portable `*.state.xml` states transfer physical data between ensembles; binary
checkpoints contain Context/RNG internals and are not assumed portable.

### 3.6. Autoimaging & Solvent Stripping (cmda merge)

#### 3.6.1. PBC Autoimaging

```bash
cd "$PROJECT_ROOT/data/$STATE/simulation/$COLABMDA_RUN_ID/$VARIANT/03_md"
cmda merge --state "$STATE" --variant "$VARIANT" --replicas r1 --root "$PROJECT_ROOT"
```

Use `--replicas r1,r2,r3` after all three replicas have finished. Repeat for each
simulated variant. `cpptraj` autoimages the completed chunks and
writes both a solute stream and a C-alpha stream. Storage savings depend on system
composition; inspect the results rather than assuming a fixed reduction.

#### 3.6.2. Reference Topologies

| Output | Matching topology | Use |
| :--- | :--- | :--- |
| `prod_*to*ps.dcd` | Original solvated topology | Raw chunk outputs, including solvent. |
| `prod_keep_ligand.dcd` | `stripped.prmtop` | Solute and retained cofactors/ions. |
| `prod_ca.dcd` | `ca.prmtop` | C-alpha-only trajectory. |
| `prod_full.log` | Not a topology | Merged production times and thermodynamic measurements. |

Never pair a stripped trajectory with the unstripped solvated PDB. Confirm atom
counts and ordering. A viewer requiring PDB needs an exported, matching topology;
verify its residue-label mapping rather than treating Amber indices as author numbers.

### 3.7. Trajectory Analysis & Comparison

#### 3.7.1. System Metrics (cmda analyze)

```bash
cd "$PROJECT_ROOT/data/$STATE/simulation/$COLABMDA_RUN_ID/$VARIANT/04_analysis"
cmda analyze --state "$STATE" --variant "$VARIANT" --root "$PROJECT_ROOT"
```

The high-level route currently analyzes r1. RMSD uses backbone fitting to the first
frame. Reported Rg is the unweighted geometric radius of protein heavy atoms.
RMSF uses aligned C-alpha coordinates, with reference identity/count gates to retain
verified reference labels and insertion codes. CDKL5 models use UniProt positions
1–303; `structure/mapping/residues.json` links them to author numbering. `analysis.json` records the reference
and implementation fingerprints. Numerical tables and figure rendering are separate
modules; figures use 7-inch width, 300 DPI and the publication typography defaults.

Direct `cmda openmm analysis` with an Amber topology requires `--reference` pointing
to the reviewed, numbered model. Its `--interval` argument supports an explicit frame
interval. The high-level route currently uses 10 ps; do not silently use it for a
trajectory sampled at another interval. PCA and contact maps are not calculated by
the standard analysis route.

#### 3.7.2. Comparative Overlays

Use directories containing already calculated `rmsd.csv`, `rg.csv` and `rmsf.csv`:

```bash
cd "$PROJECT_ROOT/data/$STATE/simulation/$COLABMDA_RUN_ID"
cmda openmm compare --series WT="$PWD/WT/04_analysis" \
  --series L119R="$PWD/L119R/04_analysis" --outdir "$PWD/compare"
```

The output PNGs are under the requested directory's `figures/`. Comma-separated
analysis directories within one series group compatible replicas. Check protocol,
residue identity and time axes before aggregation. Top-level `cmda compare` currently
reports a request without plotting, so it is not substituted here.

#### 3.7.3. Frame Counts, Sampling & Stride

With regular sampling and no extra initial frame:

$$N_\mathrm{frames}=\left\lfloor\frac{1000\,T_\mathrm{ns}}{\Delta t_\mathrm{frame,ps}}\right\rfloor$$

For durations divisible by the frame interval, 1 ns at 10 ps gives 100 frames;
100 ns at 10 ps gives 10,000. Saved frame index 0 corresponds to the first saved
frame, not necessarily physical time zero:

$$t_k=t_\mathrm{first}+k\,\Delta t_\mathrm{frame},\qquad k=0,1,\ldots,N-1$$

The lower-level `cmda openmm merge --stride s` can retain every s-th frame. Its
current implementation retains the first frame, giving:

$$N_\mathrm{kept}=\left\lceil\frac{N_\mathrm{frames}}{s}\right\rceil,\qquad
\Delta t_\mathrm{kept}=s\,\Delta t_\mathrm{frame}$$

Thus 10,000 frames with stride 10 give 1,000 frames. If the original sampling is
10 ps, the retained spacing is 100 ps. For a 1 ns trajectory starting at 10 ps,
stride 10 keeps times 10, 110, …, 910 ps; it does not reset the first frame to 0
or guarantee that the last retained frame reaches the original endpoint.

Keep the full-rate trajectory for the standard analysis and use a separately
strided copy for lighter viewing. High-level `cmda analyze` currently assumes
10 ps sampling, and lower-level analysis reconstructs time from an interval.
Before quantitative analysis of strided data, check the actual logged first-frame
time and spacing; changing an interval alone does not restore a lost time origin.
Striding reduces saved observations, not the integration time already simulated.

#### 3.7.4. Storage Calculation

For a DCD containing three 32-bit coordinate values per atom, the coordinate-data
size is approximately:

$$S_\mathrm{bytes}\approx 12\,N_\mathrm{atoms}\,N_\mathrm{frames}$$

For 35,000 atoms and 100 ns production:

| Saving interval | Frames | Approximate coordinate data per replica |
| :--- | ---: | ---: |
| 1 ps | 100,000 | 42 GB |
| 10 ps | 10,000 | 4.2 GB |

These are decimal GB estimates before file headers, box records, checkpoints,
logs or merged/stripped copies. Three replicas at 10 ps need about 12.6 GB for
these raw coordinates alone. Multiply by the number of states and variants as
well, and use the actual solvated atom count; proteins of similar residue count
can have very different solvent boxes. Coarser saving does not guarantee that
all dynamics of interest remain resolved.

#### 3.7.5. Frame & Coordinate Verification

Run this in the MD environment after merging a selected research system. It uses
the matching stripped topology, reports the saved log times, and checks frame
counts and finite coordinates without assuming that a filename proves completion:

```bash
python - <<'PY'
import os
import MDAnalysis as mda
import numpy as np
import pandas as pd
from colabmda.core.project import ProjectPaths
from colabmda.core.taxonomy import SystemSpec

state, variant = os.environ['STATE'], os.environ['VARIANT']
paths = ProjectPaths(os.environ['PROJECT_ROOT'], state=state)
prod = paths.replica_dir(SystemSpec.from_tag(state, variant=variant), 'r1')
u = mda.Universe(str(prod / 'stripped.prmtop'), str(prod / 'prod_keep_ligand.dcd'))
log = pd.read_csv(prod / 'prod_full.log')
assert len(u.trajectory) == len(log) > 0
assert all(np.isfinite(ts.positions).all() for ts in u.trajectory)
print('Atoms:', u.atoms.n_atoms, 'Frames:', len(u.trajectory))
print('First/last saved time (ps):', log['Time (ps)'].iloc[[0, -1]].tolist())
PY
```

Choose the state/variant/replica you actually ran. Check expected times, sampling
and native residue labels against the run settings. The final small-example check
also enforces time continuity and its expected 20 C-alpha atoms.

For CDKL5 compact downloads, follow the [apo](docs/apo.md#7-retrieve-compact-results)
or [ATP/Mg](docs/atpmg.md#7-retrieve-compact-results) collection commands. They use
the same state/run/variant paths as this workflow and verify SHA-256 after download.
Heavy trajectories stay on Drive or HPC scratch. Section 4.4 collects Trp-cage
outputs from its separate example layout.

### 3.8. PyMOL Visualization

#### 3.8.1. Workstation Setup

Use a separate PyMOL environment and install the same tool checkout there. Enter
the checkout first, then run locally:

```bash
conda create -y -n pymol_env -c conda-forge python=3.11 pymol-open-source pillow numpy
conda activate pymol_env
python -m pip install -e .
```

#### 3.8.2. Launching Visualization

Supply a reviewed topology/trajectory pair with matching atom order and residue
labels. The helper accepts explicit files instead of relying on old filename guesses:

```bash
: "${VIEW_TOPOLOGY:?Set the matching viewer topology path}"
: "${VIEW_TRAJECTORY:?Set the matching trajectory path}"
cmda view --topology "$VIEW_TOPOLOGY" --trajectory "$VIEW_TRAJECTORY" --resi 119
```

Inspect the generated PyMOL selection and camera. The helper retains KRAS-oriented
styling; do not assume its default residue selections annotate CDKL5 correctly. Keep native numbering consistent when highlighting a mutation.

#### 3.8.3. Snapshot Grids

```bash
: "${SNAPSHOT_CONFIG:?Set the reviewed snapshot JSON path}"
cmda snapshots --config "$SNAPSHOT_CONFIG" \
  --output "$PROJECT_ROOT/data/analysis/figures/comparison_grid.png" \
  --temp-dir "$PROJECT_ROOT/data/analysis/rendering"
```

This optional rendering description is not a project configuration required for
simulation. It contains an 18-value PyMOL `camera_view`, optional `align_ref_pdb`,
`stable_core_sel`, and a `systems` mapping. Each system supplies matching `pdb`/`dcd`
paths, `states`, `times` and optional `mut_residue`. Frame indices and physical time
are different; derive labels from the actual sampling. Use the same reviewed camera
and homologous alignment selection for comparative panels. The existing renderer
still contains KRAS-specific styling, including residues 57–75; review it before
using it for another protein.

---

## 4. Small Example & Validation

### 4.1. Create and Fetch Trp-cage

The 20-residue Trp-cage example tests software behavior quickly without replacing
the full research workflows above. **On Colab after section 1.1 setup**, run:

```bash
cmda stage --example trpcage \
  --root /content/drive/MyDrive/research/projects/colabmda-example
cd /content/drive/MyDrive/research/projects/colabmda-example
python scripts/01_structure/01_fetch.py
```

On a local GPU workstation use `~/research/projects/colabmda-example` instead.
Fetch prefers mmCIF, with PDB fallback, selecting NMR model 1, author chain A of
1L2Y. It verifies the deposited sequence, all 20 C-alpha atoms and native numbering,
then records the archive/selection hashes. As a designed peptide, it has no canonical
UniProt mutation annotation. No MODELLER key or cofactor template is needed.

### 4.2. Run the Scripts in Order

With `colabmda` active in the generated project:

```bash
export COLABMDA_RUN_ID=run01
bash scripts/02_simulation/01_prepare.sh
bash scripts/02_simulation/02_equilibrate.sh
bash scripts/02_simulation/03_production.sh
bash scripts/02_simulation/04_analysis.sh
python scripts/02_simulation/05_check.py
```

The defaults are 100 ps NVT, 100 ps NPT, 1 ns production, 2 fs steps, 10 ps sampling
and 100 ps checkpoint chunks. To choose another duration without editing scripts,
set `COLABMDA_EXAMPLE_NS` for production; the final check reads the saved request.
Use positive multiples of 0.1 ns. Choose a fresh run ID for changed inputs,
executed scripts or runtime revisions.

To test extension instead of the direct production command above, use a fresh run
ID, complete preparation/equilibration, then run:

```bash
COLABMDA_EXAMPLE_NS=0.1 bash scripts/02_simulation/03_production.sh
COLABMDA_EXAMPLE_NS=1 bash scripts/02_simulation/03_production.sh
```

Then perform analysis and the final check. To include an actual interruption in
this test, pause after the 0.1 ns smoke command, save its first chunk's checksum
from the project root, then start the extension:

```bash
sha256sum "data/02_simulation/$COLABMDA_RUN_ID/apo/WT/03_md/r1/prod_0to100ps.dcd" \
  > "data/02_simulation/$COLABMDA_RUN_ID/first_chunk.sha256"
COLABMDA_EXAMPLE_NS=1 bash scripts/02_simulation/03_production.sh
```

Press Ctrl+C while the second chunk is running. Before another command, `echo $?`
should report a nonzero exit status. Resume with the same production command, then
verify the first committed chunk:

```bash
COLABMDA_EXAMPLE_NS=1 bash scripts/02_simulation/03_production.sh
sha256sum -c "data/02_simulation/$COLABMDA_RUN_ID/first_chunk.sha256"
```

The checksum must report `OK`. Continue with analysis and the final check.
The private working project's `tests/gpu.py` automates that acceptance test using
these same generated scripts; it is not an extra file installed for every user.

### 4.3. Expected Outputs & Checks

After success, `data/02_simulation/run01/` contains the input/tool/script/runtime
provenance, production requests, analysis results, `summary.json` and `RUN_DONE`.
The default run has ten chunks and 100 merged frames, with times 10–1000 ps.
The final check requires finite coordinates/logs/tables, 20 C-alpha atoms, native
RMSF labels 1–20 and nonempty PNGs. Physical equilibration QC is reported separately.

Scripts create the log directory; save terminal output when desired. For example,
from the generated project after preparation/equilibration:

```bash
set -o pipefail
bash scripts/02_simulation/03_production.sh 2>&1 \
  | tee "logs/02_simulation/$COLABMDA_RUN_ID/production.log"
```

To rerender analysis figures from saved tables without rerunning MD:

```bash
python - <<'PY'
import os
from pathlib import Path
from colabmda.openmm.figures import plot_analysis
run = Path('data/02_simulation') / os.environ['COLABMDA_RUN_ID']
plot_analysis(run / 'apo/WT/04_analysis')
PY
```

### 4.4. Runtime, Limits & Cleanup

Our earlier T4 test took approximately 7.5 minutes for production plus its restart
test. That is an observed runtime, not a guarantee. Quick equilibration can fail
the density-slope criterion; retain the warning and do not call this a converged
production study, a mutant/cofactor validation or a full staged-protocol test.

Keep heavy trajectories and checkpoints on Drive. To return from the interactive
console to your local terminal, detach with **Ctrl-B, then D**. This does not stop
the runtime. With the local ColabMDA environment active, download and verify the
compact reports below. Set RUN_ID, STATE and VARIANT to the completed system.
The commands below retrieve the Trp-cage run. CDKL5 uses the state-based paths
in the [apo](docs/apo.md#7-retrieve-compact-results) and
[ATP/Mg](docs/atpmg.md#7-retrieve-compact-results) collection instructions:

```bash
RUN_ID=run01
STATE=apo
VARIANT=WT
SYSTEM="$STATE/$VARIANT"
REMOTE_RUN="/content/drive/MyDrive/research/projects/colabmda-example/data/02_simulation/$RUN_ID"
LOCAL_RUN="$HOME/research/projects/colabmda-example/data/02_simulation/$RUN_ID"
files=(
  summary.json provenance.json production.json RUN_DONE
  "$SYSTEM/01_prep/prep.json"
  "$SYSTEM/02_eq/equilibration_qc.json" "$SYSTEM/02_eq/figures/equilibration_qc.png"
  "$SYSTEM/02_eq/nvt.log" "$SYSTEM/02_eq/npt.log" "$SYSTEM/03_md/r1/prod_full.log"
  "$SYSTEM/04_analysis/analysis.json"
  "$SYSTEM/04_analysis/rmsd.csv" "$SYSTEM/04_analysis/rg.csv" "$SYSTEM/04_analysis/rmsf.csv"
  "$SYSTEM/04_analysis/figures/rmsd_vs_time.png"
  "$SYSTEM/04_analysis/figures/rg_vs_time.png"
  "$SYSTEM/04_analysis/figures/rmsf_per_residue.png"
)
cmda colab exec -s colabmda-example --timeout 120 \
  "cd '$REMOTE_RUN' && sha256sum ${files[*]} > reports.sha256"
for file in "${files[@]}" reports.sha256; do
  mkdir -p "$LOCAL_RUN/$(dirname "$file")"
  colab download -s colabmda-example "$REMOTE_RUN/$file" "$LOCAL_RUN/$file"
done
(cd "$LOCAL_RUN" && sha256sum -c reports.sha256)
```

Require `OK` for every downloaded report, review the PNGs and read the physical
QC result. Then stop the dedicated session **from the local terminal**:

```bash
colab stop -s colabmda-example
```

A restarted session requires Drive authentication and the recorded software again.
Historical run records remain in the working project's README; user instructions
stay here. CPU unit tests, script/CLI checks and a successful installation are
reported separately from an actual GPU integration run.

---

## Acknowledgements
- [OpenMM](https://openmm.org) & [PDBFixer](https://github.com/openmm/pdbfixer); [OpenMM 8 methods paper](https://doi.org/10.1021/acs.jpcb.3c06662)
- [Modeller](https://salilab.org/modeller/)
- [AmberTools](https://ambermd.org/AmberTools.php) (`tleap`, `cpptraj`)
- [MDAnalysis](https://www.mdanalysis.org) & [MDTraj](https://www.mdtraj.org)
- [PDB2PQR & PROPKA](https://server.poissonboltzmann.org/)
- [NumPy](https://numpy.org), [Matplotlib](https://matplotlib.org), [Biopython](https://biopython.org)
- [Google Colab](https://colab.research.google.com) & [Miniforge/Conda](https://github.com/conda-forge/miniforge)

---

## Changelog

### Unreleased Development

- Use CDKL5 WT/L119R for the root apo and ATP/one-Mg walkthrough, keeping the
  existing README design, full protocol and calculations. Trp-cage remains optional.

- Add CDKL5 apo and ATP/one-Mg command guides with shared Colab/HPC platform setup
  and unnumbered state/structure/simulation folders for new runs.
- Check archive-to-UniProt mapping, model identities, reviewed cofactor transfer,
  topology composition/charges and preparation reuse before simulation.
- Exclude Amber bulk monovalent ions from solute restraints, preserve cumulative
  stage time on ensemble transfer/recovery, and enforce research QC before MD.

- Added small FKBP12 WT/F37V and H-Ras–GDP/Mg WT/G12V workflows, replacing CDKL5
  as the main guide illustration while retaining the full protocol reference.
- Added explicit author-to-UniProt mapping, reviewed local modeling inputs and
  strict sequence/mutation checks. Raw archives and historical Trp-cage results
  remain unchanged. CPU modeling/preparation and all four short T4 workflows
  passed, including actual WT interruption/resume and reviewed analysis plots.

- Fixed ion preparation to retain the complete water-matched 12–6 potential;
  earlier diagnostic runs that overrode it with 12–6–4 LJ values have no C4 term.
  Preserve those original records and use new run IDs after the correction.

* Flat OpenMM modules and one checkpointed production implementation.
* Complete, safely staged Trp-cage example and a single installer with optional MODELLER prompting.
* Native residue-label checks, runtime/input fingerprints and separate numerical/figure modules.
* One example project name, with preserved historical results and migration records.
* Restored full user-guide organization; README is the maintained guide without ReadTheDocs.

Released history is in [CHANGELOG.md](CHANGELOG.md). Development features are not
an announced v0.3.0 release.

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
