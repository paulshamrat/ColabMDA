# ColabMDA

[![Release](https://img.shields.io/github/v/tag/paulshamrat/ColabMDA)](https://github.com/paulshamrat/ColabMDA/tags) [![CI](https://github.com/paulshamrat/ColabMDA/actions/workflows/python-test.yml/badge.svg?branch=main)](https://github.com/paulshamrat/ColabMDA/actions/workflows/python-test.yml) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![Repository](https://img.shields.io/badge/GitHub-Repository-blue?logo=github)](https://github.com/paulshamrat/ColabMDA) [![Issues](https://img.shields.io/github/issues/paulshamrat/ColabMDA)](https://github.com/paulshamrat/ColabMDA/issues) [![Black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black) [![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

**ColabMDA** helps you run molecular dynamics simulations on Google Colab and recover
from session interruptions. It saves progress to Google Drive so you can resume
from the last completed checkpoint. It supports protein modeling, simulation and
analysis through terminal commands, with the same workflow available on HPC
clusters and Linux workstations.

Start with [Quick Start](#1-platform-setup), then follow the
[Simulation Workflow](#3-simulation-protocol-cmda). For a short pilot, use the
[Short Pilot](#4-short-4bgq-pilot).

## 1. Quick Start <a id="1-platform-setup"></a>

**Runtime:** Google Colab GPU (T4 / A100 / L4) | **Persistent Storage:** Google Drive  
**Scientific Engine:** OpenMM · AmberTools · MODELLER · MDAnalysis · MDTraj · PROPKA

---

### Step 1: Connect to Google Colab GPU

Run these commands in your terminal to authenticate, create your cloud GPU instance, mount Drive, and connect:

```bash
# 1. Install Colab CLI & authenticate (one-time browser login)
pip install google-colab-cli
colab sessions

# 2. Launch Google Colab GPU session & mount Google Drive
colab new -s colabmda-example --gpu T4
# Click the OAuth URL in your browser to grant Drive access, then press Enter in terminal:
colab drivemount --session colabmda-example

# 3. Connect to the interactive Colab console
colab console -s colabmda-example
```

---

### Step 2: Install ColabMDA & Activate Environment

Now **inside the Colab console**, clone the repository and run the automated installer:

```bash
# 1. Clone ColabMDA on the GPU cloud instance
git clone --branch dev1 https://github.com/paulshamrat/ColabMDA.git /content/ColabMDA
cd /content/ColabMDA

# 2. Run the automated installer (sets up OpenMM, AmberTools, MODELLER, etc.)
bash scripts/setup.sh

# 3. Activate the scientific environment
source /root/miniforge3/etc/profile.d/conda.sh
conda activate colabmda
export MPLBACKEND=Agg

# 4. Verify installation
cmda --help
```

> 🔑 **Modeller License:** During setup, the installer will pause and prompt for an optional free MODELLER key for protein modeling. If you don't have one, register at [salilab.org/modeller/registration.html](https://salilab.org/modeller/registration.html).

<details>
<summary><b>🛠️ Developer Tip: Installing an existing checkout or custom branch</b></summary>

If you already have a custom checkout in `/content/ColabMDA` and want to install it directly in editable mode:
```bash
pip install -e /content/ColabMDA
```
</details>

---

### 💡 Session Recovery & Checkpointing

If your Colab session times out, your checkpoints and completed chunks are safely preserved on your Google Drive under `/content/drive/MyDrive/research/projects/`:
1. Re-launch the session and reconnect:
   ```bash
   colab new -s colabmda-example --gpu T4
   colab drivemount --session colabmda-example
   colab console -s colabmda-example
   ```
2. Inside the console, activate the environment (`conda activate colabmda`).
3. Re-run your simulation command (`cmda sim ...`). ColabMDA automatically detects the saved checkpoints on Drive and resumes from the last completed chunk without lost time.

> 💡 **HPC Clusters & Local Workstations:** The simulation commands in Sections 3 and 4 are identical across all platforms. Cluster (SLURM) and workstation deployment guides are documented in [`docs/platforms/hpc.md`](docs/platforms/hpc.md).


---

## 2. Project Architecture

Keep the tool checkout separate from each simulation project. The tool contains
code; the project contains inputs, results and logs. No project YAML is required.

### 2.1. Tool repository <a id="project-structure"></a>

```text
ColabMDA/
├── src/colabmda/
│   ├── cli.py             # cmda commands
│   ├── core/              # Projects, mapping and workflow checks
│   ├── modeller/          # Protein modeling and mutation
│   ├── amber/             # Parameterization and trajectory processing
│   ├── openmm/            # Simulation, QC, analysis and figures
│   └── remote/            # Colab setup and control
├── examples/              # Small workflows and HPC adapters
├── scripts/setup.sh       # Installer
├── environment.yml        # Core dependencies
├── docs/                  # State and platform guides
└── tests/                 # CPU tests
```

### 2.2. Simulation project

Use this shared layout on Colab, HPC or a local workstation:

```text
colabmda-example/                                 <- Project workspace (Google Drive or local)
├── data/{apo,atpmg}/                             <- State tracks (apo: protein in water; atpmg: ATP + Mg²⁺)
│   ├── structure/
│   │   ├── raw/                                  <- Downloaded 4BGQ.cif/pdb and uniprot.fasta
│   │   ├── mapping/                              <- 3-way coordinate/sequence map (residues.json)
│   │   └── models/{WT,L119R}/                    <- Validated models (WT.pdb, L119R.pdb, MODEL_DONE)
│   ├── simulation/
│   │   └── {WT,L119R}/
│   │       ├── 01_prep/                          <- Solvated system & Amber ff19SB/OPC prmtop/inpcrd
│   │       ├── 02_eq/                            <- Shared 24-stage equilibration hub & qc.png
│   │       └── 03_md/{r1,r2,r3}/                 <- Production trajectory chunks & checkpoints
│   └── analysis/
│       ├── {WT,L119R}/{r1,r2,r3}/                <- Single-system tables (rmsd.csv, rg.csv, rmsf.csv)
│       └── compare/figures/                      <- WT vs Mutant overlay comparison figures
└── logs/{apo,atpmg}/{WT,L119R}/
```

> 💡 **Custom Projects:** ColabMDA is directory-agnostic and executes inside whichever folder you enter (or pass via `--root <path>`). You can replace `colabmda-example` with any folder name for your own target.

---

## 3. Simulation Protocol (cmda)

### System Overview (CDKL5 Kinase Case Study)

| Parameter | Specification | Scientific Context / Motivation |
| :--- | :--- | :--- |
| **Protein Target** | Human CDKL5 Kinase Domain | Model serine/threonine protein kinase |
| **Crystal Structure** | PDB [4BGQ](https://www.rcsb.org/structure/4BGQ), Chain A | Raw experimental X-ray crystal coordinates |
| **UniProt Sequence** | [O76039](https://www.uniprot.org/uniprotkb/O76039/entry) (Residues 1–303) | Canonical reference sequence with biological numbering |
| **Experimental Gap** | 276 resolved, **27 loop residues missing** | Unresolved loops require homology completion up to residue 303 |
| **Wild Type (WT)** | Full 1–303 modeled kinase | Fills the 27 missing loop residues preserving biological numbering |
| **Pathogenic Variant** | **L119R** (Leucine 119 $\to$ Arginine) | Point mutation in ATP-binding hinge region linked to CDD |
| **Two Tracks** | **Apo** (water/salt) · **Holo** (ATP + 1 $\text{Mg}^{2+}$) | Comparative conformational dynamics |

Start with WT, then repeat with L119R. Both tracks start from public archive and
UniProt inputs and run through modeling, preparation, simulation and analysis.
See [Validation](#6-validation) for tested systems and limits.

### 3.1. Structure Modeling & Verification (WT & L119R)

Navigate to your persistent project workspace on Google Drive:

```bash
mkdir -p /content/drive/MyDrive/research/projects/colabmda-example
cd /content/drive/MyDrive/research/projects/colabmda-example
```

Fetch the public structure (4BGQ) and UniProt sequence (O76039), mapping native residue numbering (residues 1–303):

```bash
# 1. Fetch template coordinates and UniProt mapping
cmda fetch --state apo --pdb-id 4BGQ --uniprot-id O76039 --range 1 303 --mutations L119R

# 2. Build full Wild-Type model (fills 27 missing loop residues 1–303)
cmda model --state apo --variant WT

# 3. Build L119R point mutant on the verified WT backbone
cmda model --state apo --variant L119R
```

The `cmda model` command automatically calls `cmda check` upon completion, checking:
- Continuous biological numbering (1–303) matching UniProt O76039.
- Complete coordinate density (all 303 C$\alpha$ atoms present and finite).
- Verified stereochemistry and touches `MODEL_DONE` receipt.

> 💡 **Bringing Your Own Pre-Built Models (AlphaFold / SwissModel):**
> If you already have a pre-built PDB model and do not want to run MODELLER, you can place your model directly and gate it with `cmda check`:
> ```bash
> mkdir -p data/apo/structure/models/WT
> cp /path/to/my_model.pdb data/apo/structure/models/WT/WT.pdb
> cmda check --state apo --variant WT
> ```
> `cmda check` verifies sequence alignment and coordinate density, creating `MODEL_DONE` so downstream preparation can proceed immediately without MODELLER.

### 3.2. Track 1: Apo Dynamics (WT & L119R)

#### 3.2.1. Wild-Type (WT) Pipeline

Simulate the verified Wild-Type protein in physiological water and salt:

```bash
# 1. Solvate & Parameterize (Amber ff19SB + OPC water + 150 mM NaCl)
cmda prep --state apo --variant WT
# -> Outputs: data/apo/simulation/WT/01_prep/ (system.prmtop, system.inpcrd)

# 2. 24-Stage Equilibration (1.96 ns heating & density ladder) + QC check
cmda sim --state apo --variant WT --stage eq
# -> Outputs: data/apo/simulation/WT/02_eq/ (equilibration_qc.json, qc.png)

# 3. Production MD (e.g. 100 ns, saving frames every 10 ps)
cmda sim --state apo --variant WT --stage md --replica r1 --total-ns 100 --checkpoint-ps 1000 --traj-interval 10
# -> Outputs: data/apo/simulation/WT/03_md/r1/ (prod_*.dcd chunks, prod.chk)
```

#### 3.2.2. Mutant (L119R) Pipeline

Once WT is validated, repeat the identical physical pipeline for the L119R variant:

```bash
# 1. Solvate & Parameterize L119R
cmda prep --state apo --variant L119R
# -> Outputs: data/apo/simulation/L119R/01_prep/

# 2. 24-Stage Equilibration & QC
cmda sim --state apo --variant L119R --stage eq
# -> Outputs: data/apo/simulation/L119R/02_eq/

# 3. Production MD
cmda sim --state apo --variant L119R --stage md --replica r1 --total-ns 100 --checkpoint-ps 1000 --traj-interval 10
# -> Outputs: data/apo/simulation/L119R/03_md/r1/
```

### 3.3. Track 2: Holo Dynamics (ATP + Mg²⁺)

#### 3.3.1. Scope & Chemistry

Use the public WT/L119R protein models from section 3.1. The `atpmg` state contains
one ATP and one Mg, without phosphorylation or a substrate peptide.

#### 3.3.2. ATP & Mg²⁺ Placement

Download the public 8FP5 ATP/Mg donor and assemble WT:

```bash
cmda assemble --state atpmg --template 8FP5 --relax-sidechains
```

This command downloads the donor, copies the checked protein model and mapping,
fits homologous core residues, and writes `data/atpmg/structure/models/WT/complex.pdb`.
No study-owner files or manual file copying are needed.

[8FP5](https://www.rcsb.org/structure/8FP5) is CDK2 with ATP/Mg. Its transfer to
CDKL5 creates a modeled pose, not an experimental CDKL5 complex or a reproduction
of the earlier study. The fit uses aligned, experimentally resolved core C-alpha
positions, rejects outliers and records coverage and RMSD.

`--relax-sidechains` locally minimizes pocket sidechains with the protein backbone,
ATP and Mg fixed. This geometry step uses Amber parameters in vacuum; it is not
solvent equilibration or a binding-affinity calculation. Assembly fails if the fit,
fixed-coordinate, energy or final clash checks fail. Review the pocket and Mg
coordination in `complex.pdb` and the recorded checks in `complex.json`.

For L119R, after WT validation:

```bash
cmda assemble --state atpmg --template 8FP5 --variant L119R --relax-sidechains
```

The source apo backbone stays unchanged. Each holo backbone matches its apo
model; only selected pocket sidechains are relaxed. Original archives are retained.

#### 3.3.3. Holo Parameterization (cmda prep)

```bash
cmda prep --state atpmg --variant WT
cmda prep --state atpmg --variant L119R
```

Preparation selects that variant's assembled complex. It uses the Hu et al. B3
ATP parameters, ff19SB/OPC and complete water-matched Li–Merz 12–6 ions.
The actual topology must contain ATP with 43 atoms and charge −4, plus one Mg
with charge +2. `topology_check.json` also checks total charge and ion LJ values.

### 3.4. Equilibration, QC & Production

#### 3.4.1. Staged Equilibration (24-Stage Ladder)

From your project directory, launch equilibration for the prepared variant:

```bash
cmda sim --state apo --variant WT --stage eq
```

For the ATP/Mg system, use `--state atpmg`. Repeat with `--variant L119R` after WT validation.
Do not set `COLABMDA_TEST_1_24` for research. Quick/reduced runs require
`--software-test`. The 4BGQ pilot in section 4 keeps the full research protocol.

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

Numbered state/JSON records support recovery from completed stages. Optional quick
software tests use one short NVT and NPT interval; the section 4 pilot uses this
full protocol. The legacy combined `openmm run` adapter retains ff14SB/TIP3P defaults
and delegates to quick equilibration; it should not be confused with the staged route.

#### 3.4.2. Physical Settings & QC

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

#### 3.4.3. Production MD (r1, r2, r3)

After full equilibration passes QC, run a short WT production check:

```bash
cmda sim --state apo --variant WT --stage md --replica r1 \
  --total-ns 0.1 --checkpoint-ps 100 --traj-interval 10
```

This produces ten frames over 100 ps. Review analysis and cofactor geometry before
running L119R. For a scientific campaign, choose a justified duration and replica
count. The following illustrates three 100 ns targets:

```bash
cmda sim --state apo --variant WT --stage md --replica r1,r2,r3 \
  --total-ns 100 --checkpoint-ps 1000 --traj-interval 10
```

Repeat for each reviewed state/variant and its run ID. Targets are per replica.
The high-level route processes the listed replicas sequentially; parallel batch
execution needs an appropriate scheduler design. Replicas can share an equilibration
endpoint; this does not establish converged or uncorrelated sampling.

#### 3.4.4. Timestep, Steps & Saved Output

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
= **1,960 ps (1.96 ns)** before production. The optional quick example uses 100 ps NVT +
100 ps NPT = **200 ps (0.2 ns)** before its 1 ns production. Equilibration is not
part of the production target. Three 100 ns replicas provide 300 ns aggregate
production, not one continuous 300 ns trajectory.

Production currently fixes the timestep at 2 fs; there is no production timestep
CLI option. Changing `--traj-interval` changes saved output, not the integrator.
The high-level `--steps` testing override changes the requested production step
count as well as quick equilibration. Use `--total-ns` for research duration. Use intervals corresponding to whole integration steps;
the implementation converts requested intervals to integer step counts.

#### 3.4.5. Interruption & Checkpoint Recovery

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

### 3.5. Autoimaging & Solvent Stripping (cmda merge)

#### 3.5.1. PBC Autoimaging

```bash
cmda merge --state apo --variant WT --replicas r1
cmda merge --state apo --variant L119R --replicas r1
```

Use `--replicas r1,r2,r3` after all three replicas have finished. Repeat for each
simulated variant. `cpptraj` autoimages the completed chunks and
writes both a solute stream and a C-alpha stream. Storage savings depend on system
composition; inspect the results rather than assuming a fixed reduction.

#### 3.5.2. Reference Topologies

| Output | Matching topology | Use |
| :--- | :--- | :--- |
| `prod_*to*ps.dcd` | Original solvated topology | Raw chunk outputs, including solvent. |
| `prod_keep_ligand.dcd` | `stripped.prmtop` | Solute and retained cofactors/ions. |
| `prod_ca.dcd` | `ca.prmtop` | C-alpha-only trajectory. |
| `prod_full.log` | Not a topology | Merged production times and thermodynamic measurements. |

Never pair a stripped trajectory with the unstripped solvated PDB. Confirm atom
counts and ordering. A viewer requiring PDB needs an exported, matching topology;
verify its residue-label mapping rather than treating Amber indices as author numbers.

### 3.6. Trajectory Analysis & Comparison

#### 3.6.1. System Metrics (cmda analyze)

```bash
cmda analyze --state apo --variant WT --replicas r1
cmda analyze --state apo --variant L119R --replicas r1
```

The default is r1. Use `--replicas r1,r2,r3` to analyze completed replicas and
write their mean and population standard deviation under `aggregate/`. Time and
native residue axes must match; mismatches stop aggregation. RMSD uses backbone fitting to the first
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

For completed state-based runs that still contain `simulation/<run>/<variant>/04_analysis`,
set the original run ID and migrate their r1 reports explicitly:

```bash
cd "$PROJECT_ROOT"
cmda stage --state "$STATE" --variant "$VARIANT" --migrate-analysis
```

Migration checks every file hash and refuses incomplete reports, symlinks, heavy
simulation data or conflicting destinations. It verifies the new copy before
removing the old folder and records the path mapping under
`data/00_setup/analysis_migrations/`. Original provenance is retained unchanged.

#### 3.6.2. Comparative Overlays (WT vs Mutant)

Generate direct publication-grade overlay comparison figures across variants:

```bash
cmda compare --state apo --variants WT,L119R
```

This automatically aggregates all completed replicas for `WT` and `L119R`, generates 300-DPI publication plots with mean and shaded standard deviations, and writes them into `data/apo/analysis/compare/figures/` (`rmsd_avg.png`, `rg_avg.png`, `rmsf_avg.png`).

*(For custom directories or external paths, the lower-level `cmda openmm compare --series ...` remains available).*

#### 3.6.3. Frame Counts, Sampling & Stride

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

#### 3.6.4. Storage Calculation

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

#### 3.6.5. Frame & Coordinate Verification

Run this in the MD environment after merging a selected research system. It uses
the matching stripped topology, reports the saved log times, and checks frame
counts and finite coordinates without assuming that a filename proves completion:

<details>
<summary>Optional Python check of frames and coordinates</summary>

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

</details>

Choose the state/variant/replica you actually ran. Check expected times, sampling
and native residue labels against the run settings. CDKL5 has 303 C-alpha atoms.

For CDKL5 compact downloads, follow the [apo](docs/apo.md#7-retrieve-compact-results)
or [ATP/Mg](docs/atpmg.md#7-retrieve-compact-results) collection commands. They use
the same state/run/variant paths as this workflow and verify SHA-256 after download.
Heavy trajectories stay on Drive or HPC scratch.

### 3.7. CDKL5 Visualization

#### 3.7.1. Workstation or CPU rendering environment

Use a separate environment and install this same tool checkout there:

```bash
conda create -y -n pymol_env -c conda-forge python=3.11 pymol-open-source pillow numpy mdanalysis biopython
conda activate pymol_env
python -m pip install -e .
```

Run the view commands where the merged run data is available. Rendering uses CPU,
not a GPU. For HPC data, render on a suitable CPU node or transfer only selected
structure samples; full trajectories can stay in scratch.

#### 3.7.2. Structure and trajectory views

Use a numbered starting model or a saved trajectory frame:

```bash
cd "$PROJECT_ROOT"
cmda view --state apo --variant WT --structure --png
cmda view --state atpmg --variant WT --frame last --png
```

Omit `--png` to open the fitted frame in PyMOL; add `--frame all` for a movie.
Residue 119 is orange, ATP is
purple and Mg is blue. `--resi` selects another native residue. Missing data,
sequence mismatches or invalid frame numbers stop the command.

The default camera uses the activation project's full-kinase structural overview.
All frames fit to a common numbered WT reference using core C-alpha positions
30–220; mapped missing positions are excluded. The camera keeps that orientation,
reframes once to WT, then reuses the same framing across panels. Use
`--camera network` for the activation network view, including its 90° rotation.
No activation-project files are required at runtime.

#### 3.7.3. Comparable snapshot grids

Compare first/last WT frames across states:

```bash
cmda snapshots --state apo,atpmg --variant WT --frames 1,last
```

Compare the starting models without implying that L119R MD has been run:

```bash
cmda snapshots --state apo,atpmg --variant WT,L119R --structure
```

For trajectory comparisons of variants, omit `--structure` only after their
matched runs and merges have completed. Frame numbers are one-based; time labels
come from the production log, so frame 1 is not relabeled as time zero.

Outputs are 7-inch, 300-DPI PNGs under the first selected system's
`data/<state>/analysis/<run>/<variant>/<replica>/figures/`, with camera, input hashes, residue labels and frame/time
metadata beside them. Use `--output` for another filename and `--reference` for an
explicit common WT reference. Custom input JSON remains an advanced option;
the default CDKL5 workflow needs none. Custom time labels must be replaced with
a matching production log. Starting-model grids use `model_snapshots.png`;
trajectory grids use `structure_snapshots.png`. Starting-model views use the
variant's `figures/` folder, outside replica folders.

---

## 4. Short 4BGQ Pilot

### 4.1. Run WT first

Complete setup and public model/preparation steps in sections 3.1–3.3. In the MD environment on a GPU runtime or allocated GPU node:

```bash
cmda sim --state apo --variant WT --stage eq
cmda sim --state apo --variant WT --stage md --total-ns 0.1 --checkpoint-ps 100
cmda merge --state apo --variant WT --replicas r1
cmda analyze --state apo --variant WT --replicas r1
```

Set `STATE=atpmg` and repeat for its prepared WT complex. Keep full 24-stage
1.96 ns equilibration; the 100 ps production is a software pilot. Our A100 tests
took about 18 minutes for equilibration and 80 seconds for production per state.
Runtime depends on hardware and system size.

### 4.2. Check results and recover

Require passed `02_eq/equilibration_qc.json`, ten frames at 10–100 ps, finite
analysis tables and residue labels 1–303. Inspect the protein, ATP/Mg geometry
and the figures. Outputs use the shared state/run/variant layout in section 2.

Rerun the same compatible commands to resume complete stages/chunks, as described
in section 3.4.5. A 0.1 ns target can be extended by requesting a larger absolute
`--total-ns`; recalculate merge and analysis after extension. Run L119R only after
WT validation. These short runs do not establish mutation effects or convergence.

For compact collection and shutdown, use the [apo](docs/apo.md#7-retrieve-compact-results)
or [ATP/Mg](docs/atpmg.md#7-retrieve-compact-results) guide. Keep heavy trajectories
on Drive or HPC scratch. Optional older tutorials are in [examples](docs/examples.md).

---

## 5. Command Reference

| Command | Purpose |
| :--- | :--- |
| `cmda fetch` | Download structures and map UniProt positions. |
| `cmda model` | Build WT or a variant from checked inputs. |
| `cmda check` | Verify sequence, numbering and coordinates. |
| `cmda assemble` | Assemble ATP/Mg from a donor or reviewed pose. |
| `cmda stage` | Create folders or install an example. |
| `cmda prep` | Build Amber topology and coordinates. |
| `cmda sim` | Equilibrate, check QC and run production. |
| `cmda merge` | Autoimage and strip trajectories. |
| `cmda analyze` | Write RMSD, Rg and RMSF tables and plots. |
| `cmda compare` | Generate multi-variant overlay comparison plots (RMSD, Rg, RMSF). |
| `cmda openmm compare` | Plot custom series of analysis directories. |
| `cmda modeller` | Use custom build or mutation inputs. |
| `cmda view` | Open a matching topology/trajectory in PyMOL. |
| `cmda snapshots` | Render a configured snapshot grid. |
| `cmda colab` | Set up and control Colab sessions. |

Use `cmda --help` to list commands and `cmda COMMAND --help` for options, for
example `cmda sim --help`. Workflow sections show only the commands needed for
each step. Top-level `cmda compare --state apo --variants WT,L119R` automatically
renders aggregate multi-variant overlay figures.

---

## 6. Validation

The public-input example passed WT/L119R modeling and preparation checks.
WT apo and ATP/Mg each passed full 24-stage, 1.96 ns equilibration, 100 ps
production and analysis on Palmetto A100. L119R GPU runs remain untested.
These short pilots verify the workflow, not mutation effects or convergence.

The earlier study-based pilots used different starting inputs.

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

- Keep analysis beside structure and simulation, with explicit replica selection,
  checked aggregation and hash-verified migration of completed r1 reports.

- Use terminal and Colab CLI setup throughout the guides; remove the web notebook route.
- Make the CDKL5 WT/L119R apo and ATP/one-Mg walkthrough independent of private
  study files, with short public-input commands and a full 24-stage protocol.
- Record homologous ATP/Mg transfer, fixed-backbone pocket relaxation, native
  numbering, actual topology checks and continuous recovery clocks.
- Apply the requested MODELLER seed to WT construction and mutation.
- Use activation-project camera orientations with a shared WT reference for
  CDKL5 views; retain real frame times and remove protein-specific legacy styling.
- Keep architecture and command lookup consolidated. Optional templates remain
  in [examples](docs/examples.md).

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
