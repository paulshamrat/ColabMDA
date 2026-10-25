# HPC and Slurm

Use this platform guide with [apo](../apo.md) or [ATP/Mg](../atpmg.md).
The scientific commands are identical to Colab. Slurm allocates GPU resources;
active coordinates, trajectories and checkpoints belong in scratch. Software,
small manifests and compact final reports may live in project storage.

## Access and deployment from Divergent

Our Palmetto access uses the `hpc` CLI exclusively. Authenticate the configured
`hpc` session first. Inspect the scheduler through it; never run protein dynamics
on a login node. Other clusters may require different resources or partitions.

```bash
cd ~/research/tools/ColabMDA
hpc exec -s hpc -- bash -lc 'sinfo -a -h -o "%P %a %l" | grep -E "^compbio |^work1"'
```

For this account, use the owner partition `compbio` and account `compbio` for
short GPU pilots. It provides 52 A100 GPUs. The normal `work1` profile queued our
short probe while the compbio molecular smoke completed in 57 seconds. Check
all accessible partitions, not just the default view. Palmetto assigns QoS;
explicit `--qos` is rejected. The adapter requests one A100, four CPUs, 16 GB RAM
and one hour. Override to ten minutes for `smoke`; review wall time for longer
production. Other accounts must use their authorized partition/account.

Deploy a **clean committed snapshot**, excluding Git metadata and private study
material. In the local terminal, export only the tool runtime files:

```bash
cd ~/research/tools/ColabMDA
COMMIT=$(git rev-parse HEAD)
mkdir -p /tmp/colabmda-hpc-snapshot
git archive "$COMMIT" src examples scripts/setup.sh pyproject.toml environment.yml README.md docs LICENSE \
  | tar -x -C /tmp/colabmda-hpc-snapshot
printf '%s\n' "$COMMIT" > /tmp/colabmda-hpc-snapshot/COMMIT
hpc sync-code -s hpc /tmp/colabmda-hpc-snapshot \
  /scratch/shamrap/research/tools/ColabMDA --execute
```

Use a fresh temporary export directory for another snapshot. Record the commit
and file hashes with the run; never replace source underneath an active run.
The computational project stays separate from the tool. Upload only the reviewed
compact structure/mapping/model inputs needed for the selected workflow, using
`hpc upload`. Do not upload manuscripts or unrelated projects.

## Software and project setup

Run setup from **Divergent** through `hpc exec`. The checked account already has
a suitable `openmm_env`; activate it through the shared conda base and install
only the tool, leaving its scientific libraries unchanged:

```bash
hpc exec -s hpc -- bash -lc '
  set -e
  export COLABMDA_CONDA_ROOT=$(conda info --base)
  source "$COLABMDA_CONDA_ROOT/etc/profile.d/conda.sh"
  export COLABMDA_ENV=openmm_env
  conda activate "$COLABMDA_ENV"
  export COLABMDA_TOOL=/scratch/shamrap/research/tools/ColabMDA
  export PROJECT_ROOT=/scratch/shamrap/research/projects/colabmda-example
  cd "$COLABMDA_TOOL"
  python -m pip install --no-deps -e .
  cmda --help
  mkdir -p "$PROJECT_ROOT"
  cd "$PROJECT_ROOT"
'
```

On a fresh account, create `colabmda` from the tool's `environment.yml` once and
set `COLABMDA_ENV=colabmda` consistently in setup and submission commands. A shared
read-only conda base can activate a writable user environment. Do not modify the
system base or use sudo. Each `hpc exec` starts a new shell; include activation,
paths and `cd` in each operation rather than assuming variables persist.

Do not recreate an existing environment to resume a run. Check the pinned CUDA
runtime against the allocated GPU driver in a small allocation before dynamics.
The default environment currently requests CUDA 13; cluster compatibility is a
required pilot check, not assumed from Colab success. The adapter requests A100
rather than Volta/V100: CUDA 13 removed offline compilation and library support
for pre-Turing architectures ([NVIDIA release notes](https://docs.nvidia.com/cuda/archive/13.0.0/cuda-toolkit-release-notes/)).

Modeling needs the separate licensed `modeller_env`. If models are prepared on
Divergent, upload their verified coordinates, mappings and hashes instead. A
MODELLER key remains a user secret; retain the installer's optional hidden prompt
and never put the key in Slurm source or logs. See the root README for the
installer and environment separation.

Follow the scientific guide for input/model checks and `cmda stage`/`cmda prep`,
using this scratch `PROJECT_ROOT`. MODELLER/Amber preparation must complete and
produce `PREP_DONE` before GPU submission. A substantial CPU preparation job
belongs in a CPU allocation rather than a login node.

The inspected Palmetto environment `/home/shamrap/.conda/envs/openmm_env` has
OpenMM 8.5.2, Python 3.11.15, AmberTools and the analysis dependencies. Selecting
it uses `COLABMDA_ENV=openmm_env`; a fresh installation defaults to `colabmda`.
A CPU-only login-node plugin listing does not establish GPU support. Test an
actual CUDA Context in an allocated GPU job before selecting the environment.
The worker enforces `COLABMDA_REQUIRE_GPU=1` and refuses CPU fallback.

## A single WT pilot

Install the thin scheduler adapter in the execution project. It calls the same
`cmda sim` commands documented in the scientific guides; it is not another MD
implementation. There is no corresponding Colab runner because Colab uses those
commands directly.

```bash
cd "$PROJECT_ROOT"
mkdir -p scripts/simulation
cp -n "$COLABMDA_TOOL/examples/hpc/simulate.slurm" scripts/simulation/01_simulate.slurm
mkdir -p logs/atpmg/$COLABMDA_RUN_ID/WT
```

Submit through `hpc` from Divergent. The following commands set their complete
remote environment. These examples select the inspected `openmm_env`; select
`colabmda` instead if that is the environment you installed.
First run equilibration, inspect its JSON/PNG and job result, then submit a
100 ps production validation. Do not start the full panel before this succeeds.

```bash
hpc exec -s hpc -- bash -lc '
  export COLABMDA_CONDA_ROOT=$(conda info --base)
  export COLABMDA_ENV=openmm_env
  export PYTHONUNBUFFERED=1
  project=/scratch/shamrap/research/projects/colabmda-example
  run=cdkl5_public_v1
  cd "$project"
  sbatch --partition=compbio --account=compbio --output="logs/atpmg/$run/WT/eq_%j.out" --error="logs/atpmg/$run/WT/eq_%j.err" \
    scripts/simulation/01_simulate.slurm "$project" "$run" atpmg WT eq
'
```

After successful equilibration and physical review:

```bash
hpc exec -s hpc -- bash -lc '
  export COLABMDA_CONDA_ROOT=$(conda info --base)
  export COLABMDA_ENV=openmm_env
  export PYTHONUNBUFFERED=1
  project=/scratch/shamrap/research/projects/colabmda-example
  run=cdkl5_public_v1
  cd "$project"
  sbatch --partition=compbio --account=compbio --output="logs/atpmg/$run/WT/md_%j.out" --error="logs/atpmg/$run/WT/md_%j.err" \
    scripts/simulation/01_simulate.slurm "$project" "$run" atpmg WT md r1 0.1
'
```

Use `apo` for protein in water; create matching log folders.
Production reruns strict equilibration QC. Wall-time termination does not imply
completion: inspect job exit status and validated output records, then resubmit
the same command to recover a compatible run from complete stages/chunks.
Do not promise automatic indefinite continuation from this single-job adapter.

### Small molecular software check

Before full equilibration, the adapter supports a separate `smoke` stage. Give
it its own run ID and verified preparation. It performs full initial minimization,
then 50-step quick NVT/NPT and a five-frame production check, using explicit
`--software-test`. `SMOKE_DONE` requires finite coordinates and a nonempty
checkpoint. This checks the GPU/software route; it is not the 24-stage protocol.
Use a fresh run ID for research `eq` and `md`, retaining the same reviewed inputs.
Never relabel the smoke output as research equilibration.

## Monitoring, analysis and retrieval

```bash
hpc exec -s hpc -- bash -lc 'squeue -h -u "$USER" -o "%i %j %T %M %R" | head -25'
hpc exec -s hpc -- sacct -j JOB_ID --format=JobID,State,ExitCode,Elapsed
```

Replace `JOB_ID` with the actual submitted ID. Install the CPU analysis adapter
alongside the simulation adapter, then submit it after production completes:

```bash
hpc exec -s hpc -- bash -lc '
  set -e
  project=/scratch/shamrap/research/projects/colabmda-example
  tool=/scratch/shamrap/research/tools/ColabMDA
  export COLABMDA_CONDA_ROOT=$(conda info --base)
  export COLABMDA_ENV=openmm_env
  cd "$project"
  cp "$tool/examples/hpc/analyze.slurm" scripts/simulation/02_analyze.slurm
  sbatch --output=logs/atpmg/cdkl5_public_v1/WT/analysis_%j.out \
    --error=logs/atpmg/cdkl5_public_v1/WT/analysis_%j.err \
    scripts/simulation/02_analyze.slurm "$project" cdkl5_public_v1 atpmg WT
'
```

This uses two CPUs, 4 GB RAM and ten minutes on compbio, with no GPU request.
It calls the shared merge/analysis commands and writes r1 reports under
`data/<state>/analysis/<run>/<variant>/r1/`. It verifies nonempty tables/figures
before writing `ANALYSIS_DONE`. Retrieve compact reports with `hpc download`; verify SHA-256 against a
remote inventory. Heavy trajectories remain on scratch until a reviewed archival
transfer. A Slurm `COMPLETED` result alone does not prove convergence.

For a later panel campaign, derive array bounds dynamically from a validated
manifest and connect equilibration/production with Slurm dependencies. That
campaign wrapper should be tested separately after the single-system pilot;
this adapter deliberately submits no arrays or full-library campaign.

## Validation status

The earlier reviewed-input tests below are distinct from the new public-input
workflow. Both public WT pilots completed on compbio with run ID
`cdkl5_public_v1`, frozen tool snapshot `0d7f4b3`:

| State | Full 24-stage equilibration | 100 ps production | CPU analysis |
| --- | --- | --- | --- |
| apo | 18148627, 17:41 | 18148643, 1:17 | 18148795, 14 s |
| ATP/Mg | 18148635, 17:45 | 18148785, 1:19 | 18148810, 4 s |

All six jobs completed with exit 0:0. Strict timing/thermodynamic QC, ten finite
frames per state, residue labels 1–303 and ATP43/Mg1 retention passed. Compact
collection verified 59 files per run; all eight QC/analysis PNGs were reviewed.
L119R passed public modeling/preparation checks only. These are software pilots,
not evidence of biological convergence.

### Earlier reviewed-input pilots

Palmetto bridge, scheduler, scratch layout and existing scientific libraries have
been inspected. Committed tool/project snapshots and reviewed inputs are deployed
and hash-verified. The tool is installed without changing the scientific libraries.
The two normal-profile queued checks were cancelled after verifying they had not
started. Compbio smoke job 18088560 completed successfully on an A100 in 57 seconds,
with five finite frames. Its tiny software equilibration does not pass research
QC and is not reused. Full ATP/Mg equilibration job 18088721 completed in 17:04,
passing strict 24-stage timing and thermodynamic QC. Production job 18090842
completed 100 ps in 74 seconds; CPU analysis job 18091098 completed in ten seconds.
The ten finite frames retain ATP (43 atoms) and one Mg, with all 303 C-alpha atoms
and residue labels 1–303. Colab independently passed full apo QC and 200 ps
production, including an actual interruption/resume test. These are WT software
pilots, not biological convergence or GPU validation of L119R.

Palmetto currently rejects user changes to queued jobs. If a request must change,
first verify that it is your pending job, cancel it, and resubmit; never alter
another campaign. Record replacement IDs and dependencies. Source snapshots and
scientific inputs remain fixed while active simulations execute.
