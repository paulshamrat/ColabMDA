# HPC and Slurm

Use this platform guide with [apo](../tutorials/apo.md) or [ATP/Mg](../tutorials/atpmg.md).
The scientific commands are identical to Colab. Slurm allocates resources;
active coordinates, trajectories and checkpoints belong in scratch.

## Access and software

Use your cluster's approved access and deployment tools. Run GPU calculations
only in an allocation. Choose an authorized partition, account and compatible
GPU; the example adapters do not select these for you.

In a cluster shell, choose your scratch paths and obtain the public tool:

```bash
export COLABMDA_TOOL=/scratch/$USER/research/tools/ColabMDA
export PROJECT_ROOT=/scratch/$USER/research/projects/colabmda-example
mkdir -p "$(dirname "$COLABMDA_TOOL")" "$PROJECT_ROOT"
git clone --branch main https://github.com/paulshamrat/ColabMDA.git "$COLABMDA_TOOL"
cd "$COLABMDA_TOOL"
```

Replace `/scratch/$USER` if your cluster uses another scratch location. On an
existing checkout, use a reviewed commit or release tag. Record its commit/hash
with the run and keep it fixed while jobs are active.

Use the [installer](../installation.md) with this tool path, or
create the core environment with an existing Conda installation:

```bash
conda env create -f environment.yml
conda activate colabmda
python -m pip install -e .
export COLABMDA_CONDA_ROOT=$(conda info --base)
export COLABMDA_ENV=colabmda
export COLABMDA_PROJECT_DIR="$PROJECT_ROOT"
export COLABMDA_RUN_ID=cdkl5_public_v1
export MPLBACKEND=Agg
```

MODELLER needs its separate environment and your own license. Use the installer's
optional hidden prompt, then follow the scientific guide's environment activation
steps. Keep license values out of commands, Slurm source and logs. Alternatively,
prepare verified models elsewhere and transfer their coordinates, mappings and
hashes with an approved file-transfer tool.

The default environment selects CUDA 13 and needs a supported GPU. T4 and A100
were used in the pilots. CUDA 13 removed offline compilation and library support
for pre-Turing GPUs such as V100
([NVIDIA release notes](https://docs.nvidia.com/cuda/archive/13.0.0/cuda-toolkit-release-notes/)).
A login-node plugin listing does not establish GPU support; the smoke job below
creates an actual CUDA Context and refuses CPU fallback.

## Prepare one system

Follow the selected scientific guide to fetch public inputs, model WT/L119R,
assemble ATP/Mg when needed, stage folders and prepare Amber inputs. Use the
scratch project variables above instead of the Colab paths. Substantial CPU
preparation belongs in a CPU allocation.

Before submitting simulations, check `PREP_DONE` and the recorded identity,
topology, charge and composition gates. Start with one WT system. Install the
thin GPU and CPU adapters in its execution project:

```bash
cd "$PROJECT_ROOT"
mkdir -p scripts/simulation
cp "$COLABMDA_TOOL/examples/hpc/simulate.slurm" scripts/simulation/01_simulate.slurm
cp "$COLABMDA_TOOL/examples/hpc/analyze.slurm" scripts/simulation/02_analyze.slurm
export STATE=apo
export VARIANT=WT
mkdir -p "logs/$STATE/$COLABMDA_RUN_ID/$VARIANT"
```

The adapters call shared `cmda` operations. GPU defaults are four CPUs, 16 GB RAM
and one hour; CPU analysis uses two CPUs, 4 GB RAM and ten minutes. Review these
limits and any site-required Slurm flags. `COLABMDA_CONDA_ROOT` and
`COLABMDA_ENV` select the installed environment.

## Test GPU execution

Set authorized scheduler values for your cluster:

```bash
export PARTITION=your_gpu_partition
export CPU_PARTITION=your_cpu_partition
export ACCOUNT=your_account
export SMOKE_RUN="${COLABMDA_RUN_ID}_smoke"
COLABMDA_RUN_ID="$SMOKE_RUN" cmda stage --state "$STATE" --variant "$VARIANT"
COLABMDA_RUN_ID="$SMOKE_RUN" cmda prep --state "$STATE" --variant "$VARIANT"
smoke_job=$(sbatch --parsable --partition="$PARTITION" --account="$ACCOUNT" \
  --time=00:10:00 --export=ALL \
  --output="logs/$STATE/$COLABMDA_RUN_ID/$VARIANT/smoke_%j.out" \
  --error="logs/$STATE/$COLABMDA_RUN_ID/$VARIANT/smoke_%j.err" \
  scripts/simulation/01_simulate.slurm \
  "$PROJECT_ROOT" "$SMOKE_RUN" "$STATE" "$VARIANT" smoke)
sacct -j "$smoke_job" --format=JobID,State,ExitCode,Elapsed
```

Inspect the smoke log and GPU probe before continuing. Smoke uses a separately
named software test; its outputs do not count as research equilibration.

## Equilibrate and run a short pilot

After smoke acceptance, submit the full 24-stage equilibration:

```bash
eq_job=$(sbatch --parsable --partition="$PARTITION" --account="$ACCOUNT" --export=ALL \
  --output="logs/$STATE/$COLABMDA_RUN_ID/$VARIANT/eq_%j.out" \
  --error="logs/$STATE/$COLABMDA_RUN_ID/$VARIANT/eq_%j.err" \
  scripts/simulation/01_simulate.slurm \
  "$PROJECT_ROOT" "$COLABMDA_RUN_ID" "$STATE" "$VARIANT" eq)
sacct -j "$eq_job" --format=JobID,State,ExitCode,Elapsed
```

Review the physical QC and plots. Production also checks the complete staged QC
gate; a scheduler dependency alone does not establish scientific acceptance.
After QC passes, submit 100 ps for replica r1:

```bash
md_job=$(sbatch --parsable --partition="$PARTITION" --account="$ACCOUNT" --export=ALL \
  --dependency="afterok:$eq_job" \
  --output="logs/$STATE/$COLABMDA_RUN_ID/$VARIANT/md_%j.out" \
  --error="logs/$STATE/$COLABMDA_RUN_ID/$VARIANT/md_%j.err" \
  scripts/simulation/01_simulate.slurm \
  "$PROJECT_ROOT" "$COLABMDA_RUN_ID" "$STATE" "$VARIANT" md r1 0.1)
squeue -j "$md_job"
```

For interrupted jobs, restore the same run ID, inputs and compatible software/GPU,
then follow the scientific guide's checkpoint-recovery commands. Use a fresh run
ID for changed inputs or protocols.

## Analysis and retrieval

Submit the CPU adapter after production:

```bash
sbatch --partition="$CPU_PARTITION" --account="$ACCOUNT" --export=ALL \
  --dependency="afterok:$md_job" \
  --output="logs/$STATE/$COLABMDA_RUN_ID/$VARIANT/analysis_%j.out" \
  --error="logs/$STATE/$COLABMDA_RUN_ID/$VARIANT/analysis_%j.err" \
  scripts/simulation/02_analyze.slurm \
  "$PROJECT_ROOT" "$COLABMDA_RUN_ID" "$STATE" "$VARIANT"
```

This requests no GPU. It merges r1 and writes tables, metadata and figures to
`data/<state>/analysis/<run>/<variant>/r1/`, then checks the required outputs
before `ANALYSIS_DONE`. Retrieve compact reports through approved transfer tools
and verify SHA-256. Heavy trajectories remain in scratch until reviewed archival
transfer. A Slurm `COMPLETED` result alone does not establish convergence.

For a later campaign, derive array bounds from a validated manifest and test its
dependencies separately. These adapters submit no arrays or library campaign.

## Validation scope

Public-input CDKL5 WT apo and ATP/one-Mg passed full 24-stage, 1.96 ns
equilibration, 100 ps production and analysis on A100. Checks covered physical QC,
ten frames at 10–100 ps, native positions 1–303 and ATP43/Mg1 retention.
L119R has CPU modeling/preparation coverage only. The short pilots verify software
operation, not mutation effects or biological convergence. Account-specific run
records belong in the private execution project.
