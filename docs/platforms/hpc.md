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
hpc exec -s hpc -- bash -lc 'sinfo -h -o "%P %a %l %G" | head -12'
```

On the checked Palmetto configuration, `work1` offers A100 GPUs and allows up to
three days. The shipped pilot adapter requests one A100, four CPUs, 16 GB RAM and
four hours. Resource discovery is verified; a successful scheduler query is not
a completed GPU pilot.

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

Inside a terminal opened through the configured `hpc` CLI, set these paths for
your account. Use an existing **writable** conda installation or install Miniforge
under your account. Do not modify the system conda root or use sudo.

```bash
export COLABMDA_CONDA_ROOT=/path/to/your/miniforge3
source "$COLABMDA_CONDA_ROOT/etc/profile.d/conda.sh"
export COLABMDA_TOOL=/scratch/shamrap/research/tools/ColabMDA
export PROJECT_ROOT=/scratch/shamrap/research/projects/colabmda-example
export COLABMDA_RUN_ID=cdkl5_reviewed_v1
cd "$COLABMDA_TOOL"
conda env create -f environment.yml
conda activate colabmda
python -m pip install -e .
cmda --help
mkdir -p "$PROJECT_ROOT"
cd "$PROJECT_ROOT"
```

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

## A single WT pilot

Install the thin scheduler adapter in the execution project. It calls the same
`cmda sim` commands documented in the scientific guides; it is not another MD
implementation. There is no corresponding Colab runner because Colab uses those
commands directly.

```bash
cd "$PROJECT_ROOT"
mkdir -p scripts/simulation
cp -n "$COLABMDA_TOOL/examples/hpc/simulate.slurm" scripts/simulation/01_simulate.slurm
mkdir -p logs/apo/$COLABMDA_RUN_ID/WT
```

Submit through `hpc` from Divergent. The following commands set their complete
remote environment; change the conda path to the verified user installation.
First run equilibration, inspect its JSON/PNG and job result, then submit a
100 ps production validation. Do not start the full panel before this succeeds.

```bash
hpc exec -s hpc -- bash -lc '
  export COLABMDA_CONDA_ROOT=/path/to/your/miniforge3
  project=/scratch/shamrap/research/projects/colabmda-example
  run=cdkl5_reviewed_v1
  cd "$project"
  sbatch --output="logs/apo/$run/WT/eq_%j.out" --error="logs/apo/$run/WT/eq_%j.err" \
    scripts/simulation/01_simulate.slurm "$project" "$run" apo WT eq
'
```

After successful equilibration and physical review:

```bash
hpc exec -s hpc -- bash -lc '
  export COLABMDA_CONDA_ROOT=/path/to/your/miniforge3
  project=/scratch/shamrap/research/projects/colabmda-example
  run=cdkl5_reviewed_v1
  cd "$project"
  sbatch --output="logs/apo/$run/WT/md_%j.out" --error="logs/apo/$run/WT/md_%j.err" \
    scripts/simulation/01_simulate.slurm "$project" "$run" apo WT md r1 0.1
'
```

Use `atpmg` for its independently prepared state; create matching log folders.
Production reruns strict equilibration QC. Wall-time termination does not imply
completion: inspect job exit status and validated output records, then resubmit
the same command to recover a compatible run from complete stages/chunks.
Do not promise automatic indefinite continuation from this single-job adapter.

## Monitoring, analysis and retrieval

```bash
hpc exec -s hpc -- bash -lc 'squeue -h -u "$USER" -o "%i %j %T %M %R" | head -25'
hpc exec -s hpc -- sacct -j JOB_ID --format=JobID,State,ExitCode,Elapsed
```

Replace `JOB_ID` with the actual submitted ID. Run merge/analysis from the
scientific guide in a suitable allocation, using the same project/run/state
variables. Retrieve compact reports with `hpc download`; verify SHA-256 against a
remote inventory. Heavy trajectories remain on scratch until a reviewed archival
transfer. A Slurm `COMPLETED` result alone does not prove convergence.

For a later panel campaign, derive array bounds dynamically from a validated
manifest and connect equilibration/production with Slurm dependencies. That
campaign wrapper should be tested separately after the single-system pilot;
this adapter deliberately submits no arrays or full-library campaign.

## Validation status

Palmetto bridge, scheduler, scratch root and A100 resource availability have been
inspected. Shell syntax and wrapper arguments are checked locally. A new HPC GPU
pilot, environment installation and cross-platform numerical acceptance have not
yet been completed. Keep this distinct from the Colab run record.
