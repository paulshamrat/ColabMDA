# ColabMDA

ColabMDA prepares protein systems, runs checkpointed OpenMM simulations, processes
trajectories with AmberTools, and produces numerical analysis and figures.
MODELLER is optional for modeling and mutation workflows.

## Repository layout

```text
ColabMDA/
├── README.md                      # Installation and complete example guide
├── AGENTS.md                      # Development rules
├── pyproject.toml                 # Python package
├── environment.yml                # Core simulation dependencies
├── scripts/setup.sh               # Single installer
├── src/colabmda/
│   ├── cli.py
│   ├── core/                      # project.py, examples.py, taxonomy.py, ions.py
│   ├── modeller/                  # engine.py, substrate.py
│   ├── amber/                     # leap.py, cpptraj.py, params/
│   ├── openmm/
│   │   ├── commands.py
│   │   ├── prepare.py
│   │   ├── em.py
│   │   ├── nvt.py
│   │   ├── npt.py
│   │   ├── md.py
│   │   ├── qc.py
│   │   ├── merge.py
│   │   ├── analysis.py
│   │   ├── figures.py
│   │   ├── common.py
│   │   └── run.py
│   └── remote/colab.py
├── examples/trpcage/scripts/       # Runnable templates, no data or extra guides
├── tests/
└── docs/                          # Detailed scientific/API reference
```

`src/` implements the tool. `examples/` contains workflows copied into a user's
project. Project `data/` and `logs/` hold runtime products and stay outside Git.
The working example project is consistently called **colabmda-example** locally
and on Drive. Existing research projects retain their own names.

## Install

On a Linux workstation with a supported NVIDIA GPU:

```bash
git clone https://github.com/paulshamrat/ColabMDA.git
cd ColabMDA
bash scripts/setup.sh
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate colabmda
```

The installer reads `environment.yml`, installs the package and checks the runtime.
It asks whether to install MODELLER and, if selected, prompts for a license key
without displaying it. The key is configured in its separate `modeller_env`,
not written to repository files or the user's global shell startup file.
MODELLER is not required for the Trp-cage example.

For an existing Conda installation, the core environment can also be created
with `mamba env create -f environment.yml`, followed by `conda activate colabmda`
and `python -m pip install -e .`. Optional modeling installation still uses the
single installer: `WITH_MODELLER=1 bash scripts/setup.sh`.

For automated setup, supply the key through `KEY_MODELLER` or `MODELLER_LICENSE`
in the process environment. Do not put a real key in command examples, logs or
Git. Use `cmda colab setup --modeller` to request it remotely; the CLI prompts
locally if needed and transfers it through a temporary credential file rather
than embedding the key in remote command text. `--no-modeller` skips that question.

## Create the complete example

Run this where ColabMDA is installed:

```bash
cmda stage --example trpcage --root ~/research/projects/colabmda-example
cd ~/research/projects/colabmda-example
```

The command copies seven small scripts, records their hashes and tool version,
and refuses to overwrite different existing files. Repeating it with identical
scripts is safe. It does not download structures, install software or run MD.
No project configuration YAML is created. Settings are explicit in the scripts.

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
│   ├── 01_structure/trpcage/       # Created by fetch
│   └── 02_simulation/<run>/        # Created during execution
└── logs/02_simulation/<run>/
```

## Run the example

Use a CUDA GPU machine. Divergent runs standard tests, not protein GPU pilots.
Activate the environment, enter the generated project and run the scripts in order:

```bash
python scripts/01_structure/01_fetch.py
export COLABMDA_RUN_ID=run01
bash scripts/02_simulation/01_prepare.sh
bash scripts/02_simulation/02_equilibrate.sh
bash scripts/02_simulation/03_production.sh
bash scripts/02_simulation/04_analysis.sh
python scripts/02_simulation/05_check.py
```

Fetch retrieves [Trp-cage, PDB 1L2Y](https://www.rcsb.org/structure/1L2Y), preferring
mmCIF with a PDB fallback. It selects NMR model 1, author chain A, checks all 20
residues against the deposited sequence, verifies C-alpha atoms and preserves
native numbering. The raw archive and source/input hashes remain in project data.
This designed miniprotein has no canonical UniProt mutation annotation.

The default uses Amber ff19SB/OPC preparation, quick 100 ps NVT and 100 ps NPT,
then 1 ns production at 2 fs. Frames are saved every 10 ps in 100 ps checkpoint
chunks. This is a software example, not a production equilibration recommendation.
Our earlier T4 test took about 7.5 minutes for production plus the restart test;
measure performance on your own runtime. A quick density QC warning remains a
warning even if the software checks pass.

To change duration, edit the visible setting in `03_production.sh`, or set
`COLABMDA_EXAMPLE_NS` for both production and check. Use a positive duration
in multiples of 0.1 ns. Use a new run ID for changed inputs or protocols.

Keep execution logs, for example:

```bash
set -o pipefail
bash scripts/02_simulation/03_production.sh 2>&1 | tee "logs/02_simulation/$COLABMDA_RUN_ID/production.log"
```

## Interrupt and resume

Press Ctrl+C during production, then rerun the same production command with the
same project, run ID, inputs, software and settings. An interrupted command must
report failure; it resumes the last committed chunk checkpoint and replays the
unfinished chunk. Do not delete chunks to force recovery or assume binary
checkpoints are portable between different GPUs or software versions.

A short smoke test can be extended without changing the underlying protocol:

```bash
COLABMDA_EXAMPLE_NS=0.1 bash scripts/02_simulation/03_production.sh
COLABMDA_EXAMPLE_NS=1 bash scripts/02_simulation/03_production.sh
```

The private working project's `tests/gpu.py` automates actual SIGINT, unchanged
chunk hashes, resume, analysis and final checks using these same generated scripts.

## Results and checks

Each run contains `provenance.json`, production request history and
`apo/WT/{01_prep,02_eq,03_md/r1,04_analysis}`. Provenance records the input,
installed template revision, executed scripts, runtime package versions and GPU
identity.

`04_analysis/` contains RMSD, Rg and RMSF CSVs plus `figures/` PNGs. Native author
residue labels are verified against the model; direct Amber analysis through
`cmda openmm analysis` requires `--reference model.pdb`. `analysis.json` records
its reference and implementation fingerprints. Computation and figure rendering
are separate modules. Figures use 7-inch width, 300 DPI and publication typography.

The final check requires continuous frame times, finite coordinates/logs/tables,
20 C-alpha atoms, RMSF author labels 1–20 and nonempty figures. Only then does it
write `summary.json` and `RUN_DONE`. Physical equilibration QC is reported
separately. The default complete run has ten chunks and 100 merged frames.

To rerender from saved numerical outputs without MD, use the public figure API:

```bash
python -c 'import os; from pathlib import Path; from colabmda.openmm.figures import plot_analysis; plot_analysis(Path("data/02_simulation")/os.environ["COLABMDA_RUN_ID"]/"apo/WT/04_analysis")'
```

## Google Colab

From Divergent, create a dedicated GPU session and authenticate Drive:

```bash
colab new -s colabmda-example --gpu T4
colab drivemount --session colabmda-example
cmda colab setup -s colabmda-example --local-source
colab console -s colabmda-example
```

For private development, local-source setup exports committed `dev1` runtime
files and hashes them in `SNAPSHOT.json`; commit changes before deployment.
Inside Colab, activate `colabmda`, then create the example directly on Drive:

```bash
source /root/miniforge3/etc/profile.d/conda.sh
conda activate colabmda
cmda stage --example trpcage --root /content/drive/MyDrive/research/projects/colabmda-example
cd /content/drive/MyDrive/research/projects/colabmda-example
```

Follow the same Run section above. Keep jobs in the foreground with live output.
Alternatively, dispatch foreground commands from Divergent through
`cmda colab exec -s colabmda-example --timeout 21600 'YOUR_COMMAND'`.
Sessions can terminate; retain checkpoints on Drive for recovery. Retrieve and
verify compact reports and plots, leaving trajectories on Drive, then stop the
session with `colab stop -s colabmda-example`.

## Development checks

Use `dev1` for private work. Public `main` contains only approved releases.
Run CPU unit/regression tests on Divergent using `pytest tests/`, plus Black,
Ruff, Python compilation and `bash -n` for changed scripts. GPU integration tests
are separate. Use `cmda --help` for modeling and other simulation commands and
[docs/](docs/) for detailed reference.
