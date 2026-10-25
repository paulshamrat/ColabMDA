# Contributing to ColabMDA

Read the [main README](README.md) for the preserved complete baseline and repository
layout. Shared CDKL5 guides are [apo](docs/apo.md) and [ATP/Mg](docs/atpmg.md),
with [Colab](docs/platforms/colab.md) and [HPC](docs/platforms/hpc.md) setup. Follow [AGENTS.md](AGENTS.md) for development rules.

## Development

Use `dev1` for private research and active development. Public `main` is reserved
for approved, sanitized releases; do not merge private development history into it.
Keep changes focused and describe the resulting behavior and validation.

Install the core environment from `environment.yml`, activate `colabmda`, then
install the checkout with `python -m pip install -e .`. MODELLER is optional;
use the single `scripts/setup.sh` installer and its key prompt when needed.

## Architecture

- `src/colabmda/cli.py` parses commands and routes operations.
- `core/` owns project paths, example installation, state definitions and ions.
- `modeller/` owns modeling and mutation logic.
- `amber/` owns parameterization and trajectory processing.
- `openmm/` is flat: preparation, EM, NVT, NPT, MD, QC, merge, analysis and figures
  use separate operation modules. Do not recreate engine/modular subfolders.
- `remote/colab.py` owns committed snapshot deployment and Colab execution.
- `examples/<name>/scripts/` contains thin user workflows calling the public
  CLI/API. Install fkbp, hras or trpcage with `cmda stage --example NAME --root ...`.
  Do not duplicate scientific implementations across platform adapters.

Keep calculations separate from plotting. Preserve native author numbering and
validate coordinate identities. Write runtime outputs into ignored project data
and logs, never into tool source or example templates. Test overwrite protection,
checkpoint recovery and output validation rather than only successful process exit.

## Verification

Before committing Python changes, run Black, `ruff format --check`, `ruff check`,
Python compilation and `pytest tests/`. Run `bash -n` on changed shell scripts.
If Black multiprocessing stalls, check one file per process.

Generated example scripts use the tool's 100-character line length. A project
should check them with `ruff format --line-length 100 --check scripts tests`;
reformatting them differently invalidates their recorded template hashes.

Run CPU tests on Divergent. Use Palmetto compbio for current protein GPU pilots;
use Colab through its CLI for Colab-specific acceptance, including the working
project's `tests/gpu.py`. Keep Colab runs foreground, inspect physical QC separately,
retrieve verified compact reports and stop dedicated sessions. Never include license keys or credentials in
source, test fixtures using real secrets, logs or provenance.

Natural-protein acceptance uses the staged FKBP12 WT/F37V and H-Ras–GDP/Mg
WT/G12V scripts. Verify author-to-UniProt mapping, unchanged raw inputs, modeled
mutations, actual topology charges/composition and model-based analysis labels.
Keep CPU modeling/preparation results separate from real GPU simulation evidence.

Update README/AGENTS.md when behavior changes. Keep the established README
organization: platform setup, project architecture, full simulation workflow and
the small validation example. Scientific explanations and advanced command
reference belong in the appropriate README sections. Do not replace it with an
example-only quickstart or recreate a separate documentation site.

Validate guide code blocks and command arguments against the current CLI. State
which checks were executed and which GPU or cluster steps remain unexecuted.
Record changes in the unreleased changelog without inventing a public release.

## CDKL5 guides and platforms

Maintain the complete root README alongside docs/apo.md and docs/atpmg.md. Platform
setup belongs in docs/platforms/{colab,hpc}.md; common scientific commands stay in
the state guides. Use explicit cd, environment activation and direct cmda commands.
HPC alone needs the thin examples/hpc/simulate.slurm scheduler adapter.

The root walkthrough also follows CDKL5 WT/L119R for apo and ATP/one-Mg. Keep its
current design and section order. Use concise wording and matching state-based
paths throughout modeling, preparation, analysis and comparisons. Trp-cage stays
the optional quick test; retain older FKBP12/H-Ras templates and run records.

The default CDKL5 example must use public inputs, not files copied from another
study. `cmda model` delegates to the existing MODELLER build/mutate engine and
checks outputs. Public `cmda assemble --template` records homologous donor fitting
and optional fixed-backbone sidechain relaxation. A modeled CDK2-to-CDKL5 ATP/Mg
pose is not an experimental complex or reproduction of a legacy study. Keep
refinement, assembly and solvated equilibration as separate scientific steps.
Run short GPU pilots on Palmetto compbio when available; do not infer acceptance
from modeling/preparation or from the previous study-based pilots.

New runs use unnumbered state/structure/simulation folders and numbered variant
processing steps. Preserve existing numbered example paths. Require canonical
identity/density gates, frozen preparation inputs, real topology charges and
noncovalent boundaries. Research production requires complete 24-stage timing and
physical QC. Separate software tests, Colab GPU results and HPC pilot status.
Update AGENTS.md and the local ignored PLAN.md with material design/status changes.
