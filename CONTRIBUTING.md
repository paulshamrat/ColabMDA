# Contributing to ColabMDA

Read the [main README](README.md) for installation, FKBP12/H-Ras workflows, Trp-cage
and the repository layout. Follow [AGENTS.md](AGENTS.md) for development rules.

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
  Do not duplicate scientific implementations or add extra example guide files.

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

Run CPU tests on Divergent. Protein GPU acceptance runs on Colab through the CLI,
using the generated example and the working project's `tests/gpu.py`. Keep it
foreground, inspect physical QC separately, retrieve verified compact reports,
and stop the dedicated session. Never include license keys or credentials in
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
