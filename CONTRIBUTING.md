# Contributing to ColabMDA

Read the [main README](README.md) for installation, the runnable Trp-cage example
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
- `examples/trpcage/scripts/` contains thin, runnable user workflows that call the
  public CLI/API. Install them with `cmda stage --example trpcage --root ...`.
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

Update README/AGENTS.md and affected reference pages when behavior changes.
Keep user instructions in the main README. Scientific protocol and advanced
CLI/visualization reference belong in docs/protocol.md and docs/reference.md;
do not recreate duplicate installation or workflow pages. Build Sphinx with
warnings treated as errors after documentation changes. Document only implemented commands;
identify compatibility-only commands and limits explicitly. Record changes in
the unreleased changelog without inventing a public release.
