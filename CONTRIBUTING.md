# Contributing to ColabMDA

The [README](README.md) is the complete user guide. Scientific walkthroughs are available in [`docs/tutorials/`](docs/tutorials/) (Trp-Cage Quickstart, KRAS Apo, CDKL5 Apo, and CDKL5 Holo).

## Development

Keep changes focused and describe the behavior and validation. Public `main`
contains reviewed release snapshots. Keep private research history, credentials,
coordinates, trajectories and run records outside the public repository.

Create the scientific environment from `environment.yml`, activate it and install
the checkout with `python -m pip install -e .`. MODELLER is optional and requires
your own license; use the single `scripts/setup.sh` installer and hidden prompt.

Report reproducible problems through [GitHub Issues](https://github.com/paulshamrat/ColabMDA/issues).
Formatting and linting use [Black](https://github.com/psf/black) and
[Ruff](https://github.com/astral-sh/ruff).

## Credential protection

Install [Gitleaks](https://github.com/gitleaks/gitleaks#installing), then enable the
repository hooks:

```bash
git config core.hooksPath .githooks
gitleaks version
```

If the scanner is not on PATH, set its executable using the local Git setting
`colabmda.gitleaksPath`. The commit hook scans staged content. The push hook scans
all reachable history. Both redact findings and fail if the scanner is missing.
The configuration adds MODELLER license and Google session-cookie checks to
the standard credential rules. Notebook outputs must not contain captured
request headers or browser cookies.
Only explicit placeholders belong in examples; never add a real value to an
allowlist or scanner configuration.

Keep credentials in private files outside Git or use the installer's hidden
prompt. `.env`, credential files and logs are ignored; ignore rules do not remove
files already committed. GitHub Actions credentials belong in repository secrets.
If a value was published, replace/revoke it with its provider and clean affected
history. Coordinate replacement clones after a history rewrite.

For a private research checkout, set `colabmda.privateCheckout=true` and disable
its public push URL. Export reviewed code to a separate clean-history release
checkout; never merge the private branch into public `main`.

## Architecture

- `cli.py` parses commands; `core/` owns project paths, mapping and checks.
- `modeller/` owns modeling and mutation; `amber/` owns parameterization and merging.
- `openmm/` keeps preparation, EM, NVT, NPT, MD, QC, analysis and figures in flat,
  separate operation modules. Calculations never import plotting libraries.
- `remote/colab.py` deploys committed source snapshots and controls Colab execution.
- `examples/` contains thin workflows and Slurm adapters calling the public CLI/API.
  Keep one scientific implementation across platforms.

Preserve native numbering/insertion codes and check identities and coordinate
presence. Raw structures stay unchanged. Runtime data, reports and logs belong
in the execution project. New state workflows keep structure, simulation and
analysis as unnumbered siblings; only simulation processing steps are numbered.
Replica reports, aggregates and variant comparisons remain separate. Historical
numbered examples retain their paths and template hashes.

Check noncovalent component boundaries, actual charges and complete water-matched
ion parameters. Research production requires full 24-stage timing and physical
QC. Quick smoke tests and completion markers do not establish convergence.
Advance checkpoints only after complete chunks are published; recovery must
preserve continuous frame times. Distinguish model/preparation tests from GPU runs.

CDKL5 public inputs are 4BGQ/O76039 and homologous ATP/Mg transfer from 8FP5. Keep
its modeled pose distinct from an experimental complex or previous study inputs.
Visualization fits native core positions to a common WT reference and uses shared
framing. Frame labels come from matching logs, and model panels remain distinct
from trajectory panels.

## Verification

Before committing Python changes, run formatting, lint, compilation and CPU tests:

```bash
python -m pip install pytest build black ruff
black --check .
ruff format --check src tests
ruff check .
python -m compileall -q src tests
python -m pytest tests/ -q
```

If Black multiprocessing stalls, check one file per process. Run `bash -n` on
changed shell/Slurm scripts. Generated example scripts use 100-character lines;
changing their formatting also changes recorded template hashes.

Tests cover public input gates, overwrite/symlink protection, replica isolation,
checkpoint recovery, wheel resources, committed source export and hidden license
prompting. Run credential-hook tests with Gitleaks on PATH. GPU pilots belong on
an allocated cluster GPU or Colab; never claim GPU acceptance from CPU tests.
Retrieve verified compact reports and review plots before ending a pilot session.

## Documentation

Preserve the README's existing section order, protocol and timestep/frame/storage
calculations. Use simple, direct, concise language. Keep all seven badges directly
below the title, followed by Colab interruption recovery and the starting path.
Platform/dependency details belong in Platform Setup, architecture in section 2,
and one compact command reference near the end.

Shared scientific guides live in docs/tutorials/ (trpcage_apo.md, kras_apo.md, cdkl5_apo.md, cdkl5_atpmg.md). Use terminal commands for Colab; do not recreate notebook
instructions or a separate documentation site. Check local links, shell blocks
and CLI arguments after guide changes. Keep account-specific deployment records
in the private execution project. Record unreleased changes without announcing
a release or claiming unperformed simulation tests.
