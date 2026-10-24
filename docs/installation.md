# Installation reference

The [main README](../README.md#install) is the maintained installation guide. It
also contains the complete Trp-cage example, project creation, Google Colab setup,
restart instructions and expected outputs. Follow that guide for a first run.

## Installer and dependencies

`scripts/setup.sh` is the only installer. It reads the core dependency list from
`environment.yml`, provisions the `colabmda` environment, installs the package and
runs health checks. There is no separate `envs/` installer directory.

MODELLER is optional and uses `modeller_env` for Python/license compatibility.
Interactive setup asks whether to install it, then requests a key without echo.
Automated setup accepts `KEY_MODELLER` or `MODELLER_LICENSE` through the process
environment. The installer does not add a key to global shell startup files.

Remote `cmda colab setup --modeller` prompts locally if necessary; `--no-modeller`
skips the optional question. Private development deployment uses
`cmda colab setup --local-source` from a committed `dev1` checkout. It exports
runtime files and example templates with a SHA-256 snapshot record.

## Storage and recovery

Keep the installed tool separate from the execution project. The worked project
is `colabmda-example` on Divergent and Drive, with numbered scripts and runtime
`data/01_structure`, `data/02_simulation/<run>` and `logs/02_simulation/<run>`.
No project configuration YAML is required.

Store persistent Colab results under the mounted Drive project. Runtime-local
`/content/ColabMDA` holds the installed snapshot, not the persistent data. A new
session needs Drive authentication and the appropriate software before resuming;
do not assume binary checkpoints work across changed GPUs/software versions.
Use `colab --help` for the installed CLI's account/session controls.

Standard tests run on Divergent without a GPU. GPU acceptance is separate and
must not be inferred from installation health checks or a successful command exit.
