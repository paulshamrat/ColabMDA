# Google Colab

Use this platform setup with [apo](../tutorials/apo.md) or [ATP/Mg](../tutorials/atpmg.md).
The scientific commands and data layout are shared with HPC. Use your local
terminal to manage runtimes with the Colab CLI and its remote console to run
simulation commands.

## Setup

Use the checkout whose instructions you are reading. Public releases use `main`;
private development stays separate. The tool stays outside the project where
simulations run.

### Terminal and Colab CLI

Prerequisites are Linux/Bash, Git, Python 3.12 or newer for the current Colab CLI,
internet access and a Colab account. The core ColabMDA package supports Python
3.10+, but the local Colab control environment must also meet the CLI requirement.
Complete CLI authentication when requested. Existing installations can activate
their environment instead of creating another one.

**On your local computer**, clone the public tool and enter its folder.
A GPU is not required on the local computer:

```bash
git clone --branch main https://github.com/paulshamrat/ColabMDA.git
cd ColabMDA
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e . google-colab-cli
cmda --help
colab sessions
```

Use a clean, committed checkout before creating the GPU runtime. If you edit
runtime files, commit those changes before deployment. An idle runtime can expire; a replacement needs
fresh Drive authentication even when the previous session was mounted.

Then create a GPU runtime, authenticate Drive and deploy the committed tool:

```bash
colab new -s colabmda-example --gpu T4
colab drivemount --session colabmda-example
cmda colab setup -s colabmda-example --local-source --modeller
colab console -s colabmda-example
```

Local-source deployment requires a clean, committed Git checkout. Branches and
release tags are both supported.
It sends runtime files and example templates with `SNAPSHOT.json` hashes. Run
these commands serially and wait for setup to finish before entering the console.
Setup asks whether to install MODELLER; CDKL5 model building and mutation require it.

**Inside the Colab console**, activate the installed environment:

```bash
source /root/miniforge3/etc/profile.d/conda.sh
conda activate colabmda
export MPLBACKEND=Agg
cmda --help
```

The installed snapshot is `/content/ColabMDA`; persistent project data belongs under
`/content/drive/MyDrive/research/projects/`. Continue with the linked CDKL5 workflow after setup. You can also dispatch foreground shell commands from the local
terminal using `cmda colab exec -s colabmda-example --timeout 21600 'COMMAND'`.

### Pulling Results Locally (`cmda pull`)

Lightweight simulation deliverables (`01_prep/model.pdb`, `02_eq/eq.pdb`, `analysis/<variant>/02_eq/qc.png`, `analysis/<variant>/02_eq/qc.txt`, and `analysis/` figures) can be pulled from Google Drive to your local workstation at any time with a single command:

```bash
# On your local workstation:
cmda pull -s colabmda-example
```

Large binary trajectories (`.dcd`), restart checkpoints (`.chk`), and raw coordinate files are skipped by default, so deliverables download in seconds. To download full trajectories for external rendering, add `--all`:

```bash
cmda pull -s colabmda-example --all
```

### Crash Recovery & Checkpointing

Keep simulations in the foreground with live output. Sessions can still terminate.
Keep inputs, completed chunks and checkpoints on Drive. On a replacement runtime,
remount Drive and restore the recorded software before attempting recovery.
Checkpoint compatibility across different GPUs or software versions is not assumed.
The selected scientific guide explains restart behavior.

An exec connection error does not prove the remote calculation stopped. Check
the foreground console, process and log before retrying, so two writers cannot
target the same run. Use `colab console -s colabmda-example` for long foreground
commands. If the runtime connection fails after computation has ended,
`colab restart-kernel -s colabmda-example` can restore control; verify Drive is
still mounted before collection. Do not restart a kernel during an active job.

### MODELLER License

Model building and mutation require MODELLER and a valid key. Its registration page
is [salilab.org/modeller/registration.html](https://salilab.org/modeller/registration.html).
The installer preserves the optional question and hidden key prompt.

```bash
cmda colab setup -s colabmda-example --local-source --modeller
```

Automated setup accepts `KEY_MODELLER` or `MODELLER_LICENSE` through the process
environment. Do not place real keys in commands saved to Git, logs or provenance.
`--no-modeller` skips the optional remote question. MODELLER uses `modeller_env`;
the selected scientific guide activates it and installs the same tool snapshot.

---

Each remote exec call starts a new shell. Include environment activation, project
variables and `cd` in that call; shell state does not carry between invocations.
A quoted multiline shell block runs in one Bash process and propagates failure.

## Execution

Follow the selected guide in the foreground. Its production target is per replica.
Heavy data stays on Drive; retrieve and verify compact reports before stopping
the dedicated session. Restore the recorded snapshot and compatible GPU/software
when resuming. HPC scheduler wrappers are not used in Colab.
