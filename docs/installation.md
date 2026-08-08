# Installation & Setup Guide

This guide describes how to install, authenticate, and set up **ColabMDA** using the **Google Colab CLI** workflow ("The Colab CLI Way") as well as local workstation setups.

---

## 1. The Colab CLI Workflow (Recommended)

Using the `colab` command-line interface allows you to provision, configure, and monitor GPU-accelerated simulation runtimes directly from your local workstation terminal, without needing to open a web browser.

### 1.1. Local Installation & Authentication

First, install the Google Colab CLI on your local workstation/laptop:

```bash
# Install the Colab CLI
pip install google-colab-cli
```

To authenticate the CLI with your Google account:
1. Run a command that queries the backend, such as:
   ```bash
   colab sessions
   ```
2. The CLI will output an authentication URL.
3. Open the link in your web browser, log in with your Google account, and grant the required permissions.
4. Copy the authorization code provided and paste it back into your terminal.
5. The OAuth token is saved locally and will be automatically reused.

#### 🔄 Switching Google Accounts (Starting Fresh)
If you need to switch to a different Google account, you must clear your locally cached OAuth token and session records:
```bash
rm -f ~/.config/colab-cli/token.json ~/.config/colab-cli/sessions.json
```
On your next command, the Colab CLI will trigger a fresh login flow and prompt you to authorize using your new account credentials.

---

### 1.2. Storage Strategy: Google Drive as the Primary Working Directory

All simulation data is read from and written directly to your **mounted Google Drive**. This is critical — anything written to `/content/` (the local VM disk) is **lost permanently** when your Colab session ends or times out.

```text
/content/drive/MyDrive/ColabMDA/         <- your persistent workspace
  data/
    str/        <- structure templates (cleaned PDBs)
    sim/        <- simulation output folders (trajectories, checkpoints)
    analysis/   <- per-system analysis plots and metrics
```

This directory layout is **gitignored** in the repository so that large trajectory files (`.dcd`, `.chk`) are never accidentally committed to Git.

---

### 1.3. Provisioning and Bootstrapping a Remote GPU VM

Run the following commands on your **local workstation** to set up a remote GPU simulation node:

#### Step 1: Allocate a new GPU VM session
```bash
# Create a new session named 'kras-sim' with a T4 GPU
colab new -s kras-sim --gpu T4
```

#### Step 2: Mount your Google Drive
```bash
# Mounts your Drive at /content/drive/ on the remote VM
colab drivemount -s kras-sim
```

> [!WARNING]
> **Consent Screen Requirement:** When the browser window opens for authorization, you **must check the box** authorizing Google Colab to **"See, edit, create, and delete all your Google Drive files"**. If this box is left unchecked, Google will authorize the credentials but block Colab from accessing the filesystem, causing the command to hang or fail.

#### Step 3: Connect to the remote console
Open an interactive console session to run shell commands directly on the remote VM:
```bash
# Launch a remote console session
colab console -s kras-sim
```

#### Step 4: Bootstrap the remote scientific environment
Inside the remote console, run the bootstrap installer to install Miniforge (Conda), CUDA-enabled OpenMM, PDB2PQR, PROPKA, MDAnalysis, Modeller, and parameterization toolkits:

```bash
cd /content
curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/pwpl/scripts/bootstrap_colab_openmm_gpu.sh -o bootstrap_colab_openmm_gpu.sh
chmod +x bootstrap_colab_openmm_gpu.sh
export WITH_MODELLER=1
export KEY_MODELLER=MODEL_LICENSE_KEY_HERE
bash ./bootstrap_colab_openmm_gpu.sh 2>&1 | tee /content/bootstrap_colabmda.log
```

> 🔑 **Modeller License:** Replace `MODEL_LICENSE_KEY_HERE` with your license key. If you don't have one, register for a free academic license at [salilab.org/modeller/registration.html](https://salilab.org/modeller/registration.html).
>
> 🛠️ **Developer Tip:** If you are testing custom or unmerged changes from a specific git branch, you can replace `pwpl` in the `curl` URL above with your custom branch name (e.g. `main` or a feature branch name).
>
> 🧬 **Force-field note:** The standard examples use Amber19/ff19SB with OPC water, so the bootstrap installs a modern OpenMM build with Amber19 XML support. Colab T4 runtimes currently report CUDA 13.0; if your runtime requires a different CUDA package set, override it with `COLABMDA_CUDA_VERSION=12.9` or another compatible value before running the bootstrap.
>
> 🎯 **Docking option:** If you also want virtual-screening tools, add `export WITH_DOCKING=1` before running the bootstrap. This installs RDKit, Meeko, and AutoDock Vina into a separate `docking_env`.

If you want to install and test your **current local checkout** instead of a pushed GitHub branch, run this from your local workstation:

```bash
python scripts/colab_exec_bootstrap.py --session kras-sim --from-local --timeout 7200
```

For the full MODELLER install from your current local checkout:

```bash
export KEY_MODELLER='YOUR_MODELLER_LICENSE_KEY'
python scripts/colab_exec_bootstrap.py --session kras-sim --from-local --with-modeller --timeout 7200
```

For local checkout testing with docking tools:

```bash
python scripts/colab_exec_bootstrap.py --session kras-sim --from-local --with-docking --with-ligand --timeout 7200
```

This uploads a lightweight source zip to `/content/ColabMDA-upload.zip`, installs the local checkout editable at `/content/ColabMDA`, streams the transcript locally, and saves the VM log at `/content/bootstrap_colabmda.log`.

If the Colab Terminal appears to hang or prints only progress-bar lines, open a second terminal and check whether the installer is still active:

```bash
ps -eo pid,stat,etime,cmd | grep -E 'bootstrap|mamba|conda|wget|curl|python' | grep -v grep
tail -n 40 /content/bootstrap_colabmda.log
du -sh "$HOME/miniforge3" 2>/dev/null || true
```

During `mamba create` or `conda install`, Colab Terminal may redraw package progress bars poorly. If the `ps` command still shows `mamba` or `conda`, let it continue. If there is no installer process and `conda` is still unavailable, rerun from the saved script:

```bash
cd /content
bash ./bootstrap_colab_openmm_gpu.sh 2>&1 | tee -a /content/bootstrap_colabmda.log
```

### 1.4. Go to Your Drive Workspace

Once the bootstrap script completes, your scientific environments are ready. **All simulation commands should target your persistent Google Drive path** so that data is never lost when the session ends.

Inside the remote console, go to your Drive project folder. For this branch, the default examples use `/content/drive/MyDrive/ColabMDA`. You can choose another folder if your lab or project uses a different layout; just `cd` into that folder before running ColabMDA commands.

```bash
cd /content/drive/MyDrive/ColabMDA
mkdir -p data/{str,sim,analysis}
```

From this point onward, use relative paths such as `data/str/...` and `data/sim/...`, matching the main ColabMDA command style.

---

### 1.5. Executing Simulation Commands

Inside the remote console:
```bash
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate openmm_env

# Stage a new simulation directly on your Drive
colabmda openmm stage --pdb-file data/str/kras_wt/kras_wt.pdb --name kras_wt --replica r1

# Run energy minimization
colabmda openmm em --name kras_wt --workdir data/sim/kras_wt/r1

# Run production MD (trajectories go to Drive automatically)
colabmda openmm md --name kras_wt --workdir data/sim/kras_wt/r1
```

> [!IMPORTANT]
> Run these commands from your mounted Drive workspace, for example `/content/drive/MyDrive/ColabMDA`. Do not run simulations from `/content/ColabMDA`; that folder is only for the temporary editable code install and will disappear when the Colab VM ends.

To shut down and release the VM when your simulation is finished:
```bash
# Run this from your local workstation terminal
colab stop -s kras-sim
```

---

### 1.6. Syncing Drive Data Back to Your Local Workstation

After your simulations complete, sync data from Google Drive back to your local `data/` folder for visualization and analysis.

#### Option A: Selective Sync (Recommended — Analysis & Plots Only)
```bash
# Pull only lightweight analysis results, plots, and reports
rclone copy gdrive:ColabMDA/data/ ~/works/ColabMDA/data/ \
  --include "**/analysis/**" \
  --include "**/*.png" \
  --progress

```

#### Option B: Full Sync (Trajectories Included)
```bash
rclone copy gdrive:ColabMDA/data/ ~/works/ColabMDA/data/ --progress
```

> [!WARNING]
> **Bandwidth Warning:** Molecular dynamics trajectory files (`.dcd`) are very large (often 10–50 GB per replica). Only use Option B if you have a fast internet connection and sufficient local disk space. For routine figure generation and analysis, Option A (selective sync) is strongly recommended.

## 2. Local Workstation Setup (Laptop/Desktop)

Beyond the Cloud ☁️: ColabMDA works on any local Linux system with an NVIDIA GPU.

### 2.1. Environment Setup
Use the provided `environment.yml` to create a production-ready conda environment:
```bash
mamba env create -f environment.yml
conda activate colabmda
```

### 2.2. Modeller Environment Setup
To build homology models locally, install Modeller into a separate environment:
```bash
bash envs/install_modeller_env.sh
```

---

## 3. HPC Usage (SLURM)

You can easily incorporate ColabMDA into SLURM batch scripts. Since it processes trajectories in chunks, it is highly efficient for long-running jobs on cluster partitions with time limits. Refer to the [Protein-in-Water Guide](protein_water_simulation.md) for job scripts.

---

## 4. Building and Viewing Documentation Locally

To compile the user manual and guides into local HTML pages that you can browse in your web browser:

1. **Build the HTML site** using the `colabmda_test` conda environment:
   ```bash
   conda run -n colabmda_test sphinx-build -b html docs/ docs/_build/html
   ```

2. **Open the docs locally** in your browser:
   * **Direct File Open**:
     ```bash
     xdg-open docs/_build/html/index.html
     ```
   * **Development Web Server**:
     ```bash
     python -m http.server --directory docs/_build/html 8000
     ```
     Then navigate to [http://localhost:8000](http://localhost:8000) in your web browser.
