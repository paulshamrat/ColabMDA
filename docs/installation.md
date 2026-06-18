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

### 1.2. Storage Strategy: Local VM vs. Google Drive

When designing high-throughput Molecular Dynamics runs, we split file I/O into two zones:
* **Calculation Directory (`/content/work`):** Google Colab VMs are equipped with fast local NVMe SSDs mounted at `/content/`. Running the simulation step-by-step in `/content/work` prevents network lag, avoids API rate-limiting, and ensures maximum simulation speed.
* **Persistent Sync Directory (`/content/drive/MyDrive/ColabMDA/`):** Google Drive is mounted as a network share. At the completion of each trajectory chunk, ColabMDA automatically copies the completed `.dcd` and `.log` files, along with the latest state checkpoint (`prod.chk`), to your Google Drive folder. If the Colab instance shuts down or times out, your progress is safely preserved and ready to resume.

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

#### Step 3: Bootstrap the remote scientific environment
Execute the ColabMDA bootstrap script on the remote VM to install Miniforge (Conda), CUDA-enabled OpenMM, MDAnalysis, Modeller, and parameterization toolkits.
Run the bootstrap installer by piping python code to the remote kernel:

```bash
# Run the bootstrap setup remotely
echo "
import os
os.environ['WITH_MODELLER'] = '1'
os.environ['KEY_MODELLER'] = 'MODEL_LICENSE_KEY_HERE'
os.system('curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/pwpl/scripts/bootstrap_colab_openmm_gpu.sh | bash')
" | colab exec -s kras-sim
```

> 🔑 **Modeller License:** Replace `MODEL_LICENSE_KEY_HERE` with your license key. If you don't have one, register for a free academic license at [salilab.org/modeller/registration.html](https://salilab.org/modeller/registration.html).

---

### 1.4. Executing Simulation Commands

To run simulation commands inside the active virtual environment on the remote VM, you can open an interactive terminal:
```bash
# Launch a remote tmux terminal session
colab console -s kras-sim
```
Inside the console, activate the OpenMM environment and start your run:
```bash
conda activate openmm_env
cd /content/work
# Run your modular colabmda commands here
```

Or execute commands non-interactively using Python:
```bash
echo "import os; os.system('conda run -n openmm_env colabmda openmm --help')" | colab exec -s kras-sim
```

To shut down and release the VM when finished:
```bash
colab stop -s kras-sim
```

---

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
