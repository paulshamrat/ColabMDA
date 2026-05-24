# Installation & Setup Guide

This guide describes how to install and set up **ColabMDA** for both Google Colab (cloud setup) and local environments.

---

## 1. Quick Start (Google Colab Installation)

> 💡 **Terminal Access:** All bash commands should be run in the **Colab Terminal** (Open via the **⋮** menu -> **Terminal**).

### 1.1. Setup Colab Runtime & Drive
Before starting, ensure your environment is ready:
1. **Enable GPU:** Go to `Runtime` -> `Change runtime type` and select **T4 GPU** (or any active GPU).
2. **Verify GPU:** Run `nvidia-smi` in the **Colab Terminal** (⋮ → Terminal) to confirm GPU access.
3. **Mount Drive:** Click the **Folder icon** 📂 in the left sidebar, then click the **Drive icon** (Mount Drive) to mount Google Drive natively.

### 1.2. Environment & Package Installation (Required)
Run the following in the **Colab Terminal** (⋮ → Terminal). 
*Estimated time: ~3–5 minutes.*

```bash
# 1. Install the core scientific environment
cd /content
curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/pwpl/scripts/bootstrap_colab_openmm_gpu.sh -o bootstrap_colab_openmm_gpu.sh
WITH_MODELLER=1 bash bootstrap_colab_openmm_gpu.sh latest

# 2. Install ColabMDA package
python3 -m pip install --upgrade "git+https://github.com/paulshamrat/ColabMDA.git@pwpl"
```

> [!IMPORTANT]
> **Modeller License Prompt:** During Step 1, the script will **pause** and ask you to `Enter your Modeller License Key`. You must paste your key and press **Enter** to proceed. The installation will not complete without it.
>
> 🔑 **Get a Free License:** If you don't have one, register at [salilab.org/modeller/registration.html](https://salilab.org/modeller/registration.html) (Academic licenses are free and sent instantly via email).

---

### 💡 Tip: How to Resume After a Timeout
If your Google Colab session expires:
1. Re-run **Required Steps 1.2** to reinstall the environment.
2. Run the **exact same `colabmda openmm run` command** you used before.
3. The tool will automatically detect your `.chk` files and resume from where it left off.

---

## 2. Manual / Detailed Installation (Advanced)

> ⚠️ **Note:** This section is for local workstations, HPC, or advanced manual setups. For standard Google Colab runs, please use the **Quick Start (Section 1)** above.

<details>
<summary><b>🛠️ A. Manual Terminal Installation (Step-by-Step)</b></summary>

In the Colab Terminal (⋮ → Terminal), run each step one at a time:

```bash
# Step 1: Download & install Miniforge (Conda)
wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh && bash /tmp/miniforge.sh -b -p "$HOME/miniforge3"

# Step 2: Initialize Conda in this shell
export PATH="$HOME/miniforge3/bin:$PATH" && source "$HOME/miniforge3/etc/profile.d/conda.sh"

# Step 3: Install Mamba into the base environment
conda install -y -n base -c conda-forge mamba

# Step 4: Install CUDA-enabled OpenMM and OpenMMTools
mamba install -y -c conda-forge cudatoolkit=11.8 openmm openmmtools

# Step 5: Install PDBFixer (conda, fallback to pip)
conda install -y -c conda-forge pdbfixer || pip install pdbfixer

# Step 6: Install MDAnalysis, MDTraj, NumPy, Matplotlib, and Biopython
mamba install -y -c conda-forge mdanalysis mdtraj numpy matplotlib biopython

# Step 7: Verify installations
python3 -c "from openmm import Platform; print('OpenMM platforms:', [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())]); import MDAnalysis, mdtraj, Bio; print('MDAnalysis:', MDAnalysis.__version__, 'MDTraj:', mdtraj.__version__, 'Biopython:', Bio.__version__)"
```
</details>

<details>
<summary><b>📜 B. Alternative: Script-based Installation</b></summary>

```bash
cd /content
curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/pwpl/scripts/install_colabmda_release.sh -o install_colabmda_release.sh
bash install_colabmda_release.sh latest /content/colabmda
```
</details>

<details>
<summary><b>🧬 C. Modeller CPU Environment Setup</b></summary>

```bash
cd /content/drive/MyDrive/openmm/ColabMDA
bash envs/install_modeller_env.sh
```
</details>

<details>
<summary><b>💻 D. Local Workstation Setup (Laptop/Desktop)</b></summary>

Beyond the Cloud ☁️: ColabMDA works on any Linux system with an NVIDIA GPU. Use the provided `environment.yml` to create a production-ready environment:
```bash
mamba env create -f environment.yml
conda activate colabmda
```
</details>

<details>
<summary><b>🏢 E. HPC Usage (SLURM)</b></summary>

You can easily incorporate ColabMDA into SLURM batch scripts. Since it processes trajectories in chunks, it is highly efficient for long-running jobs on cluster partitions with time limits.
</details>
