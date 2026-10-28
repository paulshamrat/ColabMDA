# Installation & Environment Setup

ColabMDA uses a **Hybrid Architecture** designed to preserve free cloud compute units and avoid session timeouts. To keep workflows simple and consistent, the conda environment is named **`cmda`** across both local and cloud environments:

| Environment | Where It Runs | What Gets Installed | When To Install & Use |
| :--- | :--- | :--- | :--- |
| **1. Local Environment (`cmda`)** | Your Laptop / Workstation (macOS, Linux, WSL) | Python, MODELLER, Open-Source PyMOL, AmberTools (`tleap`), MDAnalysis, `google-colab-cli`, ColabMDA CLI | **Do this right now (one-time setup).** Used for Phase 1 (Modeling, Structure Prep, Solvation, PyMOL QA, Cloud Management) and Phase 4 (Inspection). |
| **2. Cloud GPU Environment (`cmda`)** | Google Colab Cloud VM (NVIDIA T4 / A100 / L4) | Miniforge, OpenMM (CUDA), AmberTools `tleap`, MDAnalysis, ColabMDA CLI | **Do this in Phase 3** (only after launching your GPU instance). Runs 100% unattended in ~90 seconds with zero license prompts. |

---

## 1. Local Workstation Installation

A single bootstrap script sets up Miniforge, the `cmda` conda environment, Modeller, PyMOL, AmberTools, MDAnalysis, and `google-colab-cli` automatically—**without requiring you to clone the git repository**:

```bash
# 1. Run the universal bootstrap script (runs in-memory, leaves no files behind):
bash -c "$(curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/dev1/scripts/setup.sh)"

# 2. Activate the environment
source ~/miniforge3/etc/profile.d/conda.sh
conda activate cmda

# 3. Authenticate Google Colab CLI (one-time browser login)
colab sessions
```

> [!NOTE]
> **Modeller Academic License:** Modeller is completely free for academic research ([registration](https://salilab.org/modeller/registration.html)). During setup, the installer will pause and prompt you to enter your license key, or you can provide it beforehand:
> ```bash
> KEY_MODELLER="YOUR_ACADEMIC_KEY" bash -c "$(curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/dev1/scripts/setup.sh)"
> ```

<details>
<summary><b>🛠️ Developer Alternative: Clone from Source</b></summary>

If you are developing new features or contributing code to ColabMDA, clone and install locally:
```bash
git clone --branch dev1 https://github.com/paulshamrat/ColabMDA.git ~/ColabMDA
cd ~/ColabMDA
KEY_MODELLER="YOUR_ACADEMIC_KEY" bash scripts/setup.sh
```
</details>

<details>
<summary><b>🧪 Running the Test Suite Locally</b></summary>

Verify code integrity, force field mechanics, and CLI parsing (`cmda` env required):
```bash
conda activate cmda
pytest tests/ -v
```
</details>

---

## 2. Cloud GPU Worker Provisioning (Google Colab)

> 💡 **Where to run this:** Do not run this locally. Execute this inside your Google Colab instance (via `colab console` or the web terminal) during **Phase 3** of your workflow after launching a GPU runtime.

Google Colab instances are ephemeral (they spin up on demand). When you launch your GPU instance in Phase 3, the exact same bootstrap script provisions the cloud environment in ~90 seconds:

```bash
# Provision Cloud GPU Worker (executed inside Colab console):
curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/dev1/scripts/setup.sh | bash
source /root/miniforge3/etc/profile.d/conda.sh
conda activate cmda
```

> [!NOTE]
> **Zero License Prompts on Cloud:**  
> Because homology modeling and coordinate building are already completed locally on your workstation, MODELLER and PyMOL are automatically excluded on Colab (`WITH_MODELLER=0`). Google Colab only installs Miniforge, OpenMM (CUDA), AmberTools `tleap`, MDAnalysis, and ColabMDA. Setup runs **100% unattended in ~90 seconds with zero interactive prompts**.
