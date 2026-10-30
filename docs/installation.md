# Installation & Setup

ColabMDA divides work between your local computer and Google Colab:
* **Local machine (`cmda`):** Structure modeling, protonation, and PyMOL inspection.
* **Google Colab (`cmda`):** GPU molecular dynamics simulations.

| Environment | Where It Runs | What Gets Installed | When To Use |
| :--- | :--- | :--- | :--- |
| **Local (`cmda`)** | Linux, macOS, WSL2 | OpenMM, PyMOL, Modeller, AmberTools, `cmda` CLI | Phase 1 (Prep) & Phase 4 (Analysis) |
| **Cloud (`cmda`)** | Google Colab | OpenMM (CUDA), AmberTools, `cmda` CLI | Phase 3 (Simulation) |

---

## 1. Local Workstation Setup

Run this one-time setup on your local machine:

```bash
# 1. Run the setup script:
bash -c "$(curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/dev1/scripts/setup.sh)"

# 2. Activate environment:
source ~/miniforge3/etc/profile.d/conda.sh
conda activate cmda

# 3. Log in to Google Colab (one-time browser login):
colab sessions
```

> [!TIP]
> **Modeller Key:** Free for academic research ([register here](https://salilab.org/modeller/registration.html)). Enter it when prompted, or run:  
> `KEY_MODELLER="YOUR_KEY" bash -c "$(curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/dev1/scripts/setup.sh)"`

<details>
<summary><b>🛠️ Install from Source</b></summary>

```bash
git clone --branch dev1 https://github.com/paulshamrat/ColabMDA.git ~/ColabMDA
cd ~/ColabMDA
KEY_MODELLER="YOUR_KEY" bash scripts/setup.sh
```
</details>

<details>
<summary><b>🧪 Run Test Suite</b></summary>

```bash
conda activate cmda
pytest tests/ -v
```
</details>

---

## 2. Google Colab Setup

> 💡 **Where to run:** Inside your Colab terminal or `colab console` during **Phase 3**.

Run this command once inside Colab to install OpenMM and AmberTools (~90 seconds):

```bash
curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/dev1/scripts/setup.sh | bash
source /root/miniforge3/etc/profile.d/conda.sh
conda activate cmda
```

> [!NOTE]
> Modeller and PyMOL are automatically skipped on Google Colab (`WITH_MODELLER=0`). Setup runs completely unattended with no prompts.
