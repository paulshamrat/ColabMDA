

# Modeller Workflow (Updated Protocol)

This workflow enables automated homology modeling of proteins using the Modeller software, designed for use in Colab or Linux environments. It streamlines the process of building 3D protein structures from sequence and template information, including environment setup, license configuration, model building, and result analysis. Researchers can utilize this protocol for tasks such as protein structure prediction, comparative modeling, and downstream analysis of model quality and structural features. The workflow is suitable for both educational and research purposes, supporting reproducible and efficient modeling projects.

This folder provides an updated, streamlined workflow for homology modeling with Modeller in Colab or Linux environments. The protocol now uses three main scripts:

- `install_modeller.sh`: Automated installer for Modeller and environment setup (conda/mamba).
- `modeller4.py`: Main homology modeling pipeline (CLI, logging, alignment, model building).
- `modeller_analysis.py`: Summarizes and analyzes Modeller output (metrics, RMSD, DOPE profile).


## Workflow Steps


### 1. Make Scripts Executable

```bash
chmod +x install_modeller.sh
chmod +x modeller4.py
chmod +x modeller_analysis.py
```

### 2. Install Modeller

Run the installation script in your Colab or terminal:

```bash
./install_modeller.sh
```

This sets up a dedicated conda environment (`modeller_env`) and installs Modeller.


### 3. Configure License and Activate Environment

You must obtain a valid Modeller license key from the official Modeller website:
- Visit: https://salilab.org/modeller/registration.html
- Register for a free academic license or purchase a commercial license as needed.
- The license key will be sent to your email address.

After installation, run the following commands to activate the environment and set the license key:

```bash
source ~/miniforge/etc/profile.d/conda.sh
conda activate modeller_env
# Replace YOUR_LICENSE_KEY with the license key you received from Modeller
LICENSE_KEY="YOUR_LICENSE_KEY"
CONFIG="$HOME/miniforge/envs/modeller_env/lib/modeller-10.7/modlib/modeller/config.py"
sed -i "s/^license *=.*/license = '${LICENSE_KEY}'/" "$CONFIG"
grep "^license" "$CONFIG"
python -c "import modeller; print('Modeller OK, version', modeller.__version__)"
```

This adds your license key to the Modeller configuration file (`config.py`) and confirms that Modeller is installed and licensed correctly.

### 4. Install Biopython

Before running the Modeller script, install Biopython (required for many modeling tasks):

```bash
pip install biopython
```


### 5. Run Homology Modeling Pipeline

Use the main pipeline script to build a model:

```bash
python3 modeller4.py <PDB_ID> <UNIPROT_ID>
```

**Example:**
```bash
python3 modeller4.py 4ldj P01116
```

This will:
- Create a folder named after the PDB ID
- Download PDB and UniProt FASTA files
- Clean the PDB, extract template sequence
- Write alignment.ali
- Run Modeller automodel (output in pipeline.log)
- Insert CRYST1 record into final PDB


### 6. Analyze and Summarize Results

Use the analysis script to compute metrics and plot DOPE profile:

```bash
python3 modeller_analysis.py --dir <PDB_ID>
```

**Example:**
```bash
python3 modeller_analysis.py --dir 4ldj
```

This will:
- Parse pipeline.log for molpdf, DOPE, GA341
- Compute sequence identity and coverage from alignment.ali
- Superimpose CA atoms and calculate RMSD
- Plot per-residue DOPE profile (saved as dope_profile.png)

## Files
- `install_modeller.sh`: Installs Miniforge, creates conda env, installs Modeller, sets license.
- `modeller4.py`: Homology modeling pipeline (CLI, logging, alignment, model building).
- `modeller_analysis.py`: Summarizes results, computes metrics, plots DOPE profile.




## Directory Structure: `modeller/`

```
modeller/
├── install_modeller.sh
├── modeller4.py
├── modeller_analysis.py
├── readme_modeller.md
├── obsolete/
│   └── run_modeller2.py
└── 4ldj/
    ├── 4ldj_clean.pdb
    ├── 4ldj_truncated_1-170.pdb
    ├── dope_profile.png
    ├── target.B99990001.pdb
    ├── target.B99990001_with_cryst.pdb
    ├── target.D00000001
    ├── target.ini
    ├── target.rsr
    ├── target.sch
    ├── target.V99990001
    ├── 4ldj_orig.pdb
    ├── alignment.ali
    ├── P01116.fasta
    ├── cryst1.txt
    ├── pipeline.log
    └── 4ldj.pdb
```

## Notes
- Edit `install_modeller.sh` to use your own valid Modeller license key.
- All output and logs are saved in the folder named after your PDB ID.
- This protocol is designed for reproducibility and ease of use in Colab or Linux.

---
For questions or issues, please contact the repository maintainer.
