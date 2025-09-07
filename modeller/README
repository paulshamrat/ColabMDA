# Modeller Homology Modeling & Mutation Pipeline

This folder provides a simple, unified workflow for homology modeling and mutation of proteins using Modeller in Colab or Linux environments. The protocol is designed for researchers and students who want to build 3D protein models from sequence and template, and optionally introduce mutations, all from the command line.

---

## 1. Download Required Scripts

Download the main scripts (no analysis script needed for basic modeling):

```bash
wget https://raw.githubusercontent.com/paulshamrat/ColabMDA/main/modeller/install_modeller.sh
wget https://raw.githubusercontent.com/paulshamrat/ColabMDA/main/modeller/modeller6.py
chmod +x install_modeller.sh modeller6.py
```

---

## 2. Install Modeller

Run the installation script (creates a conda environment and installs Modeller):

```bash
./install_modeller.sh
```

---

## 3. Configure Modeller License

Obtain a free academic license from https://salilab.org/modeller/registration.html
After installation, activate the environment and set your license key:

```bash
source ~/miniforge/etc/profile.d/conda.sh
conda activate modeller_env
# Replace YOUR_LICENSE_KEY with your actual key
LICENSE_KEY="YOUR_LICENSE_KEY"
CONFIG="$HOME/miniforge/envs/modeller_env/lib/modeller-10.7/modlib/modeller/config.py"
sed -i "s/^license *=.*/license = '${LICENSE_KEY}'/" "$CONFIG"
grep "^license" "$CONFIG"
python -c "import modeller; print('Modeller OK, version', modeller.__version__)"
```

---

## 4. Install Biopython

Biopython is required for many modeling tasks:

```bash
pip install biopython
```

---

## 5. Usage: modeller6.py (Command Line)

`modeller6.py` is a unified script for building homology models and introducing mutations (single or batch) with flexible options for chain, residue range, truncation, and output. All major steps are handled from the command line.

### Build a model (default mode)
```bash
python3 modeller6.py <PDB_ID> <UNIPROT_ID> [--chain A] [--range START END] [--truncate] [--mut K76E] [--list mutations.txt] [--outdir DIR] [--outdir-mut DIR]
```

**Examples:**
- Basic build:
  ```bash
  python3 modeller6.py 4bgq O76039
  ```
- Specify chain and range:
  ```bash
  python3 modeller6.py 4bgq O76039 --chain A --range 1 303
  ```
- Build and mutate:
  ```bash
  python3 modeller6.py 4bgq O76039 --chain A --range 1 303 --mut K76E
  ```
- Batch mutations:
  ```bash
  python3 modeller6.py 4bgq O76039 --chain A --range 1 303 --list mutations.txt --outdir-mut 4bgq_fix/mutants
  ```
- Truncate after build:
  ```bash
  python3 modeller6.py 4bgq O76039 --chain A --range 1 303 --truncate
  ```

### Mutate-only mode (mutate an existing PDB)
```bash
python3 modeller6.py --pdb-in <PDB_FILE> --chain A --mut K76E [--list mutations.txt] [--outdir-mut DIR]
```

**Examples:**
- Single mutation:
  ```bash
  python3 modeller6.py --pdb-in 4bgq_fix/target.B99990001_with_cryst.pdb --chain A --mut K76E
  ```
- Batch mutations:
  ```bash
  python3 modeller6.py --pdb-in 4bgq_fix/target.B99990001_with_cryst.pdb --chain A --list mutations.txt --outdir-mut 4bgq_fix/mutants
  ```

### Options
- `--chain` : Specify chain (default: A)
- `--range START END` : UniProt residue range (inclusive)
- `--truncate` : Truncate model to specified range after build
- `--mut` : Single mutation (e.g., K76E)
- `--list` : File with one mutation per line
- `--outdir` / `--outdir-mut` : Output directories
- `--seed` : Set Modeller random seed
- `--logfile` : Custom log file path
- `--verbose` : Verbose Modeller logs

See `python3 modeller6.py --help` for all options.

---

## Directory Structure

```
modeller/
├── install_modeller.sh
├── modeller6.py
├── readme_modeller.md
├── obsolete2/
└── ...
```

---

## Notes
- Edit `install_modeller.sh` to use your own valid Modeller license key.
- All output and logs are saved in the folder named after your PDB ID.
- This protocol is designed for reproducibility and ease of use in Colab or Linux.

---
For questions or issues, please contact the repository maintainer.
