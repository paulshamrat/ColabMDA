#!/usr/bin/env bash
# 00_installation.sh — Conda+Mamba installer + OpenMM stack + MDAnalysis + Biopython + PDBFixer

set -e

# 1) Install Miniforge (Conda) into $HOME/miniforge3
wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh \
     -O /tmp/miniforge.sh
bash /tmp/miniforge.sh -b -p "$HOME/miniforge3"

# 2) Initialize Conda in this shell
export PATH="$HOME/miniforge3/bin:$PATH"
source "$HOME/miniforge3/etc/profile.d/conda.sh"

# 3) Ensure 'mamba' is available in base
conda install -y -n base -c conda-forge mamba

# 4) Install all required packages into base
mamba install -y -c conda-forge \
    cudatoolkit=11.8 \
    openmm openmmtools pdbfixer \
    mdanalysis mdtraj matplotlib numpy biopython

# 5) Fallback pip install of PDBFixer (if needed)
python3 - << 'EOF'
try:
    import pdbfixer
except ImportError:
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "pdbfixer"])
EOF

# 6) Smoke-test installations
python3 - << 'EOF'
from openmm import Platform
print("OpenMM platforms:", [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())])
import MDAnalysis, mdtraj, Bio, pdbfixer
print("MDAnalysis:", MDAnalysis.__version__)
print("MDTraj:", mdtraj.__version__)
print("Biopython:", Bio.__version__)
print("PDBFixer:", pdbfixer.__version__ if hasattr(pdbfixer, "__version__") else "import OK")
EOF

echo "✅ Installation complete: OpenMM, MDAnalysis, MDTraj, Biopython, PDBFixer, and dependencies are ready."
