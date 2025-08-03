# Modeller Workflow in ColabMDA/modeller

This folder provides a streamlined workflow for running Modeller (homology modeling) in a Colab or Linux environment. The workflow includes scripts for automated installation, environment setup, license configuration, and running example jobs.

## Workflow Steps

1. **Make Scripts Executable**
   ```bash
   chmod +x install_modeller.sh
   chmod +x run_modeller2.py
   ```

2. **Install Modeller**
   Run the installation script in your Colab or terminal:
   ```bash
   ./install_modeller.sh
   ```
   This sets up a dedicated conda environment (`modeller_env`) and installs Modeller.

3. **Configure License and Activate Environment**
   After installation, run the following commands to activate the environment and set the license key:
   ```bash
   source ~/miniforge/etc/profile.d/conda.sh
   conda activate modeller_env
   # Replace YOUR_LICENSE_KEY with your actual Modeller license key
   LICENSE_KEY="YOUR_LICENSE_KEY"
   CONFIG="$HOME/miniforge/envs/modeller_env/lib/modeller-10.7/modlib/modeller/config.py"
   sed -i "s/^license *=.*/license = '${LICENSE_KEY}'/" "$CONFIG"
   grep "^license" "$CONFIG"
   python -c "import modeller; print('Modeller OK, version', modeller.__version__)"
   ```
   This confirms that Modeller is installed and licensed correctly. You must use your own valid Modeller license key.

4. **Install Biopython and Run Modeller Example**
   Before running the Modeller script, install Biopython (required for many modeling tasks):
   ```bash
   pip install biopython
   ```
   Then use the provided script to run a modeling job:
   ```bash
   python3 run_modeller2.py
   ```
   Edit or extend this script for your own homology modeling tasks.

## Files
- `install_modeller.sh`: Automated installer for Modeller and environment setup.
- `run_modeller2.py`: Example script to run Modeller jobs.

## Notes
- You must use your own valid Modeller license key (do not use the example key above).
- This workflow is designed for reproducibility and ease of use in Colab or local Linux environments.

---
For questions or issues, please contact the repository maintainer.
