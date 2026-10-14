# ColabMDA Developer & Agent Rules

## 1. Project Overview & Scope

ColabMDA is a resume-safe, modular Molecular Dynamics (MD) and structural modeling engine tailored for Google Colab GPUs, local workstations, and HPC environments.
It combines **MODELLER** (homology modeling & mutagenesis with biological residue numbering) and **OpenMM** (chunked, fail-safe MD simulations with automated Google Drive checkpointing).

---

## 2. Branching & Release Architecture

- **Active Development Branch (`dev1`):** All active development, modular refinements, tests, and bug fixes occur on `dev1`.
- **Docking Module (`docking`):** Virtual screening and AutoDock Vina workflows are isolated in the `docking` branch for future standalone updates.
- **Public Release (`main`):** Reserved strictly for clean, verified, production-ready releases on GitHub. Never commit directly to `main` or merge unvetted changes.

---

## 3. Directory Layout & Architecture

```text
tools/ColabMDA/
├── src/colabmda/
│   ├── cli.py                     # Central CLI entrypoint (cmda and colabmda)
│   ├── core/                      # System taxonomy, ion registry, project paths
│   │   ├── taxonomy.py            # SystemSpec definition & state-tag parser
│   │   ├── ions.py                # Divalent ion parameters & registry
│   │   └── project.py             # ProjectPaths layout resolver
│   ├── amber/                     # AmberTools parameterization & cpptraj engine
│   │   ├── params/                # Hu 2024 quantum nucleotide parameters
│   │   ├── leap.py                # Fast tleap solvator with addIonsRand
│   │   └── cpptraj.py             # Dual-stream trajectory autoimaging & stripper
│   ├── modeller/                  # MODELLER homology and substrate engine
│   │   ├── engine.py              # Biological numbering and mutation builder
│   │   └── substrate.py           # Mathematical PTR builder & peptide docking
│   ├── openmm/                    # Core OpenMM simulation engine
│   │   ├── commands.py            # High-level subcommands (prep, stage, run, merge, analysis)
│   │   ├── modular/               # Sequentially numbered simulation stages
│   │   │   ├── em.py              # Energy Minimization & forcefield selection
│   │   │   ├── nvt.py             # NVT thermalization
│   │   │   ├── npt.py             # NPT density equilibration
│   │   │   ├── qc.py              # Quality control & headless stability verification
│   │   │   ├── md.py              # Production MD & chunked checkpointing
│   │   │   └── utils.py           # Shared OpenMM helpers
│   │   └── engine/                # Trajectory merge, analysis & PDBFixer cleaning
│   └── remote/                    # Google Colab VM remote session manager
│       └── colab.py               # colab CLI wrapper & upload orchestrator
├── tests/                         # Unit and regression test suite (40+ passing tests)
├── docs/                          # Sphinx / Markdown documentation
├── scripts/                       # Exactly one standalone bootstrap script: setup.sh
└── AGENTS.md                      # Developer and agent rules
```

---

## 4. Architectural Rules & Invariants

1. **Zero Data in `src/`:**
   - Never commit raw structural files (`*.pdb`, `*.cif`, `*.sdf`, `*.mol2`), topologies (`*.prmtop`, `*.inpcrd`), or trajectories (`*.dcd`, `*.nc`) into `src/colabmda/`.
   - All structural inputs belong strictly in the user's project workspace:
     - Receptor structures: `<project>/structure/protein/`
     - Template structures: `<project>/structure/template/`
   - Non-standard residues like phosphotyrosine (`pY171` -> `PTR171`) are constructed mathematically via pure NumPy vectors with zero external template files.

2. **Decoupled Scripts Architecture:**
   - Never place bash or SLURM scripts inside simulation data folders (e.g., no `simulation/atp_2mg_eb2/WT/run.sh`).
   - All execution scripts live in `<project>/scripts/`:
     - Standalone local/remote runner: `<project>/scripts/run.sh`
     - HPC Palmetto scheduler: `<project>/scripts/submit.slurm`

3. **Single Bootstrap Script:**
   - Exactly one bootstrap script exists in the repository: `scripts/setup.sh`.
   - Never introduce loose or duplicate shell bootstrap files.

4. **Single-Keyword Unified CLI:**
   - The CLI is accessible via both `colabmda` and the shortcut `cmda`.
   - Core subcommands:
     - `cmda stage`: Prepares standard folder layout and decoupled scripts.
     - `cmda prep`: Mathematical assembly and Amber `tleap` parameterization.
     - `cmda sim`: Executes simulation stages (`em`, `nvt`, `npt`, `qc`, `md`).
     - `cmda merge`: Autoimages PBC and strips bulk solvent via `cpptraj`.
     - `cmda analyze`: Extracts RMSD, Rg, RMSF, PCA, and contact maps.
     - `cmda compare`: Multi-replica or cross-variant comparative plotting.
     - `cmda colab`: Provisions and manages remote Google Colab GPU sessions.

5. **Post-Verification Compute Cleanup:**
   - Always terminate Colab GPU VM sessions (`colab stop -s <session>`) immediately after use to conserve user cloud compute units.
   - Clean up temporary local tarballs and scratch files after test completion.

---

## 5. Pre-Commit Quality Gates

Before committing any changes to `dev1`:
1. Run syntax verification:
   ```bash
   python -m py_compile src/colabmda/**/*.py
   ```
2. Run code formatting & linting checks:
   ```bash
   black --check src tests
   ruff check src tests
   ```
   *(If formatting fails, auto-format using `black src tests`)*
3. Run the test suite:
   ```bash
   pytest tests/
   ```

---

## 6. Google Colab CLI & Remote Cloud Testing

The developer workstation is pre-authenticated with Google OAuth2 (`colab whoami`). When testing cloud bootstrapping or remote GPU execution:

1. **Check Authentication & Identity:**
   ```bash
   colab whoami
   ```
2. **Launch a Cloud GPU Test Instance:**
   ```bash
   colab new -s colabmda-gpu --gpu T4
   ```
3. **Execute Commands Remotely:**
   ```bash
   echo '!nvidia-smi' | colab exec -s colabmda-gpu
   ```
4. **List Active Sessions:**
   ```bash
   colab ls
   ```
5. **Terminate Sessions Promptly (Conserves Compute Units):**
   ```bash
   colab stop -s colabmda-gpu
   ```



