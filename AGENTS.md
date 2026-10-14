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
│   ├── cli.py                     # Central CLI entrypoint (lazy-loaded modules)
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
│   └── modeller/                  # MODELLER homology and mutation engine
│       └── engine.py              # Biological numbering and mutation builder
├── tests/                         # Unit and regression test suite
├── docs/                          # Sphinx / Markdown documentation
├── scripts/                       # Colab bootstrap & remote execution scripts
└── AGENTS.md                      # This developer rules file
```

---

## 4. Naming & Execution Standards

- **Single-Keyword Stages:** Simulation stages follow clean single-keyword identifiers: `em`, `nvt`, `npt`, `qc`, `md`.
- **Equilibration Directory:** The shared equilibration baseline is strictly named `equil/` (e.g. `sim/<system_id>/equil/`). Production replicas branch off as `r1/`, `r2/`, etc.
- **Headless Plotting:** All plotting code must explicitly configure `matplotlib.use('Agg')` before importing `pyplot` to prevent crashes in headless Colab or remote terminal environments.
- **Decoupled Compute & Visualization:** Trajectory metric calculations (RMSD, Rg, RMSF, PCA) write raw numerical CSVs before rendering figures.
- **Figure Standards:** Output figures follow publication dimensions (double-column 7.0-inch width, 300 DPI, Okabe-Ito / high-contrast accessible color palettes).
- **wwPDB Standards:** Ingestion pipelines must support both PDBx/mmCIF (`.cif` / `.cif.gz`) and legacy `.pdb`.

---

## 5. Pre-Commit Quality Gates

Before committing any changes to `dev1`:
1. Run syntax verification:
   ```bash
   python -m py_compile src/colabmda/**/*.py
   ```
2. Run code formatting & linting checks:
   ```bash
   black --check .
   ruff check .
   ```
   *(If formatting fails, auto-format using `black .`)*
3. Run the test suite:
   ```bash
   pytest tests/
   ```

