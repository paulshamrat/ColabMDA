# Changelog

All notable changes to **ColabMDA** will be documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased development

- Flattened OpenMM operation modules and consolidated checkpointed production.
- Added the complete Trp-cage example and safe `cmda stage --example trpcage`
  installation, including example resources in wheels and source snapshots.
- Consolidated installation into scripts/setup.sh and environment.yml, retaining
  interactive MODELLER setup and hidden key entry.
- Recorded example/input/runtime fingerprints, preserved author RMSF numbering
  with an identity gate, and separated numerical processing from figure modules.
- Consolidated the working project as colabmda-example while preserving historical
  inputs/results and verified migration records.
- Centralized user instructions in README and corrected outdated developer and
  compatibility references. These development changes are not a public release.

## [v0.2.0] - 2026-10-14 (Modular Architecture & Resilience)

### 📊 Release Summary (v0.1.0 vs. v0.2.0)

| Capability | v0.1.0 (Initial Beta) | v0.2.0 (Current Release) | Impact / Rationale |
| :--- | :--- | :--- | :--- |
| **Simulation Pipeline** | Monolithic single command | **5 Modular Stages** (`em`, `nvt`, `npt`, `qc`, `md`) | Eliminates wasted GPU hours if a run interrupts |
| **Replica Acceleration** | Re-equilibrates every replica | **Shared Hub (`eq/`)** forks `r1, r2, r3` with random seeds | **3× faster** multi-replica campaign setup |
| **Protonation Assignment** | Basic PDBFixer defaults | **PROPKA 3D pKa Titration (pH 7.4)** | Accurate AMBER Histidine states (`HIE`/`HID`/`HIP`) |
| **Solvent & Force Field** | TIP3P + ff14SB only | **OPC (4-point water) + AMBER `ff19SB`** + GAFF-2 | State-of-the-art kinase & protein accuracy |
| **Stability Validation** | Manual log inspection | **Automated Headless QC (`openmm qc`)** | Headless Matplotlib validation of density & temp |
| **Mutagenesis Engine** | Numbering reset to 1 | **Clinical / Biological Numbering (`--uniprot-numbering`)** | Preserves patient mutation indexing (e.g. L119R) |
| **Data Directory Tiers** | Flat `simulations/`, `structures/` | **Standardized Tiers**: `data/{structure,simulation,analysis}` | Clean project scaling and backward compatibility |
| **Automated Testing** | 0 tests | **Pytest Suite (`12/12 passing`)** | CI/CD regression protection |

---

### Added
- **Single-Keyword Modular Execution**: Direct CLI commands for `em`, `nvt`, `npt`, `qc`, and `md`.
- **Automated Equilibration QC**: Headless Matplotlib generation of temperature, density, and box volume convergence plots.
- **PROPKA & PDB2PQR Support**: Biological pKa titration at physiological pH 7.4 with audit logging (`protonation_summary.json`).
- **Template Re-alignment**: Modeller automatically aligns generated models back onto crystal template coordinates.
- **Comparative Multi-Series Tool**: `colabmda openmm compare` generates triplicate-averaged figures with $\text{mean} \pm \text{SD}$ bands.
- **PyMOL Live Visualization**: `colabmda openmm view` opens PyMOL centered directly on the mutation site.
- **Publication Grid Rendering**: `colabmda openmm snapshots` generates multi-panel comparative structural grids.

### Changed
- Renamed internal package `openmm_pw` to `openmm` to support protein, ligand, and multi-component systems.
- Upgraded default water model to `OPC` with automated virtual site placement (`addExtraParticles`).
- Refactored directory priority to look for `data/simulation/` first with smart fallbacks to legacy names.
- Expanded Sphinx documentation with 7 comprehensive standalone guides.

### Fixed
- Fixed AttributeError in `colabmda openmm stage` by safely resolving optional root and drive arguments.
- Fixed non-headless matplotlib crashes in Google Colab and remote SSH environments by enforcing `matplotlib.use('Agg')`.
- Fixed Modeller output filenames to cleanly match target mutation names.

---

## [v0.1.0] - 2026-06-18 (Initial Beta)

### Added
- **Modular Pipeline**: Initial CLI for EM, NVT, NPT, and Production MD.
- **Resume-Safe Engine**: Integrated checkpointing logic for fail-safe simulations on Google Colab.
- **Modeling**: Automated Wild-Type building and mutation support via Modeller.
- **Analysis**: Basic trajectory merging and comparative RMSD/Rg/RMSF analysis tools.
- **CI/CD Standards**: Added GitHub Actions workflows, Black formatting, and Ruff linting.
