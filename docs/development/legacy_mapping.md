# Developer Reference: Backward-Compatibility Mapping (v0.2.0 to v0.3.0)

> 🧑‍💻 **Internal Developer Reference**  
> This document maps legacy multi-keyword commands (`colabmda openmm ...`) to their modern single-keyword `cmda` equivalents for backward compatibility maintenance, testing, and migration of legacy scripts.

| Legacy Command (v0.2.0) | Modern Command (v0.3.0) | Key Advantages |
| :--- | :--- | :--- |
| `colabmda modeller build ...` | `cmda modeller build ...` | Direct top-level access, biological numbering by default. |
| `colabmda openmm prep ...` | `cmda prep --state <tag> --variant <var>` | Combines PDBFixer, PROPKA 3 at pH 7.4, and Amber tleap solvation into one step. |
| `colabmda openmm stage ...` | `cmda stage --state <tag> --variant <var>` | Auto-creates standard `data/simulation/<state>/<var>` layout and decoupled scripts. |
| `colabmda openmm em/nvt/npt/qc` | `cmda sim --state <tag> --stage eq` | Executes full 24-stage relaxation hub with automated QC gating in one command. |
| `colabmda openmm md ...` | `cmda sim --state <tag> --stage md --replica r1` | Resume-safe checkpointed production MD. |
| `colabmda openmm merge ...` | `cmda merge --state <tag> --variant <var>` | Automated `cpptraj` PBC autoimaging and dual-stream solvent stripping. |
| `colabmda openmm analysis ...` | `cmda analyze --state <tag> --variant <var>` | Generates time-series metrics (RMSD, Rg, RMSF, PCA). |
| `colabmda openmm compare ...` | `cmda compare --state <tag> --variants WT,L119R` | Multi-replica and cross-variant comparative overlays. |
| `colabmda openmm view ...` | `colabmda openmm view --pdb-dir <dir>` | Generates automated `visualize.pml` and launches PyMOL. |
