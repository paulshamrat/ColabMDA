# One-protein workflow

Use a prepared protein PDB with reviewed residue identities, protonation and
missing-residue handling. This example does not download or select a scientific
model for you. Activate the installed ColabMDA environment with AmberTools and a
working CUDA GPU first. Run GPU simulations on the compute machine, not Divergent.

```bash
export COLABMDA_RUN_ID=example_001
bash /path/to/ColabMDA/examples/basic/run.sh /path/to/project /path/to/protein.pdb
```

The default is apo WT. Optional `COLABMDA_STATE` and `COLABMDA_VARIANT` select a
prepared system; its input must contain the cofactors required by that state.
The runner calls `cmda prep`, quick equilibration (100 ps NVT, 100 ps NPT), 1 ns
production, merge, and analysis. `COLABMDA_EXAMPLE_NS` changes production duration;
use a duration divisible by 0.1 ns for complete 100 ps chunks. Quick equilibration
is a software check, not a production protocol recommendation.

An optional final argument selects `prepare`, `equilibrate`, `simulate`, or
`analyze` so you can inspect each stage. For example, run a 100 ps smoke test:

```bash
COLABMDA_EXAMPLE_NS=0.1 bash /path/to/ColabMDA/examples/basic/run.sh /path/to/project /path/to/protein.pdb simulate
```

Inputs stay outside the tool repository. Outputs belong to
`PROJECT_ROOT/data/02_pilot/RUN_ID/STATE/VARIANT/`:

- `01_prep/`: coordinates, Amber topology and preparation metadata.
- `02_eq/`: equilibration states, QC JSON and `figures/` PNGs.
- `03_md/r1/`: completed chunks, logs, checkpoint and merged trajectories.
- `04_check/`: RMSD, radius of gyration and RMSF CSVs, plus `figures/` PNGs.

Inspect thermodynamic diagnostics, merged frame times, numerical tables and plots.
A successful command alone does not establish convergence. Project workflows own
extra validation and their completion markers; see `test-colabmda` for the private
four-system integration campaign. Resume instructions are in [restart](../restart/README.md).
