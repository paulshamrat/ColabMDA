# Simulation protocol

Installation and the complete Trp-cage workflow are in the main repository README.
This page describes the scientific settings and recovery behavior implemented by
the current checkout. Use [the reference](reference.md) for command scope.

## Preparation routes

The high-level `cmda prep` path assembles a reviewed input model and builds an
Amber topology with LEaP. Its defaults are ff19SB protein, OPC water, 14 Å solvent
padding, neutralizing ions and nominal 0.15 M NaCl. Salt-pair counts are estimated
from the padded solute dimensions. Inspect `prep.json` and the actual topology
rather than treating the requested concentration as an exact measured value.

The lower-level native OpenMM setup can instead construct its system from XML
force fields. It supports ff19SB or ff14SB and OPC, TIP3P, TIP3P-FB or TIP4P-Ew when
those files are installed. Its auto padding is 1.0 nm for a solute up to 8 nm
across, 1.2 nm up to 12 nm, and 1.4 nm beyond that. Padding below 1.0 nm is rejected.
That native route neutralizes the system without adding extra bulk salt by default.
Its four-site OPC construction uses TIP4P-Ew starting geometry with OPC parameters.
These defaults do not override a supplied Amber topology.

`cmda openmm prep` supports PDBFixer/PDB2PQR preparation with PROPKA at the chosen
pH and records protonation decisions. Review explicit AMBER residue forms such as
HID/HIE/HIP and ASH/GLH. Do not assume high-level assembly performs the same
protonation calculation or that a missing external program completed it.

Generic ligand construction uses the selected GAFF template generator and charge
assignment. Specialized ATP/Mg parameters are a different chemical model. To retain
an existing parameterization, supply its matching Amber topology and coordinates;
inspect the recorded parameter source, cofactor identities and ion parameters.
Never describe generic ligand parameters as equivalent to specialized nucleotide
parameters merely because the coordinates look similar.

## Dynamics settings

The shared OpenMM setup uses PME with a 1.0 nm direct-space cutoff, switching
from 0.8 nm, hydrogen-bond constraints and rigid water. Center-of-mass motion is
removed every 100 steps. The thermostat is LangevinMiddle with 1/ps friction,
usually at 300 K. NPT uses an isotropic Monte Carlo barostat at 1 bar, attempting
a volume move every 25 steps. CUDA/OpenCL use mixed precision by default.

The default platform search is CUDA, OpenCL, then CPU. Set
`COLABMDA_PLATFORM=CUDA` for an explicit platform and `COLABMDA_REQUIRE_GPU=1` to
reject CPU fallback. `COLABMDA_PRECISION` accepts mixed, single or double.
The generated Trp-cage scripts require CUDA; Divergent runs CPU tests only.

## Equilibration profiles

Select `cmda sim --protocol staged` or `--protocol quick` explicitly. The lower-level
stage CLI calls the same staged profile `default`/`varmdyn`. Positional restraint
strengths below are in kcal mol⁻¹ Å⁻² and apply to solute heavy atoms; water and
bulk monovalent ions remain unrestrained.

| Stage | Operation | Conditions | Restraint |
| :--- | :--- | :--- | :--- |
| EM | Initial minimization | No physical time | None |
| 01–03 | Minimization ladder | Constant box | 100, 50, 25 |
| 04–09 | Six heating intervals | 10 ps each; 50–300 K; NVT; 1 fs | 25 |
| 10–14 | Minimization ladder | Constant box | 25, 10, 5, 2, 1 |
| 15–18 | NVT release ladder | 100 ps each; 300 K; 2 fs | 10, 5, 2, 1 |
| 19–22 | NPT release ladder | 100 ps each; 300 K; 1 bar; 2 fs | 10, 5, 2, 1 |
| 23 | Unrestrained NPT | 100 ps; 300 K; 1 bar; 2 fs | None |
| 24 | Final unrestrained NPT | 1000 ps; 300 K; 1 bar; 2 fs | None |

Completed numbered stages write portable states and JSON summaries under the
equilibration directory. Recovery finds the highest completed stage rather than
repeating earlier stages. These compact records are not equilibration trajectories.

Quick equilibration uses a short minimization, one unrestrained NVT interval and
one NPT interval. In the Trp-cage example these are 100 ps each at 2 fs. The combined
legacy `openmm run` adapter also delegates to quick equilibration while retaining
its explicit ff14SB/TIP3P defaults. Do not equate that adapter with the full staged
schedule or assume it reproduces historical monolithic trajectories bit for bit.

## QC and completion

Numerical QC writes `equilibration_qc.json`; downstream plotting writes PNG under
`figures/`. Inspect temperature, density and density drift. No QC PDF is produced.

The high-level `cmda sim` route reports QC without enforcing a production gate.
A failed QC result therefore needs an explicit scientific decision before production;
a successful process exit or an `EQ_DONE` marker does not certify equilibration.
The separate `cmda openmm qc` command exits unsuccessfully on failed QC unless
`--warn-only` is supplied. Its aliases and options are listed in the CLI help.

The example's final `RUN_DONE` means its software checks passed. Physical QC is
reported separately in `summary.json`. A short example cannot establish convergence,
validate a protein's biology, or certify a full mutant/cofactor campaign.

## States, checkpoints and replicas

Portable `*.state.xml` files transfer positions, velocities and periodic boxes
between ensembles or into independent replicas. Binary `*.chk` files also contain
Context/random-generator internals and are for continuation within the same
compatible System, hardware and software environment.

A new production replica starts from the NPT physical state, assigns fresh 300 K
velocities and stochastic streams, and resets production time to zero. Independent
replicas can share the equilibration endpoint; they do not branch from a sibling's
binary checkpoint. Inspect saved states and seed records when reproducing a run.

`--total-ns` is the absolute production target per replica, excluding equilibration.
For example, rerunning with 1 ns extends a compatible 0.1 ns run to 1 ns total.
It does not add another 1 ns. At 10 ps sampling, 1 ns produces 100 frames, with
logged times from 10 to 1000 ps.

The production implementation publishes a complete trajectory/log chunk before
advancing the checkpoint. Interruption replays the unfinished chunk from the last
committed checkpoint. Verify unchanged committed chunk hashes and merged times
without gaps or overlaps; filenames and process exit alone are insufficient.
Do not assume binary checkpoints are portable across GPU/software changes.

## Provenance and analysis

Retain the input archive/selection hashes, author numbering, preparation metadata,
force-field topology, parameterization record, runtime packages, GPU identity,
executed script fingerprints and tool/project revisions. Preserve historical
provenance verbatim during folder migration and record the old-to-new path mapping.

RMSD is measured against the first trajectory frame with backbone fitting. The
reported Rg uses protein heavy-atom coordinates and an unweighted geometric center.
RMSF uses C-alpha coordinates aligned to the first frame. Native author labels are
mapped through a matching reference structure, including insertion codes and chain
disambiguation. `analysis.json` records the analysis implementation and reference.
Numerical tables and downstream figure rendering remain separate operations.

## Further reading

- [OpenMM application guide](https://docs.openmm.org/latest/userguide/application/02_running_sims.html)
- [OpenMM theory guide](https://docs.openmm.org/latest/userguide/theory.html)
- [OpenMM checkpoint/state API](https://docs.openmm.org/latest/api-python/generated/openmm.app.simulation.Simulation.html)
- [OpenMMForceFields](https://github.com/openmm/openmmforcefields)
