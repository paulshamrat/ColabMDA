# Simulation Protocol, Physics, and Restart Semantics

This page describes what ColabMDA actually does at every stage. The default
`--equil-protocol default` uses a full staged AMBER-style equilibration schedule with
OpenMM implementations of standard minimization, NVT, NPT, thermostat, barostat, PME,
and restraint methods. It is not a claim that two different MD engines generate
bitwise-identical trajectories.

## Force field, protonation, and system construction

The native OpenMM XML path uses AMBER **ff19SB** for protein and **OPC** water. Generic
small molecules use **GAFF 2.11**. If the input ligand has no partial charges, OpenFF
assigns AM1-BCC charges; existing nonzero input charges are preserved. This generic
route is not chemically equivalent to specialized ligand parameter sets.

### Protonation and pKa Titration (PDB2PQR & PROPKA)

ColabMDA uses **PDB2PQR** with **PROPKA** titration during structure preparation (`colabmda openmm prep --ph 7.4`).
Rather than relying on generic hydrogen placement, PROPKA evaluates local 3D electrostatic microenvironments, desolvation effects, and hydrogen-bonding networks to calculate pKa values and assign explicit AMBER titration forms:

* **Histidines:** Assigned as `HID` (Nδ1-protonated), `HIE` (Nε2-protonated), or `HIP` (doubly protonated, +1 charge).
* **Glutamates & Aspartates:** Assigned as `GLH` / `ASH` (neutral protonated) when buried or titrating above target pH, vs `GLU` / `ASP` (-1 charge).
* **Lysines & Cysteines:** Assigned as `LYN` (neutral lysine) and `CYX` (disulfide-bonded cysteines).

An audit log (`protonation_summary.json`) is recorded in the preparation folder for full reproducibility. Downstream OpenMM simulation stages (`em`, `nvt`, `npt`, `md`) natively map these explicit AMBER titration forms directly to their corresponding force field parameters.

ColabMDA examples and recommended calculations use `--protein-ff ff19SB
--water-model opc`. Advanced users can choose other installed OpenMM XML combinations:

| Option | Choices | Notes |
|---|---|---|
| `--protein-ff` | `ff19SB`, `ff14SB` | `ff19SB` requires an OpenMM install with Amber19 XML files. |
| `--water-model` | `opc`, `tip3p`, `tip3pfb`, `tip4pew` | OPC is the default; TIP3P/TIP3P-FB/TIP4P-Ew are available for compatibility or comparison. |
| `--small-molecule-ff` | `gaff-2.11` by default; newer GAFF versions such as `gaff-2.2.20` if installed | Used only for `--ligand` native OpenMM setup. |

Current OpenMM documentation lists Amber19 XML files including `amber19-all.xml`,
`amber19/protein.ff19SB.xml`, and `amber19/opc.xml`. The `openmmforcefields`
project provides additional AMBER/CHARMM/OpenFF/Espaloma support and GAFF small
molecule template generation.

For specialized ATP-Mg systems, use the Hu 2024 ATP PREPI/FRCMOD parameters with
ff19SB, OPC, and OPC-compatible ions when that is the chemistry you intend to
reproduce. To preserve those exact parameters, build the solvated AMBER topology with
LEaP and give ColabMDA both files:

```bash
colabmda openmm run --name kras_atpmg --replica r1 --amber-prmtop system.prmtop --amber-inpcrd system.inpcrd --total-ns 100
```

Do not substitute generic GAFF2 for Hu 2024 ATP-Mg and describe the result as the same
chemical model.

Both native paths use PME, a 1.0 nm direct-space cutoff, force switching from 0.8 to
1.0 nm, rigid water, constraints on bonds involving hydrogen, and center-of-mass
motion removal every 100 steps. The system is neutralized with counterions. No extra
bulk salt is added by default.

OpenMM `Modeller` has no named OPC box builder. Following its documented mechanism for
four-site water models, ColabMDA places waters from the built-in TIP4P-Ew geometry,
assigns **OPC parameters** from `amber19/opc.xml`, and minimizes the resulting system.
The simulated model is OPC; TIP4P-Ew supplies only the initial four-site box geometry.

## Solvent box selection

`--padding-nm auto` measures the unsolvated solute during EM. It chooses 1.0 nm for a
compact solute up to 8 nm across, 1.2 nm for an extended 8-12 nm solute, and 1.4 nm
above 12 nm. It prints the solute spans, selected padding, and estimated box dimensions
before adding water. KRAS/4BGQ normally selects 1.0 nm, which is appropriate for the
Colab tests and substantially reduces solvent/GPU cost.

Use `--padding-nm 1.4` for a conservative 14 A LEaP-style padding. Values
below 1.0 nm are rejected because the padding would be smaller than the nonbonded
cutoff. Padding is a boundary margin, so a smaller protein does not by itself justify
going below the cutoff.

## Default staged equilibration profile

All positional restraint strengths below are in kcal mol-1 A-2 and act on solute heavy
atoms. Water and bulk Na+/Cl-/K+ ions are unrestrained. Langevin friction is 1/ps.

| Stage | Operation | Ensemble and conditions | Restraint | Purpose |
|---|---|---|---:|---|
| EM | Initial minimization | Periodic system; no physical time | none | Remove severe clashes after construction |
| 01-03 | Minimization ladder | Constant box | 100, 50, 25 | Relax solvent while holding the solute near its input structure |
| 04-09 | Six heating stages, 10 ps each | NVT; 50, 100, 150, 200, 250, 300 K; 1 fs | 25 | Heat gradually rather than shocking the system |
| 10-14 | Minimization ladder | Constant box | 25, 10, 5, 2, 1 | Relax after heating while releasing restraints |
| 15-18 | Four stages, 100 ps each | NVT; 300 K; 2 fs | 10, 5, 2, 1 | Equilibrate temperature and progressively release the solute |
| 19-22 | Four stages, 100 ps each | NPT; 300 K; 1 bar; 2 fs | 10, 5, 2, 1 | Relax density and box while releasing the solute |
| 23 | 100 ps | NPT; 300 K; 1 bar; 2 fs | none | First fully unrestrained equilibration |
| 24 | 1000 ps | NPT; 300 K; 1 bar; 2 fs | none | Final unrestrained equilibration and production branch point |

Each numbered stage writes `equilibration/stage_NN.state.xml` plus a small
`stage_NN.json` summary containing its description, final step/time, and potential
energy. These files provide a compact audit trail for the numbered equilibration
stages. ColabMDA intentionally does not write a large equilibration trajectory for
every stage by default; production trajectories are unchanged.

After stage 24, ColabMDA evaluates the final temperature, density range, and density
drift and writes `equilibration_qc.json`, PNG, and PDF. A failed QC gate blocks
production under the full profile. Inspect the system rather than bypassing the gate;
`check-equil --warn-only` is available only for deliberate diagnostics.

The OpenMM Monte Carlo barostat attempts an isotropic volume move every 25 steps. The
thermostat is `LangevinMiddleIntegrator`. These are OpenMM's standard implementations,
documented in the [OpenMM simulation guide](https://docs.openmm.org/latest/userguide/application/02_running_sims.html).

## Shared equilibration and multi-replica branching

To minimize redundant compute, ColabMDA supports shared-equilibration branching for multi-replica simulations:

* **Single Equilibration Endpoint:** Running `colabmda openmm run --equil-only` executes Energy Minimization, NVT heating, and NPT density equilibration once, outputting `system.xml`, `solvated.pdb`, and `npt.state.xml`.
* **Automatic Replica Auto-Discovery:** Subsequent replicas (`r1`, `r2`, `r3`, ...) auto-discover `npt.state.xml` from `equil/`, `r1/`, or the parent simulation directory (`../`).
* **Stochastic Velocity Resampling:** Replicas do not load Context checkpoints from sibling runs. Instead, each replica initializes a fresh Context from `npt.state.xml` and assigns velocities at target temperature (300 K) using a **unique per-replica random seed** derived via `derive_replica_seed(seed, replica_name)`.

## GPU selection and precision

ColabMDA selects OpenMM platforms in `CUDA -> OpenCL -> CPU` order and prints the
chosen platform whenever it creates a Context. CUDA/OpenCL use mixed precision by
default. For a production Colab run, set `COLABMDA_REQUIRE_GPU=1` to abort instead of
silently falling back to CPU. Advanced overrides are `COLABMDA_PLATFORM=CUDA` and
`COLABMDA_PRECISION=mixed|single|double`. CPU is appropriate for local construction
tests, not for a long production trajectory.

### Quick test profile

`--equil-protocol quick --equil-time 100` performs one unrestrained 100 ps NVT stage
and one 100 ps NPT stage. It exists for installation, GPU, and restart smoke tests. It
is not the full production equilibration protocol and should not be represented as one.

## State files versus checkpoints

Files named `*.state.xml` carry positions, velocities, and periodic box vectors between
EM, NVT, NPT, and new replicas. This is deliberate: those stages use different
OpenMM `System` objects as the barostat and restraint forces change. Binary `*.chk`
files contain complete Context and random-generator internals and are retained only for
exact continuation within an unchanged stage/System.

Production starts from `npt.state.xml`, creates a fresh Context, draws new 300 K
velocities, assigns independent Langevin and barostat streams, and resets its clock to
zero. `prod.chk` then resumes only that same replica. See OpenMM's
[checkpoint/state API](https://docs.openmm.org/latest/api-python/generated/openmm.app.simulation.Simulation.html).

## Production time and replicas

`--total-ns` is the absolute production target **per replica** and excludes all EM and
equilibration time. Reissuing `--total-ns 0.5` stops at 500 ps; it does not add 500 ps.
`--checkpoint-ps 100` divides output/recovery into 100 ps chunks.

Independent replicas share `system.xml`, `solvated.pdb`, and `npt.state.xml`, but never
another replica's full checkpoint. Each receives fresh velocities and stochastic seeds.
Later chunks of the same replica resume from its own `prod.chk`.

## Reproducibility record

EM writes `parameterization.json` with the force-field source, protein/water/ligand
models, charge source, pressure, cutoff, switching distance, and chosen solvent padding.
Archive it with `system.xml`, `solvated.pdb`, and the command line used for the run.

## Primary references

- [OpenMM application guide](https://docs.openmm.org/latest/userguide/application/02_running_sims.html)
- [OpenMM theory guide](https://docs.openmm.org/latest/userguide/theory.html)
- [OpenMMForceFields GAFF generator](https://github.com/openmm/openmmforcefields)
- Eastman P, et al. OpenMM 8. *J Phys Chem B*. 2024.
  [doi:10.1021/acs.jpcb.3c06662](https://doi.org/10.1021/acs.jpcb.3c06662)
