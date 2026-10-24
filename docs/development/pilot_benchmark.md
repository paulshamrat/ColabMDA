# GPU integration pilot

The separate `research/projects/test-colabmda` project owns the input registry,
numbered execution scripts and run reports. ColabMDA remains installed software
under `research/tools/ColabMDA`; do not copy its implementation into a project.

## Fast default

The default software check uses Trp-cage (PDB 1L2Y, first NMR model, author chain A,
20 residues) with native residue numbering. The project retrieves mmCIF, checks
the full deposited sequence and coordinate presence, and records source/input
hashes. It invokes the installed tool's `examples/basic/run.sh` through
`scripts/02_pilot/05_example.sh`. Run 100 ps NVT, 100 ps NPT and 1 ns production,
including a 100 ps smoke test and actual SIGINT/resume of the second chunk.
Require ten complete chunks, 100 frames, continuous logged times, finite
measurements and saved plots before PILOT_DONE. This one-system apo test covers
core workflow and recovery; it does not validate mutant/cofactor workflows.

The Trp-cage deposited sequence is the reference for this designed miniprotein;
a canonical UniProt mutation annotation is not applicable. Raw coordinates remain
in the project's ignored data tree, not in tool source or examples.

## Extended CDKL5 scope

Test WT and L119R in apo and ATP–Mg states. Reviewed inputs contain 303 protein
residues; holo models contain one ATP and one Mg. Preparation records Amber
ff19SB/OPC and nucleotide parameters in `prep.json`.

For each system, use quick equilibration with 100 ps NVT and 100 ps NPT, followed
by 1 ns production at 2 fs. Save frames every 10 ps and checkpoints after each
complete 100 ps chunk. These settings test software behavior; they do not establish
biological convergence or validate the full staged equilibration protocol.

## Execution

Run standard tests on Divergent. Run protein simulations only on the Colab GPU
through the Colab CLI. Deploy committed tool and project snapshots serially,
mount Drive, and select a fresh `COLABMDA_RUN_ID`. Follow the project's numbered
workflow README rather than duplicating its deployment or input-selection logic.

1. Prepare all four systems and verify a real CUDA context.
2. Equilibrate WT apo and complete a 100 ps production smoke test.
3. Equilibrate the remaining systems.
4. Interrupt WT apo during its second production chunk with SIGINT; resume to
   1 ns and require the first committed chunk's SHA-256 to remain unchanged.
5. Complete the other three 1 ns runs, merge trajectories and analyze them.

The project wrappers invoke the public `cmda prep`, `sim`, `merge` and `analyze`
commands. They write provenance and completion markers under
`data/02_pilot/<run>/` and execution logs under `logs/02_pilot/<run>/`. Heavy
trajectories and checkpoints stay on Drive. Measure actual wall time on the
allocated GPU; no fixed speed is assumed.

## Acceptance and recovery

Require ten production chunks and 100 merged frames per system. Check finite
coordinates, finite thermodynamic logs, times from 10 to 1000 ps without gaps or
overlap, correct ATP/Mg membership, finite RMSD/Rg/RMSF tables, and nonempty PNG
figures. Inspect equilibration diagnostics separately from software completion.
`PILOT_DONE` is written only after all four systems pass downstream checks.

Resume the same run ID with the recorded input and tool revisions. New settings
or inputs require a fresh ID. Binary checkpoints are not assumed portable across
GPU/software changes. Retrieve only compact reports and figures, verify their
hashes, review the plots, and stop the dedicated Colab session after completion.

The baseline environment is Python 3.12, CUDA 13.0 and OpenMM >=8.5. Record exact
installed versions and GPU information in each run's provenance.
