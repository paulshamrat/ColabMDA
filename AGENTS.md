# ColabMDA development rules

## Layout and scope

- Work on private `dev1`. Public `main` is reserved for approved clean releases.
- `src/colabmda/` contains implementation: core project/example handling, MODELLER,
  AmberTools, flat OpenMM operation modules and remote Colab orchestration.
- Keep OpenMM flat: prepare, em, nvt, npt, md, qc, merge, analysis, figures,
  common, commands and the combined-command compatibility adapter run.
- `scripts/setup.sh` is the only installer. `environment.yml` is its readable
  core dependency specification. Optional MODELLER stays separate for license
  and Python compatibility. Do not recreate envs/ or duplicate bootstrap scripts.
- `examples/{trpcage,fkbp,hras}/scripts/` contains thin user workflows, not engines.
  Main README retains its complete baseline. Shared CDKL5 scientific commands live
  in docs/apo.md and docs/atpmg.md; docs/platforms/{colab,hpc}.md owns platform
  setup and execution. Keep one scientific workflow across platforms; no user YAML.
  examples/hpc/simulate.slurm is a thin scheduler adapter, not another MD engine.
  examples/hpc/analyze.slurm runs shared merge/analysis on CPUs without a GPU.
- `cmda stage --example trpcage --root ...` installs templates and records hashes
  in data/00_setup/example.json. Ship templates in wheels and source snapshots.
  Refuse symlinks and different existing files; never overwrite user settings.
- One working project: research/projects/colabmda-example on Divergent and Drive.
  Its numbered scripts are copied templates; tests exercise those same scripts.
  Historical runtime paths are data/01_structure and data/02_simulation/<run>/,
  with logs/02_simulation/<run>/. Preserve their provenance verbatim.
- New direct CDKL5 runs use data/{apo,atpmg}/structure/{raw,mapping,models} and
  data/<state>/simulation/<run>/<variant>/{01_prep,02_eq,03_md} and
  data/<state>/analysis/<run>/<variant>/{r1,r2,r3,aggregate}, with
  logs/<state>/<run>/<variant>. State, structure and simulation stay unnumbered.
  Existing numbered runs remain readable; never migrate active scientific outputs.
- CDKL5 uses O76039, 4BGQ A, residues 1–303. Derive canonical positions from
  archive reference mappings; distinguish observed coordinates from modeled gaps.
  The standalone route fetches public inputs, builds WT/L119R, and transfers
  ATP/one-Mg from public 8FP5 with homologous sequence/CA fitting and explicit
  fixed-backbone pocket sidechain minimization. Record it as a new modeled pose,
  never as the reviewed activation-study complex. Keep the explicit legacy-pose
  route available. Check identities, fit coverage, clashes and hashes. Matched
  apo/holo models share the source protein backbone and modeling method.
- Research uses full 24-stage equilibration, continuous physical time and strict
  timing/energy/thermodynamic QC before MD. Quick/reduced profiles require explicit
  --software-test. Check solvent/Amber salt naming so bulk ions stay unrestrained.
  Validate noncovalent component boundaries before LEaP to prevent protein–ATP bonds.
- Never vendor external software or put coordinates/trajectories in src/examples.
  Scientific inputs, caches, reports, checkpoints and figures stay in ignored data.
  Compiler decks run in temporary sandboxes, not scientific data folders.

## Scientific and recovery invariants

- Preserve native author numbering/insertion codes. Structure retrieval supports
  mmCIF and PDB. Check sequence identity and C-alpha coordinate presence.
- Trp-cage is a designed peptide; verify its deposited sequence, model 1 and
  author chain A. A canonical UniProt mutation annotation is not applicable.
- The historical small-protein guide uses FKBP12 apo (1D6O A; author 1–107 maps to current P62942
  2–108, WT/F37V) and H-Ras GDP/Mg (4Q21 A; P01112 1–166, WT/G12V).
  Preserve raw author numbering; derived models use verified UniProt numbering.
  Check exact sequences/C-alpha presence, stop on mismatches, and retain mapping
  and model hashes. Holo coordinates come from the same experimental chain/frame.
- Natural-protein scripts install under scripts/<stage>/<name>/ with separate
  data/00_setup/<name>.json receipts so they coexist with historical Trp-cage.
  Use one run ID per protein/variant. Changing code/protocol requires a new run.
  Confirm actual nucleotide/ion topology composition and charges, not just files.
- The default Amber route uses water-matched Li–Merz 12–6 ions. Never load a
  12–6–4 ion frcmod without implementing and verifying the separate C4 term.
  Natural-protein acceptance checks actual OPC Na/Cl/Mg LJ values. Changing ion
  parameters requires fresh preparation/equilibration and new run IDs.
- RMSF author labels require a reference identity/count gate, including insertion
  codes and chain disambiguation. Never use Amber indices as author numbers.
- Numerical calculation modules never import plotting libraries. Plot only saved
  numerical outputs downstream; use 7-inch, 300-DPI PNGs and workspace typography.
- Keep one production implementation. Combined legacy CLI retains explicit
  ff14SB/TIP3P defaults while delegating to shared equilibration/production stages.
- Advance checkpoints only after complete trajectory/log chunks are published.
  Interrupted chunks replay from the previous checkpoint without time gaps.
- Validate nonempty outputs before completion markers. RUN_DONE means software
  checks passed; physical equilibration QC remains separately visible.
- Quick equilibration is a software test, not biological convergence or full
  staged-protocol validation. Never claim mutant/cofactor coverage from Trp-cage.

## Installation, secrets and remote execution

- Preserve the optional interactive MODELLER question and hidden key prompt.
  Automated tests may use an existing key through the process environment.
  Never record keys in source, provenance, global shell startup files or logs.
- Remote MODELLER setup prompts locally and transfers a mode-600 temporary key
  file, never a key literal in remote code; remove temporary credentials afterward.
- Use cmda CLI/official Python API for operations; do not create ad-hoc temporary
  workflow scripts. Project deployment/migration belongs in its numbered sync code.
- Deploy committed dev1 snapshots through cmda colab setup --local-source, record
  SNAPSHOT.json hashes and exclude Git metadata, build copies and private material.
- Run protein GPU pilots in a Palmetto allocation or on Colab, never Divergent.
  Prefer compbio for current pilots. Keep Colab jobs foreground with streaming
  output; use CPU-only sessions for Drive cleanup when possible.
- Stop dedicated sessions after verified retrieval; heavy trajectories stay on Drive.

## Verification

Before commits, run Black (one file per process if multiprocessing hangs),
ruff format --check, ruff check, Python compilation and pytest tests/.
Run bash -n on changed shell scripts. Test project creation, conflict/symlink
protection, argument forwarding, wheel resources and hidden MODELLER prompting.
GPU acceptance includes real SIGINT/resume, unchanged committed chunk hashes,
continuous frame times, finite measurements, native numbering and reviewed plots.
Do not claim a new GPU run when only CPU/template tests were performed.

## Documentation upkeep

README.md is the complete maintained user guide. Preserve its established project
information, numbered platform setup, project architecture, full modeling/apo/holo
simulation protocol, small example, acknowledgements, changelog and citation
sections. Make focused updates; do not replace the full guide with a quickstart
or reorganize it without an explicit user request.

The root README's main workflow uses CDKL5 WT/L119R for apo and ATP/one-Mg,
matching docs/apo.md and docs/atpmg.md. Keep the existing section order and use
simple, direct, concise language. Trp-cage remains the optional quick test;
FKBP12/H-Ras templates and historical evidence remain available. Do not claim
L119R GPU validation from its CPU model/preparation checks. The default walkthrough
must run from public inputs without another project or owner-supplied files.
Use short cmda commands; put copying, hashing and checks in package operations.
COLABMDA_PROJECT_DIR removes repeated --root arguments. Keep optional Python
diagnostics separate from the main command sequence. Preserve historical study
models and pilot evidence when adding the public-input route.

Keep architecture consolidated in README section 2: tool repository, then shared
simulation project. The optional Trp-cage tutorial and layout belong in docs/examples.md; section 4
is the short 4BGQ pilot.
State guides link to the shared layout instead of repeating its tree. Keep one
compact command reference near the end (section 5), with options available through
--help. Workflow sections retain runnable steps, not duplicate command catalogs.
Place one badge row directly below the title: release, CI, license, repository,
Issues, Black and Ruff. Keep the development/style badges visible there.
Follow it with package purpose and where to start. Keep platform/dependency details
in Platform Setup, with development/support links in CONTRIBUTING.md.
Put state-guide links with the workflow and test outcomes in section 6, Validation.

Developer rules belong in CONTRIBUTING.md and AGENTS.md. Maintain ordinary
Markdown workflow guides under docs/ with nested platform setup, as explicitly
requested by the user. Do not restore ReadTheDocs/Sphinx or duplicate apo/holo
scientific procedures separately for Colab and HPC. Update the root README,
AGENTS.md and local ignored PLAN.md together when design/status changes.
Keep commands accurate, label compatibility limits, state which machine runs each
step, and check shell syntax, CLI arguments and local links after guide changes.
Preserve the timestep/unit conversions, integration-step, frame/chunk, stride and
storage calculations and per-replica time semantics when editing the guide.
Validate their examples against implementation settings and explicit arithmetic.
Distinguish completed CPU/script checks from actual GPU integration testing.

HPC access follows workspace hpc CLI rules. Inspect resources before selecting
Slurm settings. Test a single WT allocation before a panel or arrays; dynamic
manifest bounds and dependencies belong in a later verified campaign adapter.
Record scheduler/syntax checks separately from actual GPU acceptance. Do not
modify or interrupt the active docking campaign while checking ColabMDA support.

Remote multiline shell blocks must run in one Bash process, preserving cd/export and
propagating nonzero exits. Independent exec calls do not share shell state.

HPC adapters may select an existing user environment with COLABMDA_ENV; never
assume login-node CUDA plugin visibility proves allocated-GPU compatibility.
Require a real device Context probe and a separately named software smoke before
full molecular HPC validation. Enforce COLABMDA_REQUIRE_GPU=1. Never reuse quick
smoke outputs as research equilibration. Code remains a committed runtime snapshot
under scratch tools/, with the execution project separate under scratch projects/.

Use compbio partition/account for this user's short Palmetto GPU pilots; inspect
owner partitions with sinfo -a. Do not explicitly set QoS: Palmetto assigns it.
Normal work1 is not the only available profile. Verify allocated-GPU operation,
not just plugin presence on the login node. Update the platform guide from tested
scheduler behavior. Headless plotting must override MPLBACKEND before importing
Matplotlib, since an inherited notebook backend can be unavailable in conda.

CDKL5 reviewed WT pilots passed full 24-stage QC on Colab T4 (apo) and Palmetto
A100 (ATP/one-Mg). Keep runtime source snapshots and collection receipts intact.
GPU mutant coverage and biological convergence remain separate, untested claims.

Public-input CDKL5 WT pilots also passed full 24-stage QC and 100 ps production
plus analysis on compbio (run cdkl5_public_v1, runtime0d7f4b3). Keep this evidence
separate from the earlier reviewed-input tests. Public L119R is CPU-only.

CDKL5 visualization uses the activation project's full-kinase orientation by
default, with its network orientation optional. Store numeric camera metadata
only; never require private activation-study files. Fit native core C-alpha
positions to one WT reference and reuse its framing across panels. Read physical
time from matching production logs; never invent zero-time labels. Validate
declared mutations and distinguish starting-model grids from MD snapshots.
Optional small-protein tutorials/templates remain available outside the main
4BGQ walkthrough; preserve their historical validation records.

Colab user instructions use the local terminal, Colab CLI and its remote console.
Keep account/Drive authentication, committed snapshot deployment, recovery and
MODELLER prompting in that path. Do not reintroduce a web notebook setup route
or notebook-cell instructions. Internal runtime protocol and plotting-backend
compatibility code are separate from the documented user workflow.

Keep structure, simulation and analysis as unnumbered siblings within each state.
Write new replica reports to analysis/<run>/<variant>/<replica>, model views to
the variant's figures/, replica summaries to aggregate/, and variant comparisons
to analysis/<run>/compare. Explicit migration applies only to completed flat
state-based r1 reports; verify hashes, reject conflicts/symlinks/heavy data, keep
historical provenance/receipts and never move outputs during path resolution.
Historical numbered examples keep their layouts and template hashes.
