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
  Main README provides the overview and tutorial links. Complete scientific commands and workflows live
  in docs/tutorials/colab/ (trpcage_apo.md, kras_apo.md, kras_apo_colab.md, cdkl5_apo.md, cdkl5_atpmg.md). Keep one scientific workflow across platforms; no user YAML.
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

### Standard simulation and analysis architecture

```text
research/projects/<project-example>/
└── data/
    └── <state>/                              # apo or atpmg
        ├── structure/
        │   └── <target>/                     # Biological Target (UniProt) + Structural Template (PDB)
        │       ├── raw/                      # Untouched experimental inputs (cmda fetch)
        │       │   ├── <pdb_id>.cif          # Raw crystallographic archive record
        │       │   ├── template.pdb          # Cleaned monomer (waters/buffers stripped)
        │       │   └── uniprot.fasta         # Canonical sequence for biological numbering
        │       │
        │       ├── variants/                 # Modeled coordinates (cmda model)
        │       │   ├── WT/
        │       │   │   ├── WT.pdb            # Canonical modeled Wild-Type coordinate
        │       │   │   ├── target.B99990001.pdb  # Modeller build
        │       │   │   ├── target.B99990001_with_cryst.pdb # Model with crystal headers
        │       │   │   ├── target.ini        # Starting coordinate mesh
        │       │   │   ├── target.rsr        # Spatial restraints
        │       │   │   ├── target.sch        # Optimization schedule
        │       │   │   ├── target.D00000001  # Restraint progress log
        │       │   │   ├── target.V99990001  # Restraint violation profile
        │       │   │   ├── target.ali        # PIR target sequence alignment
        │       │   │   ├── alignment.ali     # PIR template pairwise alignment
        │       │   │   ├── <uniprot>.fasta   # Sequence reference
        │       │   │   ├── <pdb>_clean.pdb   # Sanitized template
        │       │   │   ├── <pdb>_orig.pdb    # Raw template backup
        │       │   │   ├── <pdb>.pdb         # Isolated chain coordinate
        │       │   │   ├── pipeline.log      # Modeller build log
        │       │   │
        │       │   └── <mutant>/             # Point mutant
        │       │       ├── <mutant>.pdb      # Mutated coordinate
        │       │       ├── WT_<mutant>.pdb   # Direct mutation coordinate
        │       │       ├── WT_<mutant>_with_cryst.pdb # Mutant with crystal headers
        │       │       ├── mutate.log        # Rotamer sampling & energy log
        │       │
        │
        ├── simulation/
        │   ├── WT/
        │   │   ├── prep/                     # Solvated Amber system (system.xml, model.pdb, solvated.pdb, system.prmtop, system.inpcrd)
        │   │   ├── eq/                       # 24-stage equilibration
        │   │   │   ├── eq.pdb                # Final relaxed coordinates
        │   │   │   ├── system.xml            # OpenMM system topology & parameters
        │   │   │   ├── state.xml             # OpenMM thermodynamic state snapshot
        │   │   │   └── stages/               # Granular restart checkpoints (state_01.xml ... state_24.xml)
        │   │   └── md/                       # Production replicas (e.g. 100 ns)
        │   │       ├── r1/                   # Production replica 1
        │   │       │   ├── md.dcd            # Continuous merged trajectory (full system)
        │   │       │   ├── md.log            # Continuous merged thermodynamic log
        │   │       │   ├── state.xml         # Latest OpenMM simulation state snapshot
        │   │       │   ├── prod.chk          # OpenMM checkpoint binary for resume
        │   │       │   ├── top.prmtop   # Solvent-stripped topology for PyMOL & analysis
        │   │       │   ├── traj.dcd # Continuous PBC-autoimaged stripped trajectory
        │   │       │   ├── topCa.prmtop         # Cα backbone topology
        │   │       │   ├── trajCa.dcd       # Cα backbone trajectory
        │   │       │   └── chunks/           # Raw crash-recovery chunks (0001.dcd, 0001.log, ...)
        │   │       ├── r2/                   # Independent velocity seed (triplicate)
        │   │       └── r3/                   # Independent velocity seed (triplicate)
        │   └── <mutant>/
        │       ├── prep/
        │       ├── eq/
        │       └── md/{r1, r2, r3}/
        │
        └── analysis/                         # Mirrored deliverables (Zero JSON policy)
            ├── WT/
            │   ├── eq/                       # Equilibration QC deliverables (qc.png, qc.txt)
            │   └── md/
            │       ├── r1/                   # Per-replica metrics (rmsd.png, rg.png, rmsf.png, summary.txt)
            │       ├── r2/
            │       └── r3/
            ├── <mutant>/
            │   ├── eq/                       # qc.png, qc.txt
            │   └── md/{r1, r2, r3}/
            └── compare/                      # Flat cross-variant comparison (rmsd.png, rg.png, rmsf.png, summary.txt)
```

- **Documentation & Writing Standard (Simple, Concise, and Direct)**:
  All documentation, guides, comments, and commit messages must be simple, concise, and direct:
  - **Zero Buzzwords & Tech Fluff**: Never use words like "hub", "nexus", "orchestrator", "ecosystem", or "VarMDyn".
  - **Direct Nouns Only**: Use `eq/` (24-stage equilibration), `md/r1/` (Production replica 1), `prep/` (Solvated Amber system).
  - **Step-by-Step Handover Context**: When documenting commands in guides, explicitly describe: (1) command run, (2) exact files generated, and (3) which subsequent step consumes those files.
- **100% Native Biophysical Architecture Contract (Zero-JSON & Zero-TSV Policy)**:
  Output trees (`data/`) must contain zero `.json` and zero `.tsv` receipts or audit clutter. Domain files (`.pdb`, `.fasta`, `.prmtop`, `.inpcrd`, `.xml`, `.log`, `.dcd`, `.chk`, `.csv`, `.png`) and verified stage logs are the sole sources of physical truth:
  - `structure/raw/` contains only `.cif`, `template.pdb`, and `uniprot.fasta` (no `residues.json`, `input.json`).
  - `structure/variants/` contains coordinate models (`.pdb`), native Modeller outputs (`.ali`, `.log`, etc.), with model verification in logs/ (no `model_request.json`, `model_check.json`).
  - `simulation/prep/` contains strictly native Amber files (`model.pdb`, `solvated.pdb`, `system.prmtop`, `system.inpcrd`, `build.in`, `leap.log`). No `prep.json` or `topology_check.json`.
  - `simulation/eq/` contains `eq.pdb`, `system.xml`, `state.xml`, restart states in `stages/state_*.xml`. Zero JSON receipts.
  - `simulation/md/` contains `md.log`, `npt.state.xml`, `system.xml`, `solvated.pdb`, `prod.chk`, stripped/autoimaged products (`top.prmtop`, `traj.dcd`), Cα backbone files (`topCa.prmtop`, `trajCa.dcd`), raw chunks in `chunks/`. Zero JSON receipts.
  - `analysis/` contains only scientific metrics (`.csv`), figures (`.png`), summary deliverables (`summary.txt`, `qc.txt`), with stage verification under logs/. Zero JSON files in `analysis/`.
  - **Native State Derivation & Gating**: Inter-stage validation, neutrality checks, and sequence identity gating must be derived dynamically and natively from the biological domain files themselves (e.g. parsing `system.prmtop` via `check_topology(..., save_receipt=False)`, verifying point mutations directly against `WT.pdb` coordinate atoms, parsing sequence from `uniprot.fasta`).
  - **Text Completion Markers**: State markers (`STRUCTURE_DONE`, `MODEL_DONE`, `ANALYSIS_DONE`) record operational metadata (e.g., `seed`, `pdb`, `uniprot`, `range`) as clean, human-readable key-value lines without writing duplicate JSON ledgers.
- **Equilibration (`eq`)**: Final relaxed coordinates are saved exclusively as
  `eq.pdb` (duplicate `equilibrated.pdb` eliminated). Intermediate restart states
  are saved into `stages/` as `state_01.xml` ... `state_24.xml`. All QC deliverables
  belong in `analysis/<variant>/eq/` as `qc.png` and `qc.txt`.
- **Production MD (`md/<replica>`)**: Processed trajectory/topology pairs are
  `traj.dcd`/`top.prmtop` and `trajCa.dcd`/`topCa.prmtop` across every state.
  The merged measurements are `md.log`; `prod.chk` remains the production
  checkpoint, with `system.xml`, `solvated.pdb` and `npt.state.xml` retained for restart.
  Raw crash-recovery chunks are strictly isolated in `chunks/{idx:04d}.dcd` and
  `chunks/{idx:04d}.log` using 4-digit zero-padding without units (`ns`/`ps`).
  Post-merge PBC-autoimaged and solvent-stripped trajectories (`top.prmtop`,
  `traj.dcd`) and Cα backbone (`topCa.prmtop`, `trajCa.dcd`) are generated
  via `cmda merge` for PyMOL inspection and biophysical analysis.
- **Analysis and comparison (`analysis/`)**: Analysis symmetrically mirrors simulation
  (`analysis/<variant>/eq/` for `qc.png` and `qc.txt`; `analysis/<variant>/md/<replica>/`
  for replica metrics). Cross-variant comparison lives flat in `analysis/compare/`
  (`summary.txt`, `rmsd.png`, `rg.png`, `rmsf.png`). Zero `.json` files in `analysis/`.

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
- Validate nonempty outputs before recording completion in stage logs. Run verification means software
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
- Deploy clean committed snapshots through cmda colab setup --local-source, record
  SNAPSHOT.json hashes and exclude Git metadata, build copies and private material.
- Run protein GPU pilots in a Palmetto allocation or on Colab, never Divergent.
  Prefer compbio for current pilots. Keep Colab jobs foreground with streaming
  output; use CPU-only sessions for Drive cleanup when possible.
- Stop dedicated sessions after verified retrieval; heavy trajectories stay on Drive.

## Verification and CI Standards

### Zero Red Cross Policy & Remote CI Verification
- **Every push triggers GitHub Actions CI (`Python CI`)**: The agent must NEVER push and walk away.
- **Mandatory post-push verification**: Immediately after running `git push`, the agent MUST check remote CI status using `gh run list --limit 1` or `gh run watch`.
- **Immediate remediation**: If GitHub Actions shows a red cross (`X`), run `gh run view <run-id> --log-failed` immediately to inspect the failure, fix the issue, run local verification, and push the fix until a green checkmark (`✓`) is confirmed. Never leave a failed CI run on GitHub.

### Mandatory Local CI Parity Gate
Before any commit or push, run the exact suite executed by `.github/workflows/ci.yml`:
1. **Black code formatting**: Run `black --check .`. If formatting fails, run `black .` to reformat.
   *(Note: A single line exceeding 88 characters will cause `black --check` to fail with exit code 123 under xargs in CI).*
2. **Ruff linting**: Run `ruff check .` and resolve all warnings and errors.
3. **Shell script syntax**:
   ```bash
   find examples scripts .githooks -type f \( -name '*.sh' -o -name '*.slurm' -o -name 'pre-commit' -o -name 'pre-push' \) -exec bash -n {} \;
   ```
4. **Unit and regression tests**: Run `python -m pytest tests/ -q` (all tests must pass).
5. **Gitleaks secret scan**: Run `.githooks/pre-commit` (scans staged snapshot with Gitleaks, Black, Ruff, and bash -n).

### Synthetic End-User Verification & Colab Pilot Protocol
When tasked with validating workflows on Google Colab as an end user:
1. **Clean-Room End-User Persona**: Act strictly as an external researcher. Follow the tutorial guide (`docs/tutorials/colab/*.md`) verbatim from top to bottom without using internal developer shortcuts or uncommitted local paths.
2. **Ephemeral VM Storage (`/content/`)**: Run directly in `/content/<project>` on the Colab VM. **Never mount Google Drive** for automated agent pilots. Google Drive persistence is intended for human researchers preserving 100 ns runs; automated pilots must avoid interactive browser OAuth redirects and Drive FUSE I/O throttling.
3. **Zero-Leak Credential Discipline**: Never embed license keys or credentials (such as `KEY_MODELLER`) as literal strings in bash commands, `colab exec` calls, command-line arguments, or repository files. For remote Modeller setup, use the mode-600 temporary file upload pattern (`cmda colab setup --session <name> --modeller`) with an immediate trap wipe, or model coordinates locally on the workstation and upload coordinates only. Keep the cloud GPU environment 100% keyless.
4. **Biophysical Scaled Verification**: Exercise 100% of the authentic pipeline stages (`fetch` $\to$ `model` $\to$ `prep` $\to$ `eq` $\to$ `md` $\to$ `merge` $\to$ `analyze`). Scale down duration (e.g. `--software-test` for equilibration, $0.1\text{ ns}$ to $1\text{ ns}$ for production MD) to verify all computational contracts, file handovers, and figure generation within minutes.
5. **Session Lifecycle Discipline**: Provision with `colab new -s <name> --gpu T4`, execute commands via `colab exec -s <name> "..."`, and **always terminate the VM immediately** upon completion (`colab stop -s <name>`) to conserve compute units.

Test project creation, conflict/symlink protection, argument forwarding, wheel resources and hidden MODELLER prompting.
GPU acceptance includes real SIGINT/resume, unchanged committed chunk hashes,
continuous frame times, finite measurements, native numbering and reviewed plots.
Do not claim a new GPU run when only CPU/template tests were performed.

## Documentation upkeep

README.md introduces the project, lists the three Trp-cage/KRAS tutorials with serial
numbers, outlines the workflow,
and keeps a collapsed project layout, changelog, acknowledgements and citations.
Keep CDKL5 guides available but unlisted in the root README until requested. Use
KRAS for the README project layout. Keep installation commands in each tutorial,
not duplicated in the root README.
Do not introduce a shared installation guide without an explicit user request.

Each tutorial is self-contained: local prerequisites, installer, activation in each
terminal, and browser authorization instructions precede the scientific workflow.
Show a compact output tree in a collapsed details block after each stage's
commands. Keep commands visible and use a short summary label for the tree. Put
the complete workflow tree in a collapsed details block after the workflow, before
acknowledgements, citations and the visitor counter. Verify filenames and downstream
handoffs against the implementation. Preserve scientific commands, parameters,
time conversions, checkpoint behavior and historical evidence when reorganizing.
State which machine runs each step. Preserve the root badge row and project
introduction; developer setup belongs in CONTRIBUTING.md.

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

Colab user instructions support the local terminal, Colab CLI and its remote console.
The KRAS Colab tutorial also offers the explicitly requested browser-only terminal
route: select a GPU and mount Drive in the browser, then run the same cloud setup
and scientific commands in its terminal. Local installation, CLI authentication,
downloads and desktop PyMOL are optional for that route. Keep recovery instructions
usable through either route; do not duplicate scientific procedures or add notebook
simulation cells. Preserve hidden MODELLER prompting and committed snapshot rules.
Internal runtime protocol and plotting-backend compatibility code stay separate
from the documented user workflow.

Keep structure, simulation and analysis as unnumbered siblings within each state.
Write new replica reports to analysis/<run>/<variant>/<replica>, model views to
the variant's figures/, replica summaries to aggregate/, and variant comparisons
to analysis/<run>/compare. Explicit migration applies only to completed flat
state-based r1 reports; verify hashes, reject conflicts/symlinks/heavy data, keep
historical provenance/receipts and never move outputs during path resolution.
Historical numbered examples keep their layouts and template hashes.

## Public release and credentials

Keep this canonical dev1 checkout private, with its public push URL disabled and
colabmda.privateCheckout=true. Enable .githooks; use Gitleaks default rules plus
MODELLER rules before commits and on all history before public pushes. Never
store actual key values in scanner configuration, output or tracked fixtures.
Publish only whitelisted clean-history snapshots from a separate release checkout.
AGENTS.md, PLAN.md, runtime data, lab account paths and private run records are
excluded from that public snapshot. Keep private research provenance unchanged.
Colab source export accepts any clean committed branch or tag. Package all seven
nucleotide resources in wheels; verify an isolated wheel import. Public HPC
adapters leave partition/account selection to the user; private execution-project
adapters retain account-specific settings. Preserve the README design/scientific
protocol. License installation reads the key from the environment, uses private
permissions and does not place key values in subprocess arguments.

### Required checks for every commit and release

- Before committing, verify `core.hooksPath=.githooks`, both hooks are executable,
  and the configured Gitleaks binary is available. Keep the private checkout's
  push URL disabled and `colabmda.privateCheckout=true`. Never use `--no-verify`,
  unset these guards or weaken scanner rules to make a commit/push pass.
- Verify CI parity before push and inspect remote CI results post-push:
  execute Black, Ruff, bash syntax checks, and pytest locally. After pushing, run
  `gh run list --limit 1` or `gh run watch` to confirm a green checkmark (`✓`). Any red
  cross (`X`) on GitHub must be diagnosed with `gh run view <id> --log-failed` and fixed immediately.
- Review staged paths and scan the exact staged snapshot through the pre-commit
  hook, including initial commits. Pay particular attention to notebook output,
  copied terminal sessions, configuration, fixtures and generated archives.
  Do not print raw suspicious diffs; use fully redacted scanner output.
- Keep default Gitleaks rules plus the custom MODELLER and Google authentication
  cookie rules enabled. Tests use synthetic credentials or explicit placeholders,
  never historical leaked values. New credential-handling code must keep hidden
  prompts, avoid argument/log exposure and preserve private file permissions.
- A public release uses a separate reviewed checkout with clean public ancestry.
  Scan all reachable refs using the command below, then scan the actual source
  export and unpacked wheel with the same configuration.
  Verify exclusions and scan a fresh clone of the published refs after release.
  Record commit IDs, scan scope and results without credential values.
- Keep secret scanning, push protection and the CI credential scan enabled.
  They supplement local checks; a green scan is not proof that credentials were
  revoked or that external copies were removed. Do not claim checks not performed.
- Historical cleanup bundles and reports stay restricted in the owning execution
  project's ignored data tree. Never include them in source snapshots, releases,
  Colab/HPC deployments or public issue attachments. Report exposure safely and
  leave provider account actions to the owner; do not test leaked credentials.

```bash
gitleaks git --log-opts=--all --redact=100 --config .gitleaks.toml .
```

## Processed trajectories and desktop movies

New merge outputs use `traj.dcd`, `top.prmtop`, `trajCa.dcd`, `topCa.prmtop`
and one merged `md.log` across all states, variants and replicas. Readers retain
legacy filename support. `cmda merge --rename-files` validates completed pairs,
refuses conflicts, symlinks and incomplete chunks, and renames without rewriting
scientific bytes. Never rename an active run or its raw chunks/checkpoints.
`cmda view --frame all` must find actual merged trajectories even without an
explicit run ID, generate a matching viewer PDB and start playback. Document one
primary movie command and a collapsed manual alternative. Viewer provenance is
plain-text `visualize.log`, not a generated JSON receipt. Keep complete tutorial
trees before acknowledgements, citations and the visitor footer.

## Completion logs — approved marker replacement

Do not create standalone completion flags. Stage metadata and output verification
are stored in plain-text stage logs under logs/, mirroring the owning data folder.
Use StageLog for completion checks; an existing log alone does not prove success:
its final status and current output hashes must pass validation. Invalidating a
stage retains metadata in its log and prevents stale legacy flags from bypassing
the current status. Legacy flags are readable for compatibility; cmda logs
--migrate preserves their metadata and removes them after verifying outputs.
Visual provenance also belongs in these logs. Native MD measurements and trajectory
chunks stay paired in data/. Preserve raw scientific bytes during log migration.
The README and current tutorials show logs/ rather than standalone flags.
