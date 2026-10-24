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
- `examples/trpcage/scripts/` contains runnable user workflows, not engine code.
  Main README owns installation/example/Colab/restart instructions. Do not add
  separate example Markdown guides or project configuration YAML.
- `cmda stage --example trpcage --root ...` installs templates and records hashes
  in data/00_setup/example.json. Ship templates in wheels and source snapshots.
  Refuse symlinks and different existing files; never overwrite user settings.
- One working project: research/projects/colabmda-example on Divergent and Drive.
  Its numbered scripts are copied templates; tests exercise those same scripts.
  Runtime paths are data/01_structure and data/02_simulation/<run>/, with
  logs/02_simulation/<run>/. Preserve existing historical provenance verbatim.
- Never vendor external software or put coordinates/trajectories in src/examples.
  Scientific inputs, caches, reports, checkpoints and figures stay in ignored data.
  Compiler decks run in temporary sandboxes, not scientific data folders.

## Scientific and recovery invariants

- Preserve native author numbering/insertion codes. Structure retrieval supports
  mmCIF and PDB. Check sequence identity and C-alpha coordinate presence.
- Trp-cage is a designed peptide; verify its deposited sequence, model 1 and
  author chain A. A canonical UniProt mutation annotation is not applicable.
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
- Run protein GPU pilots on Colab, not Divergent. Keep Colab jobs foreground with
  streaming output. Use CPU-only sessions for Drive cleanup when possible.
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

Keep CONTRIBUTING and reference pages consistent with the flat source layout and
main README. The installation reference links to the README rather than duplicating
a second setup workflow. Document implemented behavior: standard analysis is
RMSD/Rg/RMSF; comparison plotting currently uses openmm compare --series. Legacy
path handlers and historical provenance are intentional compatibility support.
