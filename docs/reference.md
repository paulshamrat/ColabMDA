# Advanced reference

The main repository README owns installation, project creation and the complete
Trp-cage example. This page covers the implemented CLI, compatibility behavior
and visualization. Scientific settings and recovery semantics are in
[the protocol](protocol.md).

## Commands and scope

`cmda` and `colabmda` invoke the same CLI. Use each command's `--help` for its
available arguments; the table describes the current development checkout,
not a newly announced public release.

| Operation | Command | Behavior |
| :--- | :--- | :--- |
| Example creation | `cmda stage --example trpcage --root PROJECT` | Installs seven scripts with source hashes; preserves different existing files. |
| Existing state layout | `cmda stage --state TAG --variant NAME --root PROJECT` | Compatibility layout and runner generation for an existing research system. |
| Modeling | `cmda modeller build` | Homology modeling from a supplied template/UniProt selection; MODELLER key required. |
| Mutation | `cmda modeller mutate` | Mutates a supplied model while retaining its biological residue mapping. |
| Preparation | `cmda prep --state TAG --variant NAME --root PROJECT` | Assembly and Amber parameterization of reviewed inputs. |
| Equilibration | `cmda sim --state TAG --stage eq --protocol staged --root PROJECT` | EM/NVT/NPT/QC; review physical QC separately. |
| Production | `cmda sim --state TAG --stage md --replica r1 --root PROJECT` | Checkpointed production with an absolute target duration. |
| Merge | `cmda merge --state TAG --variant NAME --root PROJECT` | Autoimaging and solvent stripping with matching topologies. |
| Analysis | `cmda analyze --state TAG --variant NAME --root PROJECT` | RMSD/Rg/RMSF tables, native residue labels and downstream plots. |
| Comparison | `cmda openmm compare --series LABEL=DIR1,DIR2 --outdir OUT` | Aggregates supplied analysis tables and writes PNGs under OUT/figures. |
| Viewing | `cmda openmm view --pdb-dir DIR` | Generates a PyMOL view for a supported PDB/trajectory pair. |
| Snapshots | `cmda snapshots --config FILE --output FILE --temp-dir DIR` | Uses the custom snapshot description with the existing PyMOL renderer. |
| Remote setup | `cmda colab setup --local-source -s SESSION` | Deploys a committed development snapshot and runs the single installer. |
| Remote execution | `cmda colab exec -s SESSION 'COMMAND'` | Foreground execution with streamed output and checked completion. |

The top-level `cmda compare` currently reports the requested comparison; it does
not generate plots. Use the implemented `cmda openmm compare --series ...` command.
The standard analysis pipeline does not calculate PCA or contact maps. Docking
is not part of this active source branch or its maintained documentation.

## Modeling inputs

Build accepts a template PDB ID, UniProt ID, chain and optional residue range.
Choose those from the scientific question rather than copying a project-specific
kinase or KRAS selection. Mutation accepts an existing PDB, chain and mutation
or mutation list. Check native residue identity and coordinate presence before
using the output. Missing loops and cofactors require explicit review.

MODELLER runs in its optional environment. Its key prompt, automated key handling
and installation are documented in the main README; no real keys belong in
scripts, project records or documentation.

## Topology, timing and analysis inputs

A trajectory must be loaded with a topology containing the same atoms in the same
order. Use the solute trajectory with its stripped topology, and the C-alpha stream
with its C-alpha topology. An unstripped solvated PDB is not interchangeable with
an already stripped trajectory. Check atom/frame counts after merging.

For direct lower-level Amber analysis, supply a native-numbered reference:

```bash
cmda openmm analysis --pdb-dir PROJECT/data/02_simulation/RUN/apo/WT/03_md/r1 \
  --topology PROJECT/data/02_simulation/RUN/apo/WT/03_md/r1/stripped.prmtop \
  --trajectory PROJECT/data/02_simulation/RUN/apo/WT/03_md/r1/prod_keep_ligand.dcd \
  --reference PROJECT/data/02_simulation/RUN/apo/WT/01_prep/model.pdb \
  --interval 10 --outdir PROJECT/data/02_simulation/RUN/apo/WT/04_analysis
```

Replace PROJECT and RUN with actual paths. The high-level `cmda analyze` supplies
the model reference automatically. Its current route uses a 10 ps interval;
for a different sampling interval use the lower-level analysis argument and
check logged frame times. Native labels must not be inferred from Amber indices.

Comparison consumes existing `rmsd.csv`, `rg.csv` and `rmsf.csv` tables. Group only
replicas with compatible coordinates, residue identities, protocols and time axes.
A comparison figure does not establish biological convergence or statistical
independence. Review the supplied tables and the definition of the mean/SD bands.

## PyMOL and snapshot inputs

Use a separate environment with PyMOL, Pillow and the tool installed for local
viewing/rendering. Keep the computational environment and its checkpoints intact.
The view command accepts explicit `--topology`, `--trajectory` and `--resi` values
when automatic path detection does not match a project's layout.

The optional snapshot description is JSON for that rendering operation; it is
not a project configuration file required by the runnable example. Its fields are:

| Field | Meaning |
| :--- | :--- |
| `camera_view` | An explicit 18-value PyMOL view matrix, reviewed for the selected structure. |
| `align_ref_pdb` | Optional reference PDB for consistent inter-system alignment. |
| `stable_core_sel` | Explicit selection of homologous stable-core atoms. |
| `systems` | Mapping of system names to their rendering inputs. |
| `systems.NAME.pdb` / `.dcd` | Matching topology and trajectory paths. |
| `systems.NAME.states` | Frame/state indices to render. |
| `systems.NAME.times` | Display labels corresponding to those states. |
| `systems.NAME.mut_residue` | Residue highlighted by the helper. |

Frame/state numbering and physical time are different: derive labels from the
actual trajectory sampling and production origin. Use the same reviewed camera,
reference and homologous selection for panels intended for direct comparison.
Do not copy a camera matrix or a residue selection from an unrelated protein.

The existing renderer retains KRAS-specific styling, including residues 57–75.
A custom input description does not remove all of those assumptions. Inspect
its selections and rendered output before using it for another protein; it is
not a general protein-specific annotation engine. The Trp-cage example does not
require this helper. Numerical figure outputs use 7-inch, 300-DPI PNGs under
project `figures/` directories; review structural compositions separately.

## Existing projects and developer guidance

Legacy path resolution remains for existing projects. New example projects use
numbered `data/01_structure`, `data/02_simulation/<run>` and matching logs, without
project YAML. Keep old run provenance intact and use migration inventories to
locate moved files. Historical instructions remain available in Git history.

CONTRIBUTING.md and AGENTS.md at the repository root contain the development
architecture, formatting, verification and release rules. The working example's
`tests/gpu.py` exercises its generated scripts with an actual interruption/resume
and downstream checks. Unit tests and a successful installer are not substitutes
for a GPU integration run or scientific validation.
