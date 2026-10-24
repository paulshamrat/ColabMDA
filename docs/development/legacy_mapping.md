# CLI compatibility reference

This describes the current development checkout, not an announced release.
`colabmda` and `cmda` are aliases for the same CLI. The main README owns the
complete runnable example; this page helps maintain existing project workflows.

| Operation | Current command | Scope |
| :--- | :--- | :--- |
| Modeling | `cmda modeller build` or `mutate` | MODELLER-backed modeling; a valid key is required. |
| Complete example | `cmda stage --example trpcage --root PROJECT` | Installs the seven workflow templates without overwriting different files. |
| Existing state layout | `cmda stage --state TAG --variant NAME --root PROJECT` | Compatibility layout and runner generation for an existing research system. |
| Preparation | `cmda prep --state TAG --variant NAME --root PROJECT` | Assembly and Amber parameterization from supplied reviewed inputs. |
| Equilibration | `cmda sim --state TAG --stage eq --root PROJECT` | Select the staged or quick protocol explicitly; inspect physical QC. |
| Production | `cmda sim --state TAG --stage md --replica r1 --root PROJECT` | Checkpointed production; preserve inputs and runtime to resume. |
| Merge | `cmda merge --state TAG --variant NAME --root PROJECT` | Autoimaging and solvent stripping with matching topologies. |
| Analysis | `cmda analyze --state TAG --variant NAME --root PROJECT` | RMSD, Rg and native-numbered RMSF tables, then downstream figures. |
| Comparison plots | `cmda openmm compare --series LABEL=DIR1,DIR2 --outdir OUT` | Aggregate the supplied analysis tables; PNGs are under OUT/figures. |
| Structure viewing | `cmda openmm view --pdb-dir DIR` | PyMOL viewing for a supported topology/trajectory pair. |

The top-level `cmda compare` currently reports the requested comparison but does
not generate plots. Use the implemented `cmda openmm compare --series ...` command.
The standard analysis pipeline does not calculate PCA or contact maps.

Direct `cmda openmm analysis` with an Amber topology requires `--reference` to
recover author residue labels. Legacy path resolution remains for existing
projects; new example projects use the numbered layout documented in the README.
