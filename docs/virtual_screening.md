# Docking and Virtual Screening Guide

This guide shows the lightweight docking layer used after MD analysis identifies a candidate pocket. The typical loop is:

```text
MD protein/cofactor system → inspect pocket → dock/screen small molecules → select poses → MD protein-ligand complexes
```

Docking is not a replacement for MD. It is a fast ranking and pose-generation step. The selected docked poses should be checked chemically, then re-run with the protein-ligand MD workflow.

The practical pipeline is:

```text
1. Build protein-only WT and variant structures
2. Build native cofactor complexes such as 4ldj_gdp and 4ldj_g12d_gdp
3. Run and analyze MD for the protein/cofactor states
4. Use the selected receptor state and pocket/region for docking
5. Merge and inspect docking rankings
6. Promote only selected docked poses into data/str/<pdb>/docked/
7. Run new protein-ligand MD systems for those selected poses
```

This mirrors the lab-style workflow: docking narrows the chemical search space, while MD checks whether selected poses are physically stable in a solvated protein system.

---

## 1. Environment

Docking uses a separate environment from MD:

```bash
source "$HOME/miniforge3/etc/profile.d/conda.sh" && conda activate docking_env
```

If you are bootstrapping a fresh Colab session and want docking tools installed:

```bash
WITH_DOCKING=1 WITH_LIGAND=1 bash ./bootstrap_colab_openmm_gpu.sh 2>&1 | tee /content/bootstrap_colabmda.log
```

For local-branch testing through `colab exec`:

```bash
python scripts/colab_exec_bootstrap.py --session kras-sim --from-local --with-docking --with-ligand --timeout 7200
```

This installs **RDKit** for molecule handling, **Meeko** for PDBQT preparation and SDF export, and **AutoDock Vina** for docking in `docking_env`. Vina is CPU-based, which is fine for small tests and batch screens; the GPU remains most important for the later MD step in `openmm_env`.

---

## 2. Folder Layout

Run commands from your Drive workspace, usually:

```bash
cd /content/drive/MyDrive/ColabMDA
```

Recommended layout:

```text
data/
  str/                 <- reusable structure inputs
  sim/                 <- MD simulations, named by biological system
  analysis/            <- MD analysis and comparison plots
  docking/             <- docking campaign workspaces, usually split into shards
```

`data/` is gitignored, so ligand libraries, docking poses, and MD trajectories stay out of the code repository.

For the current 4LDJ/GDP workflow, keep the system names consistent from structure
building through docking and MD:

```text
data/
  str/
    4ldj/
      raw/
        4ldj_orig.pdb

      proteins/
        4ldj_wt/
          4ldj_wt.pdb
        4ldj_g12d/
          4ldj_g12d.pdb

      cofactors/
        gdp.sdf

      complexes/
        4ldj_gdp/
          4ldj_gdp.pdb
          gdp.sdf
        4ldj_g12d_gdp/
          4ldj_g12d_gdp.pdb
          gdp.sdf

      docked/
        4ldj_wt/
          lig001_pose.sdf
        4ldj_g12d/
          lig001_pose.sdf
        4ldj_gdp/
          lig001_pose.sdf
        4ldj_g12d_gdp/
          lig001_pose.sdf

  sim/
    4ldj_wt/
      r1/
    4ldj_g12d/
      r1/
    4ldj_gdp/
      r1/
    4ldj_g12d_gdp/
      r1/
    4ldj_gdp_lig001/
      r1/

  docking/
    4ldj_gdp/
      receptor/
      source/
      shards/
        000001/
          library/
          ligands/
          results/
      merged/
        ranked_hits.csv
        report/

  analysis/
    4ldj_gdp/
      r1/
    4ldj_gdp_lig001/
      r1/
    compare/
```

The difference between `docking/` and `docked/` is intentional:

- `data/docking/` is the full screening campaign: source libraries, shard folders,
  Vina logs, rankings, and reports.
- `data/str/4ldj/docked/` contains only the selected docked SDF poses promoted to
  structure inputs for the next MD round.

The GDP/Mg cofactor is a native coordinate transfer, not a docking result: extract it
from the original crystal structure while preserving coordinates, align the modeled WT
or mutant protein back to that crystal frame, and then combine the cofactor with the
protein for `4ldj_gdp` or `4ldj_g12d_gdp`.

For production MD, native nucleotide/cofactor systems such as GDP/Mg should use a
validated cofactor parameterization route. The simple docking-to-SDF-to-GAFF handoff is
intended for ordinary screened small molecules that pass chemical inspection, not for
claiming curated nucleotide/cofactor chemistry.

---

## 3. Build or Provide a Ligand Library

For a real screen, download molecules from a named source such as ZINC, ChEMBL, or PubChem, then save the exact SMILES/SDF file and query metadata with the project. ColabMDA does not silently choose a database for you.

For a small local test only:

```bash
colabmda dock library --outdir data/docking/4ldj_wt/shards/000001/library --preset demo
```

For your own SMILES table:

```bash
colabmda dock library --smiles-csv data/docking/4ldj_wt/source/zinc_subset_000001.csv --outdir data/docking/4ldj_wt/shards/000001/library
```

The CSV can contain columns named `smiles` and `id`. The command writes:

```text
data/docking/4ldj_wt/shards/000001/library/library.sdf
data/docking/4ldj_wt/shards/000001/library/library_manifest.csv
data/docking/4ldj_wt/shards/000001/library/library_metadata.json
```

For million-scale screens, split the source library into shard folders such as `000001`, `000002`, `000003`, etc. Each shard can run independently on a different Colab session or CPU node, then the ranked results can be merged later.

---

## 4. Prepare the Receptor and Docking Box

Use a cleaned receptor PDB from your structure or MD workflow. The box center and size are in Ångström. Choose them from the pocket you identified during analysis, for example by inspecting a cofactor/ligand centroid in PyMOL.

```bash
colabmda dock prep-receptor --receptor data/str/4ldj/proteins/4ldj_wt/4ldj_wt.pdb --name 4ldj_wt --center 10.0 20.0 30.0 --size 22.0 22.0 22.0 --outdir data/docking/4ldj_wt/receptor
```

This writes:

```text
data/docking/4ldj_wt/receptor/4ldj_wt.pdbqt
data/docking/4ldj_wt/receptor/vina_box.txt
data/docking/4ldj_wt/receptor/receptor_metadata.json
```

---

## 5. Prepare Ligands

Convert the 3D SDF library into one PDBQT file per ligand:

```bash
colabmda dock prep-ligands --library-sdf data/docking/4ldj_wt/shards/000001/library/library.sdf --outdir data/docking/4ldj_wt/shards/000001/ligands
```

This uses Meeko so ligand bond orders and formal charges come from the SDF rather than being guessed from PDB.

---

## 6. Run Docking

Run Vina over the prepared ligand folder:

```bash
colabmda dock run --receptor-pdbqt data/docking/4ldj_wt/receptor/4ldj_wt.pdbqt --ligands-dir data/docking/4ldj_wt/ligands/pdbqt --config data/docking/4ldj_wt/receptor/vina_box.txt --outdir data/docking/4ldj_wt/results --exhaustiveness 8 --num-modes 9
```

For the sharded layout:

```bash
colabmda dock run --receptor-pdbqt data/docking/4ldj_wt/receptor/4ldj_wt.pdbqt --ligands-dir data/docking/4ldj_wt/shards/000001/ligands/pdbqt --config data/docking/4ldj_wt/receptor/vina_box.txt --outdir data/docking/4ldj_wt/shards/000001/results --exhaustiveness 8 --num-modes 9
```

Important outputs:

```text
data/docking/4ldj_wt/results/ranked_hits.csv
data/docking/4ldj_wt/results/docking_metadata.json
data/docking/4ldj_wt/results/poses_pdbqt/
data/docking/4ldj_wt/results/logs/
```

The score is Vina’s docking score in kcal/mol. More negative is better within the same receptor, box, scoring function, and ligand preparation protocol. Do not compare scores across unrelated setups as if they are experimental binding free energies.

For multiple shards, merge the shard rankings:

```bash
colabmda dock merge-results --campaign-dir data/docking/4ldj_wt
```

Then make a plot and short report:

```bash
colabmda dock report --results data/docking/4ldj_wt/merged --top-n 20
```

This writes:

```text
data/docking/4ldj_wt/merged/ranked_hits.csv
data/docking/4ldj_wt/merged/report/docking_scores.png
data/docking/4ldj_wt/merged/report/docking_report.md
data/docking/4ldj_wt/merged/report/top_hits.csv
```

---

## 7. Export Selected Poses for MD

Export the top-ranked poses back to SDF:

```bash
colabmda dock export-top --results data/docking/4ldj_wt/merged --outdir data/str/4ldj/docked/4ldj_wt --top-n 3
```

Then inspect the pose visually. If the chemistry and placement look sensible, run MD:

```bash
source "$HOME/miniforge3/etc/profile.d/conda.sh" && conda activate openmm_env
colabmda openmm run --name 4ldj_wt_top1 --replica r1 --ligand data/str/4ldj/docked/4ldj_wt/01_benzene_pose.sdf --total-ns 10.0 --traj-interval 10 --checkpoint-ps 1000
```

For publication-level protein-ligand simulations, keep the exact ligand source, protonation/tautomer assumptions, docking box, docking version, and MD force-field choices with the results.

---

## 8. What Is Happening Under the Hood?

1. **Library step:** RDKit reads SMILES or uses the tiny demo list, adds hydrogens, generates 3D conformers, and writes SDF.
2. **Receptor preparation:** Meeko writes a receptor PDBQT and the Vina box file.
3. **Ligand preparation:** Meeko converts each SDF record into a PDBQT ligand with torsions and atom types.
4. **Docking:** AutoDock Vina docks each ligand independently and writes docked PDBQT poses plus logs.
5. **Export:** Meeko converts selected docked poses back to SDF so ColabMDA’s OpenMM protein-ligand workflow can use them.

References: AutoDock Vina documents receptor/ligand preparation, docking, batch screening, and SDF export with Meeko; Meeko documents ligand/receptor PDBQT preparation and SDF export; RDKit documents molecule handling and 3D conformer generation.
