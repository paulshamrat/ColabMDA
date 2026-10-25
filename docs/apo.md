# CDKL5 protein in water

This guide uses CDKL5 (UniProt **O76039**), chain A of **4BGQ**, residues
**1–303**, with WT and **L119R** as the worked pair. It retains the full
research equilibration protocol. A short GPU validation does not establish
mutation effects or biological convergence.

The [main README](../README.md) remains the full baseline guide. This file is a
command-by-command CDKL5 workflow. The companion [ATP/Mg guide](atpmg.md) uses
the same folder design. Neither guide requires project YAML or generated runners.

## 1. Choose the execution platform

Follow [Google Colab](platforms/colab.md) for a Drive-backed GPU runtime, or
[HPC](platforms/hpc.md) for software setup, scratch storage and Slurm execution.
Both run the same scientific `cmda` commands below. Install MODELLER only when
building or mutating models; its installer retains the optional hidden key prompt.

The commands below use the Colab project path. On HPC set `PROJECT_ROOT` to the
scratch project and activate your installed environments instead. The platform
guide explains the scheduler wrappers; do not run dynamics on a login node.

## 2. Project setup

Use the shared [simulation project layout](../README.md#22-simulation-project).
Both states use the same folders; only variant processing steps are numbered.

**Inside Colab**, choose an empty project for this public-input walkthrough. If
the example folder contains older inputs, choose another name. Define the project
once for the current shell:

```bash
source /root/miniforge3/etc/profile.d/conda.sh
conda activate colabmda
export COLABMDA_TOOL=/content/ColabMDA
export PROJECT_ROOT=/content/drive/MyDrive/research/projects/colabmda-example
export COLABMDA_PROJECT_DIR="$PROJECT_ROOT"
export COLABMDA_RUN_ID=cdkl5_public_v1
mkdir -p "$PROJECT_ROOT"
cd "$PROJECT_ROOT"
```

On a local installation, change `PROJECT_ROOT` to your local project path.
Divergent can fetch, model, prepare and analyze; run protein GPU dynamics on Colab
or a suitable GPU machine. Export the same project and run variables again after
opening a new shell.

## 3. Structure, numbering and variants

Start from public 4BGQ and UniProt O76039. This workflow does not need another
project, a reviewed study model or manually copied coordinates.

Fetch the public structure and UniProt sequence:

```bash
cmda fetch --state apo --pdb-id 4BGQ --uniprot-id O76039 --range 1 303 --mutations L119R
```

Fetch keeps the raw mmCIF archive, with PDB fallback. The mapping records native
author positions and insertion codes. In 4BGQ A, 276 residues have coordinates and
27 are missing within UniProt 1–303. The requested WT amino acid must match both
UniProt and its mapped coordinates, with a C-alpha present.

Activate the MODELLER environment, install this same tool checkout there once,
then build WT and make L119R:

```bash
conda activate modeller_env
python -m pip install -e "$COLABMDA_TOOL"
cd "$PROJECT_ROOT"
cmda model --state apo
cmda model --state apo --variant L119R
conda activate colabmda
```

The commands use the fetched template and sequence, save models under
`data/apo/structure/models/{WT,L119R}/`, and check the full sequence, numbering
and C-alpha coordinates before `MODEL_DONE`. Derived models use UniProt 1–303;
raw coordinates stay unchanged. Seed 13031 is recorded and applied to modeling.
Completed models are reused only when inputs and settings match.

MODELLER requires a license. Installation keeps its optional question and hidden
key prompt; no key belongs in a README command. The lower-level `cmda modeller
build` and `mutate` commands remain available for custom inputs.

Review alignment, modeled loops, stereochemistry and scores. Exact sequence and
numbering checks do not establish biological accuracy. A new model is distinct
from the earlier activation-study pilot.


## 4. Prepare the system

Start with WT. After its short validation succeeds, repeat with `VARIANT=L119R`.
The run ID groups variants; each variant has independent preparation and replicas.

```bash
conda activate colabmda
export STATE=apo
export VARIANT=WT
cd "$PROJECT_ROOT"
cmda stage --state "$STATE" --variant "$VARIANT"
cd "$PROJECT_ROOT/data/$STATE/simulation/$COLABMDA_RUN_ID/$VARIANT/01_prep"
cmda prep --state "$STATE" --variant "$VARIANT"
```

`cmda stage` creates folders; direct commands perform the work. It creates no
runner scripts unless the compatibility option `--scripts` is requested.

Default Amber preparation uses **ff19SB, OPC, 14 Å padding**, neutralization and
an estimated **0.15 M NaCl** pair count. The salt count uses an estimated volume;
inspect the resulting box and ion count rather than claiming exact final molarity.
Inspect protonation, termini, modeled loops and charge before dynamics. This
route does not perform automatic pKa titration.

`topology_check.json` verifies neutral total charge, expected solute composition,
actual water-matched **Li–Merz 12–6** ion LJ values and nonempty outputs.
`PREP_DONE` is written after these checks. `preparation_input.json` freezes the
receptor, mapping, state and builder hashes. Identical completed preparation is
reused; changed inputs require a fresh run ID. `--no-amber` does not certify a
prepared simulation.

## 5. Equilibration, QC and production

### 5.1. Run all 24 stages

```bash
cd "$PROJECT_ROOT/data/$STATE/simulation/$COLABMDA_RUN_ID/$VARIANT/02_eq"
export COLABMDA_REQUIRE_GPU=1
cmda sim --state "$STATE" --variant "$VARIANT" --stage eq --protocol staged
```

Do not set `COLABMDA_TEST_1_24` for research. The shortened internal test and
`--protocol quick` require `--software-test`; they do not validate this schedule.

Restraints act on solute heavy atoms; bulk water and monovalent ions remain
unrestrained. Strengths below are in kcal mol⁻¹ Å⁻²:

| Stage | Operation | Conditions | Restraint |
| :--- | :--- | :--- | :--- |
| EM | Initial minimization | No physical time | None |
| 01–03 | Minimization ladder | Constant box | 100, 50, 25 |
| 04–09 | Heating ladder | Six 10 ps NVT intervals, 50–300 K, 1 fs | 25 |
| 10–14 | Minimization ladder | Constant box | 25, 10, 5, 2, 1 |
| 15–18 | NVT release ladder | Four 100 ps intervals, 300 K, 2 fs | 10, 5, 2, 1 |
| 19–22 | NPT release ladder | Four 100 ps intervals, 300 K, 1 bar, 2 fs | 10, 5, 2, 1 |
| 23 | Unrestrained NPT | 100 ps, 300 K, 1 bar, 2 fs | None |
| 24 | Final unrestrained NPT | 1000 ps, 300 K, 1 bar, 2 fs | None |

Numbered state/JSON records support recovery from completed stages. An explicitly selected quick
software test instead uses one short NVT and NPT interval; it is not this
full protocol. The legacy combined `openmm run` adapter retains ff14SB/TIP3P defaults
and delegates to quick equilibration; it should not be confused with the staged route.

### 5.2. Review QC before production

The research command fails when QC fails. It checks all 24 stage pairs, finite
stage energies and cumulative physical time of **1,960 ps**, alongside final
mean temperature within 10 K of 300 K, mean density 0.90–1.10 g/mL and absolute
final density slope at most 1e-4 g/mL/ps. Read `equilibration_qc.json` and
`figures/equilibration_qc.png`. Passing these checks does not establish biological
convergence; inspect the protein and cofactors as well.

Shared dynamics use PME, 1.0 nm cutoff, switching from 0.8 nm, H-bond constraints,
rigid water, LangevinMiddle friction 1/ps and an isotropic 1 bar barostat every
25 steps in NPT. Bulk water and monovalent ions are unrestrained; retained solute
heavy atoms, including a bound Mg where present, follow the restraint ladder.

### 5.3. Run a short production validation

After the WT equilibration passes and its plots are reviewed:

```bash
cd "$PROJECT_ROOT/data/$STATE/simulation/$COLABMDA_RUN_ID/$VARIANT/03_md"
cmda sim --state "$STATE" --variant "$VARIANT" --stage md --replica r1 \
  --protocol staged --total-ns 0.1 --checkpoint-ps 100 --traj-interval 10
```

The command reruns strict equilibration QC before production. This 100 ps target
checks the complete software route after full equilibration. It is too short to
interpret mutation effects. For a reviewed scientific campaign, choose a justified
duration and independent replicas; for example `--replica r1,r2,r3 --total-ns 100`.
That is 100 ns per replica, processed sequentially by this command. It is not an
automatic instruction to launch a full campaign.

### Timestep, Steps & Saved Output

The **integration timestep**, **frame-saving interval** and **checkpoint interval**
are different settings:

| Quantity | Current implementation | Meaning |
| :--- | :--- | :--- |
| Integration timestep | 1 fs during staged heating; 2 fs during subsequent dynamics, quick equilibration and production | Physical time advanced by one integration step. |
| Frame interval | `--traj-interval`, in ps; example uses 10 ps | How often coordinates and production measurements are saved. |
| Checkpoint/chunk interval | `--checkpoint-ps`, in ps; example uses 100 ps | Completed interval published before the production checkpoint advances. |
| Production target | `--total-ns`, in ns, per replica | Total requested production time, excluding equilibration. |

Unit conversions are **1 ps = 1,000 fs** and **1 ns = 1,000 ps = 1,000,000 fs**.
At the production timestep of 2 fs = 0.002 ps:

$$N_\mathrm{steps} = \frac{T_\mathrm{ns}\times 10^6}{\Delta t_\mathrm{fs}}$$

$$N_\mathrm{steps/frame} = \frac{\Delta t_\mathrm{frame,ps}\times 10^3}{\Delta t_\mathrm{fs}},\qquad
N_\mathrm{steps/chunk} = \frac{\Delta t_\mathrm{chunk,ps}\times 10^3}{\Delta t_\mathrm{fs}}$$

| Duration or interval | Calculation at 2 fs | Integration steps |
| :--- | :--- | ---: |
| 1 ps | 1,000 / 2 | 500 |
| 10 ps frame interval | 10,000 / 2 | 5,000 |
| 100 ps equilibration interval or chunk | 100,000 / 2 | 50,000 |
| 1 ns production | 1,000,000 / 2 | 500,000 |
| 100 ns production | 100,000,000 / 2 | 50,000,000 |

A 1 ns run with 10 ps frames and 100 ps chunks therefore has **100 frames**,
**10 completed chunks** and **10 frames per chunk**. Its first saved production
frame is at 10 ps, not at 0 ps. A 10 ps heating interval at 1 fs instead requires
10,000 steps. Minimization iterations do not represent physical simulation time.

The normal staged schedule integrates 60 ps heating + 400 ps NVT + 1,500 ps NPT
= **1,960 ps (1.96 ns)** before production. The quick example uses 100 ps NVT +
100 ps NPT = **200 ps (0.2 ns)** before its 1 ns production. Equilibration is not
part of the production target. Three 100 ns replicas provide 300 ns aggregate
production, not one continuous 300 ns trajectory.

Production currently fixes the timestep at 2 fs; there is no production timestep
CLI option. Changing `--traj-interval` changes saved output, not the integrator.
The `--steps` option is a software-testing override; use `--total-ns` for research production. Use intervals corresponding to whole integration steps;
the implementation converts requested intervals to integer step counts.

### Recovery

Rerun the original commands with the same run ID and compatible inputs, software
and GPU. Completed equilibration stages transfer positions, velocities, box,
time and step count. Gaps in completed stage records stop recovery for inspection.

Production publishes complete trajectory/log chunks before advancing `prod.chk`.
An interrupted chunk replays from the previous committed checkpoint. Targets are
absolute: extending 0.1 to 1 ns requests 1 ns total. New replicas start from the
portable NPT state with fresh velocities and reset their production clock to zero.
Do not use a sibling replica's binary checkpoint. Checkpoint portability across
hardware/software versions is not assumed.

## 6. Trajectories and analysis

```bash
cd "$PROJECT_ROOT/data/$STATE/simulation/$COLABMDA_RUN_ID/$VARIANT/03_md"
cmda merge --state "$STATE" --variant "$VARIANT" --replicas r1
cd "$PROJECT_ROOT/data/$STATE/simulation/$COLABMDA_RUN_ID/$VARIANT/04_analysis"
cmda analyze --state "$STATE" --variant "$VARIANT"
```

Merge verifies contiguous completed chunks and uses cpptraj autoimaging. It writes
`prod_keep_ligand.dcd` with `stripped.prmtop`, and `prod_ca.dcd` with `ca.prmtop`.
Raw solvated chunks require the original solvated topology. Never mix these pairs.

The high-level analysis currently uses **r1**, at **10 ps spacing**. It calculates
backbone RMSD relative to the first frame, unweighted geometric protein-heavy-atom
Rg, and aligned C-alpha RMSF with checked model identity/residue labels. Saved CSVs
and `analysis.json` precede the separate figure module, which writes 7-inch,
300-DPI PNGs. Missing trajectories cause an error. PCA and contact maps are not
part of this command. For other replicas or sampling, use explicit paths and
`--interval` with `cmda openmm analysis --help`.

For matched WT/L119R protocols, compare already calculated tables:

```bash
cd "$PROJECT_ROOT/data/$STATE/simulation/$COLABMDA_RUN_ID"
cmda openmm compare --series WT="$PWD/WT/04_analysis" \
  --series L119R="$PWD/L119R/04_analysis" --outdir "$PWD/compare"
```

Review actual time ranges and residue identities before combining replicas or
comparing states. Protein/cofactor geometry, equilibration and sufficient sampling
remain scientific checks beyond successful software execution.

### Frame Counts, Sampling & Stride

With regular sampling and no extra initial frame:

$$N_\mathrm{frames}=\left\lfloor\frac{1000\,T_\mathrm{ns}}{\Delta t_\mathrm{frame,ps}}\right\rfloor$$

For durations divisible by the frame interval, 1 ns at 10 ps gives 100 frames;
100 ns at 10 ps gives 10,000. Saved frame index 0 corresponds to the first saved
frame, not necessarily physical time zero:

$$t_k=t_\mathrm{first}+k\,\Delta t_\mathrm{frame},\qquad k=0,1,\ldots,N-1$$

The lower-level `cmda openmm merge --stride s` can retain every s-th frame. Its
current implementation retains the first frame, giving:

$$N_\mathrm{kept}=\left\lceil\frac{N_\mathrm{frames}}{s}\right\rceil,\qquad
\Delta t_\mathrm{kept}=s\,\Delta t_\mathrm{frame}$$

Thus 10,000 frames with stride 10 give 1,000 frames. If the original sampling is
10 ps, the retained spacing is 100 ps. For a 1 ns trajectory starting at 10 ps,
stride 10 keeps times 10, 110, …, 910 ps; it does not reset the first frame to 0
or guarantee that the last retained frame reaches the original endpoint.

Keep the full-rate trajectory for the standard analysis and use a separately
strided copy for lighter viewing. High-level `cmda analyze` currently assumes
10 ps sampling, and lower-level analysis reconstructs time from an interval.
Before quantitative analysis of strided data, check the actual logged first-frame
time and spacing; changing an interval alone does not restore a lost time origin.
Striding reduces saved observations, not the integration time already simulated.

### Storage Calculation

For a DCD containing three 32-bit coordinate values per atom, the coordinate-data
size is approximately:

$$S_\mathrm{bytes}\approx 12\,N_\mathrm{atoms}\,N_\mathrm{frames}$$

For 35,000 atoms and 100 ns production:

| Saving interval | Frames | Approximate coordinate data per replica |
| :--- | ---: | ---: |
| 1 ps | 100,000 | 42 GB |
| 10 ps | 10,000 | 4.2 GB |

These are decimal GB estimates before file headers, box records, checkpoints,
logs or merged/stripped copies. Three replicas at 10 ps need about 12.6 GB for
these raw coordinates alone. Multiply by the number of states and variants as
well, and use the actual solvated atom count; proteins of similar residue count
can have very different solvent boxes. Coarser saving does not guarantee that
all dynamics of interest remain resolved.

## 7. Retrieve compact results

For Colab, keep trajectories and checkpoints on Drive. For HPC, use `hpc download`
for compact reports and preserve active trajectories in scratch according to the
cluster storage policy. The collection commands below are for Colab. Detach the console with Ctrl-B then D;
this leaves the runtime running. **In the local terminal**, activate the local
control environment and collect a completed variant using matching variables:

```bash
cd ~/research/tools/ColabMDA
source .venv/bin/activate
STATE=apo
VARIANT=WT
RUN_ID=cdkl5_public_v1
REMOTE_RUN="/content/drive/MyDrive/research/projects/colabmda-example/data/$STATE/simulation/$RUN_ID/$VARIANT"
LOCAL_RUN="$HOME/research/projects/colabmda-example/data/$STATE/simulation/$RUN_ID/$VARIANT"
files=(
  01_prep/prep.json 01_prep/preparation_input.json 01_prep/topology_check.json
  01_prep/PREP_DONE 02_eq/equilibration_qc.json 02_eq/nvt_protocol.json 02_eq/npt_protocol.json
  02_eq/nvt.log 02_eq/npt.log 02_eq/figures/equilibration_qc.png
  03_md/r1/prod_full.log 04_analysis/analysis.json
  04_analysis/rmsd.csv 04_analysis/rg.csv 04_analysis/rmsf.csv
  04_analysis/figures/rmsd_vs_time.png 04_analysis/figures/rg_vs_time.png
  04_analysis/figures/rmsf_per_residue.png
)
cmda colab exec -s colabmda-example --timeout 120 \
  "cd '$REMOTE_RUN' && sha256sum ${files[*]} > reports.sha256"
for file in "${files[@]}" reports.sha256; do
  mkdir -p "$LOCAL_RUN/$(dirname "$file")"
  colab download -s colabmda-example "$REMOTE_RUN/$file" "$LOCAL_RUN/$file"
done
(cd "$LOCAL_RUN" && sha256sum -c reports.sha256)
```

Require `OK` for every file. Retain the model/mapping records, source snapshot and
all stage records on Drive. Review the downloaded figures and QC, then stop the
dedicated session:

```bash
colab stop -s colabmda-example
```

Installation, CPU tests and topology checks are distinct from real GPU validation.
See the working project's run record for completed tests and current limitations.
Scientific references and software acknowledgements remain in the
[main README](../README.md#acknowledgements).
