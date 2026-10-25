# Optional examples and states

The main README follows 4BGQ CDKL5 apo and ATP/Mg. These older templates and
advanced states remain available separately. FKBP12 and H-Ras can be installed
with `cmda stage --example fkbp` and `cmda stage --example hras`; each has its own
structure/model/preparation checks. They are not steps in the CDKL5 workflow.

## Trp-cage quick test

### 1. Create and Fetch Trp-cage

The 20-residue Trp-cage example tests software behavior quickly without replacing
the full research workflows in the main README. **On Colab after the [main setup](../README.md#11-google-colab)**, run:

```bash
cmda stage --example trpcage \
  --root /content/drive/MyDrive/research/projects/colabmda-example
cd /content/drive/MyDrive/research/projects/colabmda-example
python scripts/01_structure/01_fetch.py
```

On a local GPU workstation use `~/research/projects/colabmda-example` instead.
Fetch prefers mmCIF, with PDB fallback, selecting NMR model 1, author chain A of
1L2Y. It verifies the deposited sequence, all 20 C-alpha atoms and native numbering,
then records the archive/selection hashes. As a designed peptide, it has no canonical
UniProt mutation annotation. No MODELLER key or cofactor template is needed.

Staging copies seven scripts and records their tool version and hashes. It does
not fetch inputs or run a simulation. Repeating staging preserves identical files
and stops if existing files differ.

```text
colabmda-example/
├── scripts/
│   ├── 01_structure/01_fetch.py
│   └── 02_simulation/
│       ├── common.sh
│       ├── 01_prepare.sh
│       ├── 02_equilibrate.sh
│       ├── 03_production.sh
│       ├── 04_analysis.sh
│       └── 05_check.py
├── data/
│   ├── 00_setup/example.json
│   ├── 01_structure/trpcage/           # Raw archive, selected model, input hashes
│   └── 02_simulation/<run>/
│       ├── provenance.json
│       ├── production.json
│       └── apo/WT/
│           ├── 01_prep/
│           ├── 02_eq/
│           ├── 03_md/r1/
│           └── 04_analysis/figures/
└── logs/02_simulation/<run>/
```

### 2. Run the Scripts in Order

With `colabmda` active in the generated project:

```bash
export COLABMDA_RUN_ID=run01
bash scripts/02_simulation/01_prepare.sh
bash scripts/02_simulation/02_equilibrate.sh
bash scripts/02_simulation/03_production.sh
bash scripts/02_simulation/04_analysis.sh
python scripts/02_simulation/05_check.py
```

The defaults are 100 ps NVT, 100 ps NPT, 1 ns production, 2 fs steps, 10 ps sampling
and 100 ps checkpoint chunks. To choose another duration without editing scripts,
set `COLABMDA_EXAMPLE_NS` for production; the final check reads the saved request.
Use positive multiples of 0.1 ns. Choose a fresh run ID for changed inputs,
executed scripts or runtime revisions.

To test extension instead of the direct production command above, use a fresh run
ID, complete preparation/equilibration, then run:

```bash
COLABMDA_EXAMPLE_NS=0.1 bash scripts/02_simulation/03_production.sh
COLABMDA_EXAMPLE_NS=1 bash scripts/02_simulation/03_production.sh
```

Then perform analysis and the final check. To include an actual interruption in
this test, pause after the 0.1 ns smoke command, save its first chunk's checksum
from the project root, then start the extension:

```bash
sha256sum "data/02_simulation/$COLABMDA_RUN_ID/apo/WT/03_md/r1/prod_0to100ps.dcd" \
  > "data/02_simulation/$COLABMDA_RUN_ID/first_chunk.sha256"
COLABMDA_EXAMPLE_NS=1 bash scripts/02_simulation/03_production.sh
```

Press Ctrl+C while the second chunk is running. Before another command, `echo $?`
should report a nonzero exit status. Resume with the same production command, then
verify the first committed chunk:

```bash
COLABMDA_EXAMPLE_NS=1 bash scripts/02_simulation/03_production.sh
sha256sum -c "data/02_simulation/$COLABMDA_RUN_ID/first_chunk.sha256"
```

The checksum must report `OK`. Continue with analysis and the final check.
The private working project's `tests/gpu.py` automates that acceptance test using
these same generated scripts; it is not an extra file installed for every user.

### 3. Expected Outputs & Checks

After success, `data/02_simulation/run01/` contains the input/tool/script/runtime
provenance, production requests, analysis results, `summary.json` and `RUN_DONE`.
The default run has ten chunks and 100 merged frames, with times 10–1000 ps.
The final check requires finite coordinates/logs/tables, 20 C-alpha atoms, native
RMSF labels 1–20 and nonempty PNGs. Physical equilibration QC is reported separately.

Scripts create the log directory; save terminal output when desired. For example,
from the generated project after preparation/equilibration:

```bash
set -o pipefail
bash scripts/02_simulation/03_production.sh 2>&1 \
  | tee "logs/02_simulation/$COLABMDA_RUN_ID/production.log"
```

To rerender analysis figures from saved tables without rerunning MD:

```bash
python - <<'PY'
import os
from pathlib import Path
from colabmda.openmm.figures import plot_analysis
run = Path('data/02_simulation') / os.environ['COLABMDA_RUN_ID']
plot_analysis(run / 'apo/WT/04_analysis')
PY
```

### 4. Runtime, Limits & Cleanup

Our earlier T4 test took approximately 7.5 minutes for production plus its restart
test. That is an observed runtime, not a guarantee. Quick equilibration can fail
the density-slope criterion; retain the warning and do not call this a converged
production study, a mutant/cofactor validation or a full staged-protocol test.

Keep heavy trajectories and checkpoints on Drive. To return from the interactive
console to your local terminal, detach with **Ctrl-B, then D**. This does not stop
the runtime. With the local ColabMDA environment active, download and verify the
compact reports below. Set RUN_ID, STATE and VARIANT to the completed system.
The commands below retrieve the Trp-cage run. CDKL5 uses the state-based paths
in the [apo](apo.md#7-retrieve-compact-results) and
[ATP/Mg](atpmg.md#7-retrieve-compact-results) collection instructions:

```bash
RUN_ID=run01
STATE=apo
VARIANT=WT
SYSTEM="$STATE/$VARIANT"
REMOTE_RUN="/content/drive/MyDrive/research/projects/colabmda-example/data/02_simulation/$RUN_ID"
LOCAL_RUN="$HOME/research/projects/colabmda-example/data/02_simulation/$RUN_ID"
files=(
  summary.json provenance.json production.json RUN_DONE
  "$SYSTEM/01_prep/prep.json"
  "$SYSTEM/02_eq/equilibration_qc.json" "$SYSTEM/02_eq/figures/equilibration_qc.png"
  "$SYSTEM/02_eq/nvt.log" "$SYSTEM/02_eq/npt.log" "$SYSTEM/03_md/r1/prod_full.log"
  "$SYSTEM/04_analysis/analysis.json"
  "$SYSTEM/04_analysis/rmsd.csv" "$SYSTEM/04_analysis/rg.csv" "$SYSTEM/04_analysis/rmsf.csv"
  "$SYSTEM/04_analysis/figures/rmsd_vs_time.png"
  "$SYSTEM/04_analysis/figures/rg_vs_time.png"
  "$SYSTEM/04_analysis/figures/rmsf_per_residue.png"
)
cmda colab exec -s colabmda-example --timeout 120 \
  "cd '$REMOTE_RUN' && sha256sum ${files[*]} > reports.sha256"
for file in "${files[@]}" reports.sha256; do
  mkdir -p "$LOCAL_RUN/$(dirname "$file")"
  colab download -s colabmda-example "$REMOTE_RUN/$file" "$LOCAL_RUN/$file"
done
(cd "$LOCAL_RUN" && sha256sum -c reports.sha256)
```

Require `OK` for every downloaded report, review the PNGs and read the physical
QC result. Then stop the dedicated session **from the local terminal**:

```bash
colab stop -s colabmda-example
```

A restarted session requires Drive authentication and the recorded software again.
Historical run records remain in the working project's README; user instructions
stay here. CPU unit tests, script/CLI checks and a successful installation are
reported separately from an actual GPU integration run.

---


## Additional chemical states

ATP, ADP and GTP resources remain available for appropriate structures. Choose
state tags and ion counts for the actual experimental/model composition. Existing
phosphorylation and substrate-peptide conventions remain advanced capabilities;
they are not required by the small tutorials or validated by their acceptance tests.
For kinase projects, the previously documented composition tags remain available:

| State tag | Requested composition |
| :--- | :--- |
| `adp_2mg_py171` | ADP, two Mg ions and phosphorylated Tyr171. |
| `atp_2mg_py171_peb2` | ATP, two Mg ions, phosphorylated Tyr171 and the EB2 peptide preset. |

These tags request components; review phosphorylation parameters, peptide
placement and coordination before interpreting an activation state.
Generic GAFF parameters are not equivalent to the specialized nucleotide route.
Preserve a matching Amber topology/coordinate pair for independently parameterized
systems, and review protonation, charge, ligand atom names and coordination.

