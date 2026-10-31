# KRAS (4LDJ) Apo State: Google Colab CLI Walkthrough

This tutorial provides the complete simulation workflow for the **Apo state** of the **KRAS GTPase** (UniProt **P01116**, PDB **4LDJ** Chain A, residues **1–169**), comparing the Wild-Type (WT) protein against the oncogenic point mutant **G12D**.

Everything is executed on a **Google Colab cloud GPU via the `colab` CLI**, with all coordinate files, checkpoints, and analysis plots continuously persisted directly inside your **Google Drive**.

---

## 📁 Project Directory Layout (Google Drive)

All files reside inside your mounted Google Drive under a standardized hierarchy:

```text
/content/drive/MyDrive/research/projects/kras-example/
└── data/
    └── apo/
        ├── structure/
        │   └── P01116_4LDJ/              # Canonical UniProt P01116 + Template PDB 4LDJ
        │       ├── raw/                  # 4ldj.cif, template.pdb, uniprot.fasta, residues.json
        │       ├── variants/             # Modeled coordinates (WT/WT.pdb, G12D/G12D.pdb)
        │       └── STRUCTURE_DONE        # Verification QC gate marker
        │
        ├── simulation/
        │   ├── WT/
        │   │   ├── 01_prep/              # system.xml, model.pdb, solvated.pdb, system.prmtop, system.inpcrd
        │   │   ├── 02_eq/                # eq.pdb, system.xml, state.xml, stages/ (24-stage ladder)
        │   │   └── 03_md/                # Production replicas (100 ns)
        │   │       ├── r1/               # md.dcd, md.log, state.xml, prod.chk, chunks/
        │   │       ├── r2/               # Independent velocity seed
        │   │       └── r3/               # Triplicate simulation
        │   └── G12D/
        │       ├── 01_prep/
        │       ├── 02_eq/
        │       └── 03_md/{r1, r2, r3}/
        │
        └── analysis/
            ├── WT/                       # 02_eq/ (qc.png, qc.txt), 03_md/ (r1, r2, r3)
            ├── G12D/
            └── compare/                  # Flat comparison: rmsd.png, rg.png, rmsf.png, summary.txt
```

---

## ⚡ 1. Connect via Google Colab CLI (Local Terminal)

On your local computer (Linux, macOS, or Windows WSL), launch a Colab GPU instance and enter the console:

```bash
# 1. Launch a Colab GPU instance:
colab new -s cmda --gpu T4

# 2. Mount Google Drive (visit the printed authorization link once, then press Enter):
colab drivemount --session cmda

# 3. Enter the Colab Console:
colab console -s cmda
```

> [!NOTE]
> **No Local OpenMM or Modeller Needed:** Your local workstation only needs `google-colab-cli` (`pip install google-colab-cli`). All heavy biophysical dependencies and GPU compilation run inside the Google Colab instance.

---

## 🛠️ 2. Cloud Environment Setup (~90 Seconds)

Inside the active **Colab console terminal**, configure your Modeller license, provision the environment, and navigate to Google Drive:

```bash
# 1. Set your Modeller key & install cmda + OpenMM + AmberTools + Modeller:
export KEY_MODELLER="YOUR_MODELLER_KEY"
curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/dev1/scripts/setup.sh | bash

# 2. Activate environment:
source /root/miniforge3/etc/profile.d/conda.sh
conda activate cmda

# 3. Navigate to your persistent project folder on Google Drive:
mkdir -p /content/drive/MyDrive/research/projects/kras-example
cd /content/drive/MyDrive/research/projects/kras-example
```

> [!IMPORTANT]
> **Why `KEY_MODELLER` Is Required:**
> On Google Colab, `setup.sh` skips Modeller by default to allow unattended installation. Providing `export KEY_MODELLER="YOUR_KEY"` automatically activates the `salilab` channel, installs Modeller, and injects your license key.  
> 🔑 **Get a Free Academic License:** Register instantly at [salilab.org/modeller/registration.html](https://salilab.org/modeller/registration.html).

---

## 📋 3. Pipeline at a Glance (Cheat-Sheet)

Run these commands in order directly inside the Colab console:

```bash
# 1. Fetch & Canonical Residue Mapping
cmda fetch --state apo --pdb-id 4LDJ --uniprot-id P01116 --chain A --range 1 169

# 2. Homology Modeling & Mutations (Modeller)
cmda model --state apo --variant WT
cmda model --state apo --variant G12D

# 3. 3D pKa Titration & Solvation (PROPKA at pH 7.4 + OPC water)
cmda prep --state apo --variant WT
cmda prep --state apo --variant G12D

# 4. Stage Hubs & 24-Stage Equilibration Ladder
cmda stage --state apo --variant WT
cmda stage --state apo --variant G12D
cmda sim eq --state apo --variant WT
cmda sim eq --state apo --variant G12D

# 5. Production MD across Triplicates (100 ns)
cmda sim md --state apo --variant WT --replica r1 --total-ns 100
cmda sim md --state apo --variant G12D --replica r1 --total-ns 100

# 6. Merge Trajectory Chunks & Autoimage PBC
cmda merge --state apo --variant WT
cmda merge --state apo --variant G12D

# 7. Audit Status & Frame Counts
cmda status --state apo --variant WT --replica r1
cmda status --state apo --variant G12D --replica r1

# 8. Single-Variant Analysis (RMSD, Rg, RMSF)
cmda analyze --state apo --variant WT
cmda analyze --state apo --variant G12D

# 9. Direct Comparative Overlay (WT vs G12D with Mean ± SD band)
cmda compare --state apo --variants WT,G12D
```

---

## 🔬 4. Step-by-Step Scientific Walkthrough

### 4.1. Structural Fetching & Biological UniProt Mapping
Downloads crystallographic coordinates, detects unresolved loops, and aligns author residues to canonical UniProt indexing:

```bash
cmda fetch --state apo --pdb-id 4LDJ --uniprot-id P01116 --chain A --range 1 169
```

*What `cmda fetch` produces in `data/apo/structure/P01116_4LDJ/raw/`:*
* `4ldj.cif`: Raw crystallographic record.
* `template.pdb`: Stripped monomer coordinates (Chain A).
* `uniprot.fasta`: Canonical sequence fetched directly from UniProtKB.
* `residues.json`: Complete author-to-canonical residue index map.
* `STRUCTURE_DONE`: Verification gate marker.

---

### 4.2. Homology Modeling & Point Mutations (Modeller)
Closes missing loop coordinates on the canonical backbone, then introduces the G12D oncogenic point mutation:

```bash
# 1. Build canonical Wild-Type model (residues 1–169):
cmda model --state apo --variant WT

# 2. Build Point Mutant (G12D) on the verified WT backbone:
cmda model --state apo --variant G12D
```

---

### 4.3. Structure Preparation & 3D pKa Titration (PROPKA at pH 7.4)
Evaluates 3D electrostatic microenvironments using **PROPKA 3**, assigns explicit AMBER protonation states (`HID`/`HIE`/`HIP`, `ASH`/`GLH`, `LYN`), solvates with OPC 4-point water, and neutralizes with counterions:

```bash
cmda prep --state apo --variant WT
cmda prep --state apo --variant G12D
```

*Key outputs in `data/apo/simulation/<variant>/01_prep/`:*
* `solvated.pdb`: Solvated coordinate PDB with protein, OPC water box, and neutralized counterions.
* `system.prmtop`: Amber topology containing all `ff19SB` physical parameters.
* `system.inpcrd`: Solvated, charge-neutralized initial coordinates.

---

### 4.4. Stage Simulation Hubs & 24-Stage Equilibration Ladder
Initializes simulation hubs and executes the progressive restraint release ladder (EM $\to$ NVT heating to 300 K $\to$ NPT density relaxation $\to$ QC check):

```bash
# 1. Initialize simulation directories:
cmda stage --state apo --variant WT
cmda stage --state apo --variant G12D

# 2. Run automated equilibration (EM -> NVT -> NPT -> QC):
cmda sim eq --state apo --variant WT
cmda sim eq --state apo --variant G12D
```

*Verify Equilibration Stability in 1 Second:*
```bash
cat data/apo/analysis/WT/02_eq/qc.txt
cat data/apo/analysis/G12D/02_eq/qc.txt
```
*Expected Output:* `Status: PASSED`, `Cα RMSD: ~0.83 Å` (Confirms temperature $\approx 300\text{ K}$, density $\approx 1.01\text{ g/mL}$, fold intact).

---

### 4.5. Production MD across Triplicates (`03_md/{r1, r2, r3}`)
Replicas fork from `02_eq/` with distinct random velocity seeds:

```bash
# Wild-Type Triplicates (100 ns production runs):
cmda sim md --state apo --variant WT --replica r1 --total-ns 100
cmda sim md --state apo --variant WT --replica r2 --total-ns 100
cmda sim md --state apo --variant WT --replica r3 --total-ns 100

# Mutant Triplicates (100 ns production runs):
cmda sim md --state apo --variant G12D --replica r1 --total-ns 100
cmda sim md --state apo --variant G12D --replica r2 --total-ns 100
cmda sim md --state apo --variant G12D --replica r3 --total-ns 100
```

---

### 4.6. Merge Trajectories & Autoimage PBC
Concatenates 1-ns chunk files into a continuous trajectory, applies periodic boundary condition (PBC) autoimaging, and strips bulk solvent:

```bash
# Merge Wild-Type replicas:
cmda merge --state apo --variant WT

# Merge Mutant replicas:
cmda merge --state apo --variant G12D
```

*What `cmda merge` outputs in `data/apo/simulation/<variant>/03_md/<replica>/`:*
* `md.dcd`: Continuous merged trajectory.
* `md.log`: Continuous merged thermodynamic log table.
* `state.xml`: Final OpenMM restart state.

---

### 4.7. Trajectory Integrity Auditing (`cmda status`)
Verify chunk continuity, frame counts, and resume readiness:

```bash
cmda status --state apo --variant WT --replica r1
```

*Example status report:*
```text
[STATUS]
  Workdir          : data/apo/simulation/WT/03_md/r1
  Chunks (DCD/log) : 100 / 100
  Topology File    : solvated.pdb
                     └─ 27273 atoms, 8388 residues
                     └─ (169 protein, 8170 water, 49 ions)
  Trajectory File  : md.dcd (10000 frames)
  Log File         : md.log (10000 frames)
  Frames (from logs): 10000
  Sim time (ps/ns) : 100000.00 ps / 100.0000 ns
  Resume-ready     : YES (needs prod.chk, system.xml, solvated.pdb)
```

---

### 4.8. Single-Variant Biophysical Analysis
Calculate backbone RMSD, Radius of Gyration, and per-residue RMSF for each variant:

```bash
# Analyze Wild-Type and Mutant trajectories:
cmda analyze --state apo --variant WT
cmda analyze --state apo --variant G12D
```

Outputs saved in `data/apo/analysis/<variant>/03_md/<replica>/`:
* `rmsd.png`: Single-system backbone RMSD over 100 ns.
* `rg.png`: Single-system Radius of Gyration compactness.
* `rmsf.png`: Single-system per-residue fluctuation profile.

---

### 4.9. Direct Multi-Variant Comparison
Generate direct comparative overlay figures (WT vs. G12D with Mean ± SD bands):

```bash
cmda compare --state apo --variants WT,G12D
```

*Key deliverables generated in `data/apo/analysis/compare/`:*
* `summary.txt`: Executive statistical summary (Mean ± SD for RMSD, Rg, Temp).
* `rmsd.png`: Comparative RMSD time-series with standard deviation shading.
* `rg.png`: Comparative Radius of Gyration compactness.
* `rmsf.png`: Per-residue fluctuation profiles highlighting Switch I (residues 30–38) and Switch II (residues 59–67) loop dynamics.

---

## 📥 5. Pull Deliverables to Local Machine (Local Terminal)

Exit the Colab console (`exit`), then in your local workstation terminal, pull the generated figures, summary tables, and reports to your computer:

```bash
mkdir -p ~/research/projects/kras-example
cd ~/research/projects/kras-example

# Pull lightweight deliverables (~2 MB):
cmda pull -s cmda

# To also pull post-processed stripped trajectories for PyMOL inspection:
cmda pull -s cmda --traj
```

---

## 💡 Practical Reference & Troubleshooting

### How to Resume After a Session Timeout
Production MD continuously saves periodic checkpoints (`prod.chk`) directly to Google Drive every 100 ps. If your Colab session expires:

1. Re-launch and enter console from your local terminal:
   ```bash
   colab new -s cmda --gpu T4
   colab drivemount --session cmda
   colab console -s cmda
   ```
2. Re-activate the environment:
   ```bash
   source /root/miniforge3/etc/profile.d/conda.sh
   conda activate cmda
   cd /content/drive/MyDrive/research/projects/kras-example
   ```
3. Re-run the exact same `cmda sim md` command:
   ```bash
   cmda sim md --state apo --variant WT --replica r1 --total-ns 100
   ```
   *ColabMDA automatically detects `prod.chk` on Google Drive and resumes dynamics seamlessly.*

---

### In-Context Reference: Simulation Math & Storage Budgeting
Calculate expected trajectory frames:
$$\text{Expected Frames} = \frac{\text{Total Time (ps)}}{\text{Trajectory Saving Interval (ps)}}$$

* A $100.0\text{ ns}$ ($100,000\text{ ps}$) simulation with the default $10.0\text{ ps}$ saving interval:
  $$\text{Expected Frames} = \frac{100,000\text{ ps}}{10\text{ ps}} = 10,000\text{ frames}$$
* **Default 10 ps Interval (`--traj-interval 10`):** In `dev1`, `--traj-interval` defaults to 10 ps, reducing total raw storage from **~36 GB to ~3.6 GB (a 90% space reduction)** so it comfortably fits within Google Drive limits.
* **Stripped Trajectories (`cmda merge`):** Removing solvent water retains only protein coordinates (~2,600 atoms), reducing merged trajectory file sizes to **<200 MB** per replica.

---

### In-Context Reference: Topology Pairing Guidelines

| Trajectory File | Description | Matching Topology File | Atom Count | Intended Usage |
| :--- | :--- | :--- | :--- | :--- |
| `prod_raw.dcd` | Raw chunked coordinates | `01_prep/solvated.pdb` or `system.prmtop` | ~27,000 atoms | Raw simulation checkpointing & resume. |
| `prod_keep_ligand.dcd` | Merged, PBC autoimaged, solvent-stripped | `03_md/r1/stripped.prmtop` | ~2,600 atoms | Standard biophysical analysis (`cmda analyze`, `cmda compare`). |
| `prod_ca.dcd` | Merged Cα backbone only | `03_md/r1/ca.prmtop` | ~169 atoms | Rapid PCA and essential dynamics. |

> [!WARNING]
> **Topology Pairing Requirement:**
> Loading `prod_keep_ligand.dcd` (~2,600 atoms) with `solvated.pdb` (~27,000 atoms) in PyMOL will cause an **atom mismatch error**. Always pair stripped trajectories with `stripped.prmtop`.

---

## Acknowledgements
- [OpenMM](https://openmm.org) & [PDBFixer](https://github.com/openmm/pdbfixer)
- [Modeller](https://salilab.org/modeller/)
- [MDAnalysis](https://www.mdanalysis.org) & [MDTraj](https://www.mdtraj.org)
- [PDB2PQR & PROPKA](https://server.poissonboltzmann.org/)
- [NumPy](https://numpy.org), [Matplotlib](https://matplotlib.org), [Biopython](https://biopython.org)
- [Google Colab](https://colab.research.google.com) & [Miniforge/Conda](https://github.com/conda-forge/miniforge)

---

## Citation

Initially, this repository was based on a GROMACS-on-Colab workflow, but was transitioned to OpenMM to improve performance and stability for automated simulation recovery on cloud platforms.

If you use this tool, please cite the underlying studies:

> Paul SK, Saddam M, Rahaman KA, Choi JG, Lee SS, Hasan M. **Molecular modeling, molecular dynamics simulation, and essential dynamics analysis of grancalcin: An upregulated biomarker in experimental autoimmune encephalomyelitis mice.** *Heliyon*. 2022 Oct 23;8(10):e11232. doi: [10.1016/j.heliyon.2022.e11232](https://doi.org/10.1016/j.heliyon.2022.e11232). PMID: 36340004; PMCID: PMC9626934.
>
> [![Cited by](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fapi.openalex.org%2Fworks%2Fhttps%3A%2F%2Fdoi.org%2F10.1016%2Fj.heliyon.2022.e11232&query=%24.cited_by_count&label=cited%20by&color=blue)](https://doi.org/10.1016/j.heliyon.2022.e11232)
>
> Paul SK, Saddam M, Tabassum N, Hasan M. **Molecular dynamics simulation of wild and mutant proteasome subunit beta type 8 (PSMB8) protein: Implications for restoration of inflammation in experimental autoimmune encephalomyelitis pathogenesis.** *Heliyon*. 2024 Dec 15;11(1):e41166. doi: [10.1016/j.heliyon.2024.e41166](https://doi.org/10.1016/j.heliyon.2024.e41166). PMID: 39802026; PMCID: PMC11719297.
>
> [![Cited by](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fapi.openalex.org%2Fworks%2Fhttps%3A%2F%2Fdoi.org%2F10.1016%2Fj.heliyon.2024.e41166&query=%24.cited_by_count&label=cited%20by&color=blue)](https://doi.org/10.1016/j.heliyon.2024.e41166)

---

<p align="center">
  <img src="https://komarev.com/ghpvc/?username=paulshamrat&repo=ColabMDA&label=Visitors&color=blue&style=flat-square" alt="Visitors" />
</p>
