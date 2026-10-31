# KRAS (4LDJ) Apo State: Google Colab CLI Walkthrough

This tutorial provides a streamlined simulation walkthrough for the **Apo state** of the **KRAS GTPase** (UniProt **P01116**, PDB **4LDJ** Chain A, residues **1–169**), comparing the Wild-Type (WT) protein against the oncogenic point mutant **G12D**.

Everything is executed on a **Google Colab cloud GPU via the `colab` CLI**, with all outputs, checkpoints, and figures continuously persisted inside your **Google Drive**.

---

## 📁 Directory Architecture (Google Drive)

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

## 1. Connect via Google Colab CLI (Local Terminal)

On your local workstation (Linux, macOS, or Windows WSL), launch a Colab GPU instance and enter the console:

```bash
# 1. Launch a Colab GPU instance:
colab new -s cmda --gpu T4

# 2. Mount Google Drive (visit the printed link once, then press Enter):
colab drivemount --session cmda

# 3. Enter the Colab Console:
colab console -s cmda
```

> [!TIP]
> **How to Resume if Session Expires:**
> If your Colab session times out during a run, simply repeat the 3 commands above to re-enter the console. When you restart your simulation command, ColabMDA detects `prod.chk` on Google Drive and **resumes dynamics from the exact last saved timestep**.

---

## 2. Cloud Environment Setup (~90 Seconds)

Inside the active **Colab console terminal**, provide your Modeller license key, run the bootstrap installer, and navigate to Google Drive:

```bash
# 1. Export Modeller key & run automated installer:
export KEY_MODELLER="YOUR_MODELLER_KEY"
curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/dev1/scripts/setup.sh | bash

# 2. Activate environment:
source /root/miniforge3/etc/profile.d/conda.sh
conda activate cmda

# 3. Navigate to persistent project folder on Google Drive:
mkdir -p /content/drive/MyDrive/research/projects/kras-example
cd /content/drive/MyDrive/research/projects/kras-example
```

> [!NOTE]
> Setting `KEY_MODELLER` tells `setup.sh` to install Modeller from the official `salilab` channel and inject your license. Academic keys are free at [salilab.org/modeller/registration.html](https://salilab.org/modeller/registration.html).

---

## 3. Homology Modeling & Point Mutations (Modeller)

Fetch experimental coordinates from the PDB archive, align to canonical UniProt indexing, and model the full WT domain and G12D point mutation:

```bash
# 1. Download 4LDJ (Chain A, residues 1–169) & verify canonical mapping:
cmda fetch --state apo --pdb-id 4LDJ --uniprot-id P01116 --chain A --range 1 169

# 2. Model canonical Wild-Type coordinates (closes missing crystallographic loops):
cmda model --state apo --variant WT

# 3. Build Point Mutant (G12D) on the verified WT backbone:
cmda model --state apo --variant G12D
```

---

## 4. Structure Preparation & 3D pKa Titration (PROPKA at pH 7.4)

ColabMDA runs **PROPKA 3** to calculate 3D electrostatic microenvironments, assign explicit AMBER titratable protonation states, solvate with OPC 4-point water, and neutralize net charge:

```bash
cmda prep --state apo --variant WT
cmda prep --state apo --variant G12D
```

| Residue | State | Physical Criteria |
| :--- | :--- | :--- |
| **Histidine** | `HID` / `HIE` / `HIP` | Neutral $\text{N}_\delta$, neutral $\text{N}_\epsilon$ (default), or positively charged doubly protonated ($pK_a > 7.4$). |
| **Acidic** | `ASP` / `GLU` | Standard carboxylate anions (`ASH`/`GLH` assigned if buried with $pK_a > 7.4$). |
| **Basic** | `LYS` / `ARG` | Positively charged physiological cations (`LYN` if neutral). |
| **Cysteine** | `CYS` / `CYX` | Free neutral thiol (`CYS`) or disulfide bridge (`CYX`). |

*Outputs in `data/apo/simulation/<variant>/01_prep/`:* `solvated.pdb`, `system.prmtop`, `system.inpcrd`, and `system.xml`.

---

## 5. 24-Stage Equilibration Ladder & QC Gating

Initialize simulation directories and run the progressive 24-stage VarMDyn relaxation ladder. This heats the system from 50 K to 300 K and gradually releases positional restraints:

```bash
# 1. Initialize simulation directories:
cmda stage --state apo --variant WT
cmda stage --state apo --variant G12D

# 2. Execute 24-stage progressive equilibration (EM -> NVT -> NPT -> QC):
cmda sim eq --state apo --variant WT
cmda sim eq --state apo --variant G12D
```

### Exact Stage Durations (VarMDyn Protocol, $\approx 1.96\text{ ns}$ Total):

| Stages | Operation | Duration | Timestep & Restraints |
| :--- | :--- | :--- | :--- |
| **01–03** | Restrained Minimization | 5,000 steps each | $k = 100, 50, 25\text{ kcal/mol/\AA}^2$ |
| **04–09** | NVT Progressive Heating ($50\text{ K} \to 300\text{ K}$) | **$6 \times 10\text{ ps}$ ($60\text{ ps}$)** | $1\text{ fs timestep}$, $k = 25$ |
| **10–14** | Post-Heating Minimization | 5,000 / 1,000 steps | $k = 25 \to 1\text{ kcal/mol/\AA}^2$ |
| **15–18** | NVT Restraint Release at $300\text{ K}$ | **$4 \times 100\text{ ps}$ ($400\text{ ps}$)** | $2\text{ fs timestep}$, $k = 10, 5, 2, 1$ |
| **19–22** | NPT Density Relaxation ($300\text{ K}, 1\text{ bar}$) | **$4 \times 100\text{ ps}$ ($400\text{ ps}$)** | $2\text{ fs timestep}$, $k = 10, 5, 2, 1$ |
| **23** | Unrestrained NPT Relaxation | **$100\text{ ps}$** | $2\text{ fs timestep}$, unrestrained ($k = 0$) |
| **24** | Final Unrestrained NPT Stabilization | **$1,000\text{ ps}$ ($1.0\text{ ns}$)** | $2\text{ fs timestep}$, unrestrained ($k = 0$) |

> [!TIP]
> **Fast Smoke Test Alternative:** For a 2-minute quick check rather than full 24-stage production relaxation, append `--protocol quick` (runs 100 ps NVT + 100 ps NPT).

#### 1-Second QC Verification:
```bash
cat data/apo/analysis/WT/02_eq/qc.txt
cat data/apo/analysis/G12D/02_eq/qc.txt
```
*Expected Output:* `Status: PASSED`, `Cα RMSD: ~0.83 Å` (Confirms temperature $\approx 300\text{ K}$, density $\approx 1.01\text{ g/mL}$).

---

## 6. Production MD Across Triplicates (100 ns)

Fork production replicas directly from `02_eq/` with distinct random velocity seeds:

```bash
# Wild-Type Triplicates (100 ns each):
cmda sim md --state apo --variant WT --replica r1 --total-ns 100
cmda sim md --state apo --variant WT --replica r2 --total-ns 100
cmda sim md --state apo --variant WT --replica r3 --total-ns 100

# Mutant Triplicates (100 ns each):
cmda sim md --state apo --variant G12D --replica r1 --total-ns 100
cmda sim md --state apo --variant G12D --replica r2 --total-ns 100
cmda sim md --state apo --variant G12D --replica r3 --total-ns 100
```

> [!TIP]
> **Storage Budgeting & Trajectory Interval:**
> By default, `cmda sim md` saves snapshots every **$10\text{ ps}$** (`--traj-interval 10`).  
> $$\text{Expected Frames} = \frac{100,000\text{ ps}}{10\text{ ps}} = 10,000\text{ frames}$$
> This provides high temporal resolution while reducing raw storage from **~36 GB to ~3.6 GB (a 90% space reduction)**, comfortably fitting within Google Drive limits.

---

## 7. Merge, Status Audit & Comparative Analysis

Concatenate raw 1-ns chunks, apply periodic boundary condition (PBC) autoimaging, strip bulk water, audit frame continuity, and generate comparative figures:

```bash
# 1. Merge chunks & autoimage PBC (produces continuous md.dcd and stripped trajectories):
cmda merge --state apo --variant WT
cmda merge --state apo --variant G12D

# 2. Audit frame counts & chunk continuity:
cmda status --state apo --variant WT --replica r1
cmda status --state apo --variant G12D --replica r1

# 3. Calculate single-variant biophysical metrics (RMSD, Rg, RMSF):
cmda analyze --state apo --variant WT
cmda analyze --state apo --variant G12D

# 4. Generate direct comparative overlay figures (WT vs G12D with Mean ± SD bands):
cmda compare --state apo --variants WT,G12D
```

> [!WARNING]
> **Topology Pairing Requirement:**
> Trajectory stripping reduces file size to **<200 MB** by removing bulk solvent (~2,600 protein atoms remaining). When viewing `prod_keep_ligand.dcd` in PyMOL or MDAnalysis, **always pair it with `stripped.prmtop`**. Loading with `solvated.pdb` (~27,000 atoms) will trigger a fatal atom-mismatch error.

*Key Deliverables Generated in `data/apo/analysis/compare/`:*
* `summary.txt`: Executive statistical table with Mean ± SD for RMSD, Radius of Gyration, and Temperature.
* `rmsd.png`: Comparative backbone RMSD time-series.
* `rg.png`: Comparative compactness profiles.
* `rmsf.png`: Per-residue fluctuation profiles highlighting Switch I (residues 30–38) and Switch II (residues 59–67) loop flexibility.

---

## 8. Pull Deliverables to Local Machine (Local Terminal)

Exit the Colab console (`exit`), then in your local workstation terminal, download the analyzed figures, reports, and stripped trajectories:

```bash
mkdir -p ~/research/projects/kras-example
cd ~/research/projects/kras-example

# Pull lightweight figures, logs, and summary tables (~2 MB):
cmda pull -s cmda

# To also pull solvent-stripped trajectories for PyMOL inspection:
cmda pull -s cmda --traj
```

View the generated publication comparison in your favorite image viewer:
```bash
# View comparative overlay plot locally:
xdg-open data/apo/analysis/compare/rmsd.png  # Linux
# or open in PyMOL:
pymol data/apo/simulation/WT/03_md/r1/stripped.prmtop data/apo/simulation/WT/03_md/r1/prod_keep_ligand.dcd
```

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
