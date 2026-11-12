# CDKL5 Holo State (ATP + 1 Mg²⁺): Complete Step-by-Step Walkthrough

This tutorial provides the complete, end-to-end simulation workflow for **Track 2: Holo Dynamics**, simulating the catalytically competent active-site complex of human cyclin-dependent kinase-like 5 (**CDKL5**, UniProt **O76039**, PDB **4BGQ**, residues **1–303**) bound to adenosine triphosphate (**ATP**) and a coordinating divalent magnesium ion ($\text{Mg}^{2+}$).

To assemble this complex with crystallographic accuracy, ColabMDA transfers cofactor coordinates from a high-resolution homologous kinase donor structure (**CDK2**, PDB **8FP5** at 1.45 Å resolution) onto your modeled CDKL5 structures, parameterizes the active site using quantum-derived nucleotide force fields, and executes unrestrained molecular dynamics across triplicates.

---

Outputs are shown after each stage. The complete directory tree is collapsed at
the end of this tutorial. Paths are relative to the project root unless stated otherwise.

---

## 💻 Phase 1: Local Holo Assembly & Structure Preparation (Local Machine)
*Execute this phase on your local laptop or workstation CPU. No GPU is required.*

### 1.0. Local Machine Setup (One-Time)
Run this one-time setup in your terminal on **Linux (Ubuntu)** or Windows via **[WSL2](https://learn.microsoft.com/en-us/windows/wsl/install)**:

```bash
# 1. Install required system tools (Ubuntu/WSL2):
sudo apt-get update
sudo apt-get install -y git curl wget ca-certificates

# 2. Install cmda, OpenMM, AmberTools, PyMOL, and Modeller:
bash -c "$(curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/main/scripts/setup.sh)"

# 3. Activate the environment in this terminal:
source ~/miniforge3/etc/profile.d/conda.sh
conda activate cmda

# 4. Log in to Google Colab (one-time browser login):
colab sessions

# 5. Navigate to project directory:
cd ~/research/projects/cmda-example
```

Run the activation commands in every new terminal or SSH session. The installer
runs in a child shell, so it cannot activate the environment in your current terminal.

When `colab sessions` prints an authorization URL:

1. Keep that terminal open and open the printed URL in a browser. Over SSH, use
   a browser on your own computer.
2. Sign in with the Google account you use for Colab and approve access.
3. Copy the authorization code Google displays.
4. Paste the code into the same terminal at `Enter the authorization code:` and
   press Enter. Keep the code private.

> [!TIP]
> **Modeller Key:** Free for academic research ([register here](https://salilab.org/modeller/registration.html)). Enter it when prompted, or run:
> `KEY_MODELLER="YOUR_KEY" bash -c "$(curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/main/scripts/setup.sh)"`

### 1.1. Structural Assembly: Homologous Cofactor Transfer from CDK2 (8FP5)

ColabMDA transfers cofactor coordinates from human CDK2 (PDB **8FP5**, resolved at 1.45 Å with ATP and $\text{Mg}^{2+}$ bound in a catalytically competent Michaelis pose) onto your modeled CDKL5 structures.

The tool performs a structural superposition of conserved catalytic core residues (residues 30–220, covering the glycine-rich P-loop, the catalytic Lys42/Glu60 pair, and the DFG motif), transfers the cofactors, and performs fixed-backbone local sidechain relaxation to relieve any steric clashes:

```bash
# 1. Assemble Holo Wild-Type from CDK2 donor (PDB 8FP5) with pocket sidechain relaxation
cmda assemble --state atpmg --variant WT --template 8FP5 --relax-sidechains
```

<details>
<summary><b>View stage directory tree</b></summary>

```text
data/atpmg/structure/O76039_4BGQ/
├── donors/
│   └── 8FP5/
│       └── 8fp5.cif                                 # Cofactor donor archive; PDB fallback is possible
└── variants/
    └── WT/
        ├── WT.pdb                                   # Assembled protein, ATP and Mg
        ├── complex.pdb                              # Copy of the assembled complex
        └── ASSEMBLY_DONE                            # Verified assembly and donor fit
```

</details>

<details open>
<summary><b>📄 Expected Terminal Output: Catalytic Superposition & Relaxation</b></summary>

```text
==============================================================================
 ColabMDA Holo Structural Assembly & Catalytic Alignment
==============================================================================
 Target Receptor:     CDKL5 WT (data/apo/structure/O76039_4BGQ/variants/WT/WT.pdb)
 Cofactor Donor:      CDK2 (PDB 8FP5 Chain A, Resolution: 1.45 Å)
 Core Alignment:      Residues 30–220 superposed (191 Cα pairs)
 Alignment RMSD:      0.684 Å across catalytic core
 Cofactors Transferred: ATP (Residue 401, 43 atoms, Net Charge: -4)
                       MG  (Residue 402,  1 atom,  Net Charge: +2)
 Initial Clashes:     3 sidechain steric overlaps detected (Met86, Leu137, Phe157)
 Pocket Relaxation:   Fixed-backbone local gradient minimization (500 steps)
 Final Clashes:       0 steric overlaps (min interatomic distance: 2.38 Å)
 Catalytic Gate:      ✓ Mg²⁺–O(β,γ) coordination distance: 2.04 Å (Target: 1.9–2.1 Å)
                      ✓ Lys42(Nζ)–O(α) distance: 2.82 Å (Salt bridge formed)
                      ✓ Asp156(Oδ)–Mg²⁺ distance: 2.08 Å (Octahedral sphere)
 Outputs Saved:       data/atpmg/structure/O76039_4BGQ/variants/WT/WT.pdb
                      data/atpmg/structure/O76039_4BGQ/donors/8FP5/8fp5.cif
==============================================================================
```
</details>

#### PyMOL Active-Site Inspection & Catalytic Coordination Gate

Visually inspect the assembled active site in PyMOL before moving to solvation:

```bash
pymol data/atpmg/structure/O76039_4BGQ/variants/WT/WT.pdb
```

In the PyMOL command line, highlight the catalytic pocket:
```text
select atpmg_pocket, resn ATP or resn MG or (resi 30-220 and within 5.0 of (resn ATP or resn MG))
show sticks, atpmg_pocket
zoom atpmg_pocket
```

> [!TIP]
> **Active-Site Coordination Verification Checklist:**
> 1. **Adenine Hinge Anchor:** The adenine ring must insert into the hydrophobic cleft between the N- and C-lobes (near CDKL5 hinge residues Glu81–Leu83).
> 2. **Magnesium Chelation:** Divalent $\text{Mg}^{2+}$ must coordinate oxygens of the $\beta$- and $\gamma$-phosphate groups with distances between **1.9 Å and 2.1 Å**.
> 3. **Catalytic Asp156 (DFG Motif):** Carboxylate oxygens of Asp156 must orient toward the $\text{Mg}^{2+}$ ion for octahedral sphere coordination.
> 4. **Catalytic Salt Bridge (Lys42):** The $\epsilon$-amino group of Lys42 must form a hydrogen bond / salt bridge with oxygen atoms of the $\alpha$- and $\beta$-phosphates.

Once verified, assemble the disease-associated mutant (**L119R**) onto the exact same cofactor pose:

```bash
# 2. Assemble Holo Mutant (L119R) from the verified donor pose
cmda assemble --state atpmg --variant L119R --template 8FP5 --relax-sidechains
```

<details>
<summary><b>View mutant assembly directory tree</b></summary>

```text
data/atpmg/structure/O76039_4BGQ/variants/L119R/
├── L119R.pdb               # Assembled mutant, ATP and Mg
├── complex.pdb             # Copy of the assembled complex
└── ASSEMBLY_DONE           # Verified assembly and donor fit
```

</details>

---

### 1.2. Force Field Parameterization & Amber `tleap` Solvation

Kinase catalytic complexes require specialized nucleotide and divalent ion force fields. Standard generic GAFF parameters often cause artificial distortion or premature dissociation of triphosphate groups.

ColabMDA applies a dedicated biophysical force field matrix:

| Component | Force Field Model | Parameters / Citation | Biophysical Rationale |
| :--- | :--- | :--- | :--- |
| **Protein** | Amber ff14SB | Maier et al. (*JCTC* 2015) | Validated benchmark for kinase backbones & sidechains |
| **Solvent** | OPC (4-site) | Izadi et al. (*JPCB* 2014) | Accurate liquid density, bulk dielectric constant, and water structure |
| **Divalent $\text{Mg}^{2+}$** | Li–Merz 12–6 (OPC set) | Li, Roberts, Chakravorty, Merz (*JCTC* 2013) | Prevents artificial over-coordination and matches OPC solvation free energy |
| **Counterions ($\text{Na}^+, \text{Cl}^-$)** | Li–Merz 12–6 (OPC set) | Li & Merz (*Chem. Rev.* 2014) | Water-matched neutralizing counterions |
| **Nucleotide (ATP)** | Quantum RESP (ATP-B3) | Hu et al. (*J. Chem. Inf. Model.* 2024) | High-level B3LYP quantum-derived RESP charges preserving triphosphate geometry |
| **Ribose Dihedrals** | ParmEd CMAP Correction | Darden et al. | Refines dihedral energy surface to prevent artificial sugar ring puckering |

#### Holo Electrostatics & Net Charge Breakdown

Understanding system electrostatics is vital for evaluating parameterization integrity:
* **Protein Net Charge:** CDKL5 kinase domain (residues 1–303) at pH 7.4 has a net charge of **$-4$** (24 Asp, 21 Glu, 18 Lys, 23 Arg, with neutral Cys/Tyr and monoprotonated His).
* **Nucleotide (ATP) Net Charge:** Fully deprotonated triphosphate tail contributes **$-4$**.
* **Divalent Magnesium ($\text{Mg}^{2+}$):** Contributes **$+2$**.
* **Total Holo Complex Charge:** $-4 + (-4) + (+2) = \mathbf{-6}$.
* **System Neutralization:** To neutralize the complex and establish $150\text{ mM}$ physiological ionic strength, `tleap` adds **$44\text{ Na}^+$** and **$38\text{ Cl}^-$** ions ($44 - 38 = +6$).

Prepare the solvated systems locally:

```bash
# 1. Parameterize and solvate Wild-Type Holo system
cmda prep --state atpmg --variant WT

# 2. Parameterize and solvate Mutant Holo system
cmda prep --state atpmg --variant L119R
```

<details>
<summary><b>View stage directory tree</b></summary>

```text
data/atpmg/simulation/
├── WT/
│   └── prep/
│       ├── model.pdb                                # Prepared solute coordinates
│       ├── solvated.pdb                             # Protein, water and ions; includes ATP/Mg
│       ├── system.prmtop                            # Amber topology passed to equilibration
│       ├── system.inpcrd                            # Amber coordinates passed to equilibration
│       └── PREP_DONE                                # Verified preparation outputs
└── L119R/
    └── prep/
        ├── model.pdb                                # Prepared solute coordinates
        ├── solvated.pdb                             # Protein, water and ions; includes ATP/Mg
        ├── system.prmtop                            # Amber topology passed to equilibration
        ├── system.inpcrd                            # Amber coordinates passed to equilibration
        └── PREP_DONE                                # Verified preparation outputs
```

</details>

#### Atomic Breakdown: Solvated CDKL5 Holo Complex (WT)

| Metric | Value | Biophysical Description |
| :--- | :--- | :--- |
| **Total Atoms** | **45,862** | Full solvated orthogonal simulation box |
| **Protein Atoms** | 4,792 | CDKL5 residues 1–303 (all loops closed) |
| **ATP Ligand Atoms** | 43 | Adenosine triphosphate |
| **Magnesium Ions** | 1 | Coordinated $\text{Mg}^{2+}$ cation |
| **Water Molecules** | 13,648 | Rigid 4-site OPC waters ($13,648 \times 3 = 40,944$ atoms) |
| **Neutralizing / Salt Ions** | 44 $\text{Na}^+$, 38 $\text{Cl}^-$ | $150\text{ mM}$ physiological ionic strength |
| **Box Dimensions** | $83.2 \times 75.1 \times 89.0\text{ Å}$ | Minimum $14\text{ Å}$ solvent buffer to periodic boundaries |
| **Net System Charge** | **0.000** | Fully neutralized |

Visually confirm the solvated water boxes in PyMOL (WT and Mutant):
```bash
pymol data/atpmg/simulation/WT/prep/solvated.pdb
pymol data/atpmg/simulation/L119R/prep/solvated.pdb
```

---

## 🚀 Phase 2: Launch Colab GPU & Stream Data to Drive
*Now that your Holo models and solvated Amber parameters are built and verified locally, launch your Colab GPU instance and stream the coordinates to Google Drive.*

In your local terminal:

```bash
# 1. Launch a persistent Colab GPU instance
colab new -s cmda --gpu T4

# 2. Mount Google Drive (visit the printed OAuth URL once, then press Enter)
colab drivemount --session cmda

# 3. Stream local data/ directly to Google Drive in ONE single command:
cmda push -s cmda
```

<details>
<summary><b>View transferred directory tree</b></summary>

```text
/content/drive/MyDrive/research/projects/cmda-example/
└── data/atpmg/
    ├── structure/          # Transferred structural inputs and models
    └── simulation/         # Transferred prepared systems
```

</details>

> [!TIP]
> **Data Destination on Google Drive:**
> * **Google Drive Filesystem Path:** `/content/drive/MyDrive/research/projects/cmda-example/data/atpmg/`
> * **Mechanism:** `cmda push` packages your local `data/` folder, streams it to your active Colab VM (`-s cmda`), and extracts it directly into your mounted Google Drive project directory.

---

## ⚡ Phase 3: Cloud GPU Simulation & Analysis (Google Colab)
*Connect to Google Colab to run OpenMM molecular dynamics simulations on an NVIDIA GPU.*

### 3.1. Connect to Colab GPU Console
In your local terminal:
```bash
colab console -s cmda
```

Inside the Colab console terminal, run the cloud setup (~90 seconds) and navigate to the project directory:

```bash
# 1. Install cmda, OpenMM (CUDA/CPU), and AmberTools:
curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/main/scripts/setup.sh | bash

# 2. Activate environment & navigate to Drive project directory:
source /root/miniforge3/etc/profile.d/conda.sh
conda activate cmda
cd /content/drive/MyDrive/research/projects/cmda-example
```

---

### 3.2. Stage & Run 24-Stage Equilibration

In holo complexes, thermalization must restrain the nucleotide and divalent magnesium ion alongside the protein backbone to prevent active-site disruption before the solvent box reaches target density.

ColabMDA executes a **24-stage progressive equilibration ladder**:

| Stages | Purpose | Restraint Schedule ($k_{\text{bb}}$, $k_{\text{ATP}}$, $k_{\text{Mg}}$) | Ensemble | Duration & Targets |
| :--- | :--- | :--- | :--- | :--- |
| **1–6** | Energy Minimization & Initial Heating | $10.0\text{ kcal/mol}\cdot\text{Å}^2$ on Backbone, ATP, and $\text{Mg}^{2+}$ | NVT (Langevin) | 5,000 steps EM; heat 0 K $\rightarrow$ 100 K |
| **7–12** | Intermediate Thermalization | Lowered to $5.0\text{ kcal/mol}\cdot\text{Å}^2$ on active site | NVT | 100 K $\rightarrow$ 250 K (25 ps per stage) |
| **13–18** | Target Temperature & Density Matching | Lowered to $2.5\text{ kcal/mol}\cdot\text{Å}^2$ | NPT (1 bar Monte Carlo) | Reach 300 K; density converges toward $1.00\text{ g/cm}^3$ |
| **19–21** | Solvent Relaxation & Gradual Release | Lowered to $1.0\text{ kcal/mol}\cdot\text{Å}^2$ | NPT (300 K, 1 bar) | Solvent molecules fully organize in pocket |
| **22–24** | Unrestrained Final Relaxation | **$0.0\text{ kcal/mol}\cdot\text{Å}^2$** (Fully unrestrained) | NPT (300 K, 1 bar) | 100 ps unrestrained MD; ready for production |

```bash
# 1. Stage directories and run full 24-stage equilibration for WT (includes QC gating)
cmda stage --state atpmg --variant WT
cmda sim eq --state atpmg --variant WT

# 2. Stage directories and run full 24-stage equilibration for L119R (includes QC gating)
cmda stage --state atpmg --variant L119R
cmda sim eq --state atpmg --variant L119R

# (Optional) Re-run standalone QC gating without re-simulating:
# cmda sim qc --state atpmg --variant WT
# cmda sim qc --state atpmg --variant L119R
```

<details>
<summary><b>View stage directory tree</b></summary>

```text
data/atpmg/
├── simulation/
│   ├── WT/
│   │   └── eq/
│   │       ├── solvated.pdb                         # Topology coordinates used by OpenMM
│   │       ├── system.xml                           # Serialized OpenMM system
│   │       ├── em.state.xml                         # Minimized state
│   │       ├── nvt.state.xml                        # State after NVT equilibration
│   │       ├── npt.state.xml                        # Equilibrated state passed to production
│   │       ├── eq.pdb                               # Final equilibrated coordinates
│   │       └── stages/
│   │           ├── state_01.xml                     # First staged restart
│   │           └── state_24.xml                     # Last of 24 staged restart states
│   └── L119R/
│       └── eq/
│           ├── solvated.pdb                         # Topology coordinates used by OpenMM
│           ├── system.xml                           # Serialized OpenMM system
│           ├── em.state.xml                         # Minimized state
│           ├── nvt.state.xml                        # State after NVT equilibration
│           ├── npt.state.xml                        # Equilibrated state passed to production
│           ├── eq.pdb                               # Final equilibrated coordinates
│           └── stages/
│               ├── state_01.xml                     # First staged restart
│               └── state_24.xml                     # Last of 24 staged restart states
└── analysis/
    ├── WT/
    │   └── eq/
    │       ├── qc.txt                               # Equilibration diagnostics
    │       └── qc.png                               # Thermodynamic and structural QC figure
    └── L119R/
        └── eq/
            ├── qc.txt                               # Equilibration diagnostics
            └── qc.png                               # Thermodynamic and structural QC figure
```

</details>

---

### 3.2.1. Intermediate Inspection Gate: Minimal Verification Checklist (Local Machine)
*Before committing cloud compute to 100 ns production runs, verify that equilibration completed stably in 4 minimal steps:*

#### 1. Pull Deliverables to Local Machine (Terminal)
```bash
cd ~/research/projects/cmda-example

# Pull lightweight equilibration deliverables (~2 MB):
cmda pull -s cmda --stage eq

# Or pull only a specific variant:
# cmda pull -s cmda --stage eq --variant WT
```

<details>
<summary><b>View transferred directory tree</b></summary>

```text
data/atpmg/
├── simulation/             # Retrieved simulation deliverables
└── analysis/               # Retrieved equilibration QC
```

</details>

> [!TIP]
> **GPU vs. CPU Session for Retrieval:**
> * **Active Session:** If your GPU session is still running, pull directly.
> * **Session Ended / Timed Out:** You **do not** need a GPU to pull results. Simply launch an instant, free CPU runtime (`colab new -s cmda`), mount Drive (`colab drivemount --session cmda`), and run `cmda pull`.

#### 2. Verify Status in 1 Second (Terminal)
```bash
cat data/atpmg/analysis/WT/eq/qc.txt
cat data/atpmg/analysis/L119R/eq/qc.txt
```
*Expected Output:* `Status: PASSED`, `Cα RMSD: ~0.83 Å` (Confirms temperature $\approx 300\text{ K}$, density $\approx 1.01\text{ g/mL}$, fold intact).

#### 3. Glance at 6-Panel QC Plot in 5 Seconds
Open the unified 7-inch publication plot in any image viewer:
* `data/atpmg/analysis/WT/eq/qc.png`
* `data/atpmg/analysis/L119R/eq/qc.png`

*Expected Visual:*
* **Panels A–D (Thermodynamics):** Temperature is flat at **$300\text{ K}$**; density is flat at **$\sim 1.01\text{ g/mL}$**; potential energy and volume are stable.
* **Panel E (Per-Residue C$\alpha$ Drift):** Core secondary structures stay flat at **$< 0.5\text{ \AA}$**; loop peaks remain below threshold ($< 1.2\text{ \AA}$).
* **Panel F (3D Backbone Alignment):** Solid cyan equilibrated fold matches semi-transparent salmon starting model ($\text{C}_\alpha\text{ RMSD} < 1.5\text{ \AA}$).

#### 4. Quick 3D Fold & Active Site Check in PyMOL in 30 Seconds
Launch PyMOL with both starting and equilibrated coordinates:
```bash
pymol data/atpmg/simulation/WT/prep/model.pdb data/atpmg/simulation/WT/eq/eq.pdb
```
In the PyMOL right panel:
1. Next to **`(all)`**: Click **`H`** (Hide) $\to$ **`waters`**, then **`S`** (Show) $\to$ **`as`** $\to$ **`cartoon`**
2. Next to **`(all)`**: Click **`S`** (Show) $\to$ **`organic`** $\to$ **`sticks`** (ATP), then **`S`** (Show) $\to$ **`spheres`** ($\text{Mg}^{2+}$)
3. Next to **`eq`**: Click **`A`** (Action) $\to$ **`align`** $\to$ **`to molecule`** $\to$ choose **`model`**

*Expected Visual:* Cyan kinase backbone (`eq`) overlays closely onto starting model (`model`). ATP and $\text{Mg}^{2+}$ remain stably bound in the catalytic cleft without displacement ($\text{C}_\alpha\text{ RMSD} < 1.5\text{ \AA}$). Repeat for L119R.

<details>
<summary><b>🔍 Advanced PyMOL GUI Inspection (Periodic Box & Ions)</b></summary>

#### Fast PyMOL Inspection Commands
Type these commands directly into the PyMOL command prompt (`PyMOL>`):

```pml
# 1. Clean view & align backbones (salmon = starting model, cyan = equilibrated)
hide (resn HOH,WAT)
show cartoon
color salmon, model_
color cyan, eq
align eq, model_

# 2. Display the periodic unit cell boundary box wireframe
show cell

# 3. Display counter-ions and cofactors/ligands (ATP + Mg2+)
show spheres, (resn NA,CL,MG)
show sticks, (organic or resn ATP)
```

> [!NOTE]
> **Why are waters seen outside the cyan box in `eq.pdb`?**
> OpenMM and AMBER integrate in continuous Cartesian space. Molecules diffusing across boundaries remain unwrapped in raw coordinate snapshots. In **Phase 4.2**, **`cmda merge`** automatically runs PBC autoimaging (folding all molecules back inside the unit cell) and solvent stripping for clean trajectory animation.
> *(Note: `eq.pdb` is a static 3D coordinate snapshot; multi-frame animation is recorded in Phase 3.3 and viewed via `cmda view`).*

</details>

*Once verified, proceed to launch Production MD (Section 3.3).*

---

### 3.3. Production MD Across Triplicates (`md/{r1, r2, r3}`)

Independent replicates fork from the relaxed `eq/` state using unique random velocity seeds, ensuring true statistical independence across triplicates:

```bash
# Wild-Type Holo Triplicates (100 ns production runs)
cmda sim md --state atpmg --variant WT --replica r1 --total-ns 100
cmda sim md --state atpmg --variant WT --replica r2 --total-ns 100
cmda sim md --state atpmg --variant WT --replica r3 --total-ns 100

# Mutant Holo Triplicates (100 ns production runs)
cmda sim md --state atpmg --variant L119R --replica r1 --total-ns 100
cmda sim md --state atpmg --variant L119R --replica r2 --total-ns 100
cmda sim md --state atpmg --variant L119R --replica r3 --total-ns 100
```

<details>
<summary><b>View stage directory tree</b></summary>

```text
data/atpmg/simulation/
├── WT/
│   └── md/
│       └── r1/
│           ├── system.xml                           # Copied OpenMM system
│           ├── solvated.pdb                         # Copied topology coordinates
│           ├── npt.state.xml                        # Starting state copied from equilibration
│           ├── prod.chk                             # Checkpoint used to resume production
│           └── chunks/
│               ├── 0001.dcd                         # First committed trajectory chunk
│               └── 0001.log                         # Matching committed thermodynamic log
└── L119R/
    └── md/
        └── r1/
            ├── system.xml                           # Copied OpenMM system
            ├── solvated.pdb                         # Copied topology coordinates
            ├── npt.state.xml                        # Starting state copied from equilibration
            ├── prod.chk                             # Checkpoint used to resume production
            └── chunks/
                ├── 0001.dcd                         # First committed trajectory chunk
                └── 0001.log                         # Matching committed thermodynamic log
```

</details>

Each replica has its own folder. Committed chunks are numbered `0001`, `0002`,
and so on; `.tmp` files belong to the chunk currently being written. Replicas
`r2` and `r3` use the same layout when run.

#### Simulation Parameters & Frame Saving Cadence

Production dynamics uses standard biophysical integration parameters optimized for cloud stability and storage efficiency:

| Parameter | Value | Biophysical Rationale & Details |
| :--- | :--- | :--- |
| **Timestep ($\Delta t$)** | **$2\text{ fs}$ ($0.002\text{ ps}$)** | Standard Langevin middle integrator with SHAKE/SETTLE constraints on all bonds involving hydrogen. |
| **Frame Saving Interval** | **Every $10\text{ ps}$ ($5,000\text{ steps}$)** | Writes coordinates to `.dcd` trajectory every 10 ps (`--traj-interval 10`), yielding high temporal resolution. |
| **Thermodynamic Log Interval** | **Every $10\text{ ps}$ ($5,000\text{ steps}$)** | Records instantaneous temperature, potential energy, kinetic energy, volume, and density to `.log`. |
| **Chunk Partitioning** | **$1\text{ ns}$ ($500,000\text{ steps}$)** | Simulations are partitioned into $1\text{ ns}$ self-contained chunks ($100\text{ frames}$ per chunk file). |
| **Checkpoint Cadence** | **Every $1\text{ ns}$** | Serializes binary state checkpoint `prod.chk` to Google Drive at the end of each chunk. |
| **Total Frames Generated** | **$10,000\text{ frames}$** | $\frac{100,000\text{ ps}}{10\text{ ps}} = 10,000\text{ frames}$ per $100\text{ ns}$ replica. |
| **Storage Footprint** | **~5.5 GB $\to$ ~580 MB** | Raw solvated chunks total ~5.5 GB; stripped holo trajectories (`cmda merge`) reduce this by >90% to ~580 MB. |

#### Downstream Handover

> **Downstream Handover:** These 100 raw 1-ns chunk files are handed over to the merge stage (`cmda merge`), which joins them, eliminates periodic boundary condition (PBC) box jumps, and strips bulk water into `traj.dcd` (preserving CDKL5 + ATP + Mg²⁺) for active-site and conformational analysis.

| Strategy | Trajectory Contents | Saving Frequency | Frames per 100 ns | Size per Replica |
| :--- | :--- | :--- | :--- | :--- |
| **Unoptimized Raw** | Full System (Water + Ions, 45,862 atoms) | 1 ps (100,000 frames) | 100,000 | **~55.0 GB** |
| **Standard Raw** | Full System (Water + Ions, 45,862 atoms) | 10 ps (10,000 frames) | 10,000 | **~5.5 GB** |
| **ColabMDA Stripped Holo** | **Protein + ATP + $\text{Mg}^{2+}$ (4,836 atoms)** | **10 ps (10,000 frames)** | **10,000** | **~580 MB (90% savings)** |

*Result:* An entire 3-replica triplicate study (WT + L119R = 6 trajectories) requires **under 3.5 GB of Google Drive storage**, fitting comfortably within free 15 GB tiers!

---

### 3.4. Continuous Monitoring & Timeout Recovery Runbook

Monitor simulation progress from your local workstation:

```bash
cmda status --state atpmg --variant WT --replica r1
```

<details>
<summary>🔍 <b>Click to view example Holo status report</b></summary>

```text
[STATUS]
  Workdir          : data/atpmg/simulation/WT/md/r1
  Chunks (DCD/log) : 100 / 100
  Topology File    : solvated.pdb
                     └─ 45862 atoms, 13952 residues
                     └─ (303 protein, 1 ATP, 1 MG, 13648 water, 82 ions)
  Trajectory File  : traj.dcd (10000 frames)
  Log File         : md.log (10000 frames)
  Frames (from logs): 10000
  Sim time (ps/ns) : 100000.00 ps / 100.0000 ns (from chunk names)
  Resume-ready     : YES (needs prod.chk, system.xml, solvated.pdb)
```

</details>

#### 💡 How to Resume After a Google Colab Timeout
Production MD writes periodic checkpoints (`prod.chk`) directly to Google Drive. If your session expires:

1. Re-launch the session and re-mount Drive from your local terminal:
   ```bash
   colab new -s cmda --gpu T4
   colab drivemount --session cmda
   colab console -s cmda
   ```
2. Re-run the exact same production command:
   ```bash
   cd /content/drive/MyDrive/research/projects/cmda-example
   cmda sim md --state atpmg --variant WT --replica r1 --total-ns 100
   ```
   *ColabMDA detects the existing `prod.chk` checkpoint on Google Drive and resumes from the last saved checkpoint.*

#### Python Trajectory & Catalytic Coordination Health Check

To quickly verify trajectory health and check that $\text{Mg}^{2+}$ remained chelated to ATP during a chunk, run this snippet inside the Colab console:

```python
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array

# Load topology and trajectory chunk
u = mda.Universe("data/atpmg/simulation/WT/prep/system.prmtop",
                 "data/atpmg/simulation/WT/md/r1/chunks/0001.dcd")
print(f"Loaded {len(u.trajectory)} frames successfully.")

# Select magnesium and phosphate oxygens
mg = u.select_atoms("resname MG and name MG")
beta_gamma_oxygens = u.select_atoms("resname ATP and name O1B O2B O1G O2G O3G")

# Calculate minimum distance in final frame
u.trajectory[-1]
coord_dist = distance_array(mg.positions, beta_gamma_oxygens.positions).min()
print(f"Final Frame Mg²⁺–phosphate distance: {coord_dist:.2f} Å (Expected: 1.9–2.2 Å)")
assert 1.8 <= coord_dist <= 2.4, "ALERT: Magnesium drifted from ATP triphosphate chain!"
```

---

## 📥 Phase 4: Retrieve Results, Trajectory Concatenation & Publication Figures (Local Workstation)
*When simulations finish, pull the analysis tables, logs, and publication figures back to your laptop in seconds.*

### 4.1. Retrieve Results in ONE Command

In your local terminal:
```bash
# Pull analysis results, figures, and logs in seconds:
cmda pull -s cmda

# To also pull post-processed stripped trajectories for PyMOL inspection:
# cmda pull -s cmda --traj

# (Optional: to pull full raw simulation chunks as well):
# cmda pull -s cmda --all
```

<details>
<summary><b>View transferred directory tree</b></summary>

```text
data/atpmg/
├── simulation/             # Retrieved simulation deliverables
└── analysis/               # Retrieved QC, numerical tables and figures
```

</details>

---

### 4.2. Multi-Segment Trajectory Concatenation & Autoimaging

Concatenate chunked production trajectories and center the catalytic complex:

```bash
# Merge WT Holo trajectories
cmda merge --state atpmg --variant WT

# Merge Mutant Holo trajectories
cmda merge --state atpmg --variant L119R
```

<details>
<summary><b>View stage directory tree</b></summary>

```text
data/atpmg/simulation/
├── WT/
│   └── md/
│       └── r1/
│           ├── traj.dcd                 # Autoimaged trajectory without bulk solvent
│           ├── top.prmtop                      # Matching topology; use with the trajectory above
│           ├── trajCa.dcd                          # Cα trajectory
│           ├── topCa.prmtop                            # Matching Cα topology
│           ├── md.log                               # Copy of the merged production log
│           └── MERGE_DONE                           # Verified merged outputs
└── L119R/
    └── md/
        └── r1/
            ├── traj.dcd                 # Autoimaged trajectory without bulk solvent
            ├── top.prmtop                      # Matching topology; use with the trajectory above
            ├── trajCa.dcd                          # Cα trajectory
            ├── topCa.prmtop                            # Matching Cα topology
            ├── md.log                               # Copy of the merged production log
            └── MERGE_DONE                           # Verified merged outputs
```

</details>

---

### 4.3. ⚠️ CRITICAL: Holo Trajectory vs. Topology Pairing Safety Rule

> [!WARNING]
> ### Do Not Pair Holo Trajectories with Water Topologies or Apo Topologies!
> During `cmda merge --state atpmg`, ColabMDA strips bulk water and counterions but **strictly preserves the catalytic complex: Protein + ATP + $\text{Mg}^{2+}$**.
>
> 1. **Stripped Atom Count:** The resulting trajectory `traj.dcd` contains atoms for **CDKL5 (4,792 atoms) + ATP (43 atoms) + $\text{Mg}^{2+}$ (1 atom) = 4,836 atoms**.
> 2. **Correct Topology Pairing:** Always pair `traj.dcd` with:
>    ```text
>    data/atpmg/simulation/WT/md/r1/top.prmtop
>    ```
> 3. **Common Fatal Errors:**
>    * **Crash A:** Pairing `traj.dcd` with `system.prmtop` (which contains ~45,000 atoms including water) will cause MDAnalysis/VMD/PyMOL to throw:
>      `ValueError: Topology and trajectory atom counts do not match (4836 vs 45862)`.
>    * **Crash B:** Pairing `traj.dcd` with Apo's `top.prmtop` (which contains only 4,792 protein atoms and lacks ATP/Mg) will cause:
>      `ValueError: Topology and trajectory atom counts do not match (4836 vs 4792)`.

---

### 4.4. Automated Triplicate Comparative Analysis & Biophysical Interpretation Guide

Generate publication-grade comparative figures evaluating how the disease mutation affects active-site stability:

```bash
# 1. Generate single-variant analysis tables
cmda analyze --state atpmg --variant WT
cmda analyze --state atpmg --variant L119R
```

<details>
<summary><b>View stage directory tree</b></summary>

```text
data/atpmg/analysis/
├── WT/
│   └── r1/
│       ├── rmsd.csv                                 # RMSD values
│       ├── rg.csv                                   # Radius of gyration values
│       ├── rmsf.csv                                 # Per-residue fluctuations
│       ├── figures/
│       │   ├── rmsd_vs_time.png                     # RMSD plot
│       │   ├── rg_vs_time.png                       # Radius of gyration plot
│       │   └── rmsf_per_residue.png                 # Fluctuation plot
│       └── ANALYSIS_DONE                            # Verified replica tables and figures
└── L119R/
    └── r1/
        ├── rmsd.csv                                 # RMSD values
        ├── rg.csv                                   # Radius of gyration values
        ├── rmsf.csv                                 # Per-residue fluctuations
        ├── figures/
        │   ├── rmsd_vs_time.png                     # RMSD plot
        │   ├── rg_vs_time.png                       # Radius of gyration plot
        │   └── rmsf_per_residue.png                 # Fluctuation plot
        └── ANALYSIS_DONE                            # Verified replica tables and figures
```

</details>

```bash
# 2. Generate publication-grade comparative overlay plots (Mean ± SD)
cmda compare --state atpmg --variants WT,L119R
```

<details>
<summary><b>View stage directory tree</b></summary>

```text
data/atpmg/analysis/compare/
├── rmsd.png                                         # Variant RMSD comparison
├── rg.png                                           # Variant compactness comparison
├── rmsf.png                                         # Variant fluctuation comparison
├── figures/
│   ├── rmsd_avg.png                                 # Copy of RMSD comparison
│   ├── rg_avg.png                                   # Copy of compactness comparison
│   └── rmsf_avg.png                                 # Copy of fluctuation comparison
└── summary.txt                                      # Comparison statistics
```

</details>

This generates publication-ready 300 DPI figures under `data/atpmg/analysis/compare/`:
* `summary.txt` (Executive summary table with Mean ± SD for RMSD, Rg, Temp)
* `rmsd.png` (Backbone and ligand stability time-series overlay)
* `rg.png` (Radius of Gyration domain closure comparison)
* `rmsf.png` (Per-residue fluctuation profiles highlighting P-loop stiffening)

#### How to Interpret the Comparative Figures:

1. **Backbone & Ligand RMSD Profiles (`rmsd.png`):**
   * **Active Site Fidelity:** Assesses whether ATP remains firmly bound within the catalytic pocket. A stable Holo complex displays a plateaued ligand RMSD ($< 1.8\text{ Å}$ relative to the initial superposed pose).
   * **Comparative Drift:** If the L119R mutant destabilizes the binding cleft, the ligand RMSD will exhibit increased drift or large standard-deviation fluctuations.

2. **Radius of Gyration Compaction (`rg_avg.png`):**
   * **Domain Closure:** Protein kinases undergo a characteristic conformational transition where the N-lobe and C-lobe rotate toward each other to clamp over the nucleotide. In Holo simulations, the average $R_g$ is typically **$0.3–0.5\text{ Å}$ lower** than in Apo simulations, confirming stable catalytic domain closure.

3. **Per-Residue Dynamic Flexibility (`rmsf_avg.png`):**
   * **P-Loop Stiffening:** In Apo simulations, the glycine-rich P-loop (residues 18–26) is highly flexible ($\text{RMSF} \sim 2.0–2.8\text{ Å}$). In Holo simulations, coordination with the triphosphate tail dampens these fluctuations ($\text{RMSF} \sim 0.8–1.2\text{ Å}$).
   * **Activation Loop & $\alpha\text{C}$-Helix:** Identifies whether the mutation alters dynamics in the activation segment (residues 163–180) or weakens the conserved salt bridge between Lys42 and Glu60.

---

### 4.5. 3D Trajectory Viewing & Pocket-Focused Snapshot Grids in PyMOL

Run these commands in the project directory on a computer with desktop PyMOL.
For a Colab run, first download the processed trajectories using `cmda pull --traj`.
Open one variant at a time; close its PyMOL window before running the next command.

```bash
cmda view --state atpmg --variant WT --replica r1 --frame all
cmda view --state atpmg --variant L119R --replica r1 --frame all
```

PyMOL opens and starts playing all saved trajectory frames. Use its playback
controls to pause, play or step through frames. `--frame all` displays the movie;
`--frame 1` displays only the first frame.

The command reads the matching `md/r1/top.prmtop` and `md/r1/traj.dcd` for the
selected state and variant. It generates `view.pdb`, `visualize.pml` and
`visualize.log` in that replica's analysis folder. The PDB has the same atoms as
the processed trajectory, so loading is handled for you.

<details>
<summary><b>Manual loading in PyMOL</b></summary>

After the command has generated the viewing files, these commands can also be
entered in PyMOL from the project directory:

```pml
load data/atpmg/analysis/WT/r1/view.pdb, protein
load_traj data/atpmg/simulation/WT/md/r1/traj.dcd, protein, state=1
hide everything
show cartoon, protein
mplay
```

For another variant or replica, change both paths together. Use `solvated.pdb`
only with the original full-system chunks, whose atom count includes water and
salt. The processed `traj.dcd` matches `top.prmtop` and the generated `view.pdb`.

</details>

#### Optional snapshot image

To save a comparison of the first and last frames instead of playing a movie:

```bash
cmda snapshots --state atpmg --variant WT,L119R --replica r1 --frames 1,last --resi 119
```

This writes a PNG image and a plain-text provenance log in the analysis figures
folder. Each variant forms one row; the first and last frames form two columns.

<details>
<summary><b>Existing runs with older filenames</b></summary>

After production has finished, rename existing processed files without rerunning
or rewriting the trajectories:

```bash
cmda merge --state atpmg --variant WT --replica r1 --rename-files
```

Repeat for each completed variant and replica. The command validates matching
atom and frame counts before renaming, rejects conflicting files, and leaves raw
chunks, equilibration states and restart checkpoints in place. Older names remain
readable when working with historical runs.

</details>

## Complete project directory structure

<details>
<summary><b>View the complete workflow directory tree</b></summary>

The tree covers the workflow outputs; MODELLER and preparation tools may also
leave auxiliary files in their owning folders. Raw archives may use `.pdb`
instead of `.cif` when the archive fallback is used. The `r2`, `r3` and
`aggregate` entries apply after all three replicas have been run and analyzed. The prerequisite
apo models remain under `data/apo/` as described in the apo tutorial.

```text
~/research/projects/cmda-example/
└── data/
    └── atpmg/
        ├── structure/
        │   └── O76039_4BGQ/
        │       ├── donors/
        │       │   └── 8FP5/
        │       │       └── 8fp5.cif                 # Cofactor donor archive; PDB fallback is possible
        │       └── variants/
        │           ├── WT/
        │           │   ├── WT.pdb                   # Assembled protein, ATP and Mg
        │           │   ├── complex.pdb              # Copy of the assembled complex
        │           │   └── ASSEMBLY_DONE            # Verified assembly and donor fit
        │           └── L119R/
        │               ├── L119R.pdb                # Assembled protein, ATP and Mg
        │               ├── complex.pdb              # Copy of the assembled complex
        │               └── ASSEMBLY_DONE            # Verified assembly and donor fit
        ├── simulation/
        │   ├── WT/
        │   │   ├── prep/
        │   │   │   ├── model.pdb                    # Prepared solute coordinates
        │   │   │   ├── solvated.pdb                 # Protein, water and ions; includes ATP/Mg
        │   │   │   ├── system.prmtop                # Amber topology passed to equilibration
        │   │   │   ├── system.inpcrd                # Amber coordinates passed to equilibration
        │   │   │   └── PREP_DONE                    # Verified preparation outputs
        │   │   ├── eq/
        │   │   │   ├── solvated.pdb                 # Topology coordinates used by OpenMM
        │   │   │   ├── system.xml                   # Serialized OpenMM system
        │   │   │   ├── em.state.xml                 # Minimized state
        │   │   │   ├── nvt.state.xml                # State after NVT equilibration
        │   │   │   ├── npt.state.xml                # Equilibrated state passed to production
        │   │   │   ├── eq.pdb                       # Final equilibrated coordinates
        │   │   │   └── stages/
        │   │   │       ├── state_01.xml             # First staged restart
        │   │   │       └── state_24.xml             # Last of 24 staged restart states
        │   │   └── md/
        │   │       ├── r1/
        │   │       │   ├── system.xml               # Copied OpenMM system
        │   │       │   ├── solvated.pdb             # Copied topology coordinates
        │   │       │   ├── npt.state.xml            # Starting state copied from equilibration
        │   │       │   ├── prod.chk                 # Checkpoint used to resume production
        │   │       │   ├── chunks/
        │   │       │   │   ├── 0001.dcd             # First committed trajectory chunk
        │   │       │   │   └── 0001.log             # Matching committed thermodynamic log
        │   │       │   ├── traj.dcd     # Autoimaged trajectory without bulk solvent
        │   │       │   ├── top.prmtop          # Matching topology; use with the trajectory above
        │   │       │   ├── trajCa.dcd              # Cα trajectory
        │   │       │   ├── topCa.prmtop                # Matching Cα topology
        │   │       │   ├── md.log                   # Copy of the merged production log
        │   │       │   └── MERGE_DONE               # Verified merged outputs
        │   │       ├── r2/
        │   │       │   ├── system.xml               # Copied OpenMM system
        │   │       │   ├── solvated.pdb             # Copied topology coordinates
        │   │       │   ├── npt.state.xml            # Starting state copied from equilibration
        │   │       │   ├── prod.chk                 # Checkpoint used to resume production
        │   │       │   ├── chunks/
        │   │       │   │   ├── 0001.dcd             # First committed trajectory chunk
        │   │       │   │   └── 0001.log             # Matching committed thermodynamic log
        │   │       │   ├── traj.dcd     # Autoimaged trajectory without bulk solvent
        │   │       │   ├── top.prmtop          # Matching topology; use with the trajectory above
        │   │       │   ├── trajCa.dcd              # Cα trajectory
        │   │       │   ├── topCa.prmtop                # Matching Cα topology
        │   │       │   ├── md.log                   # Copy of the merged production log
        │   │       │   └── MERGE_DONE               # Verified merged outputs
        │   │       └── r3/
        │   │           ├── system.xml               # Copied OpenMM system
        │   │           ├── solvated.pdb             # Copied topology coordinates
        │   │           ├── npt.state.xml            # Starting state copied from equilibration
        │   │           ├── prod.chk                 # Checkpoint used to resume production
        │   │           ├── chunks/
        │   │           │   ├── 0001.dcd             # First committed trajectory chunk
        │   │           │   └── 0001.log             # Matching committed thermodynamic log
        │   │           ├── traj.dcd     # Autoimaged trajectory without bulk solvent
        │   │           ├── top.prmtop          # Matching topology; use with the trajectory above
        │   │           ├── trajCa.dcd              # Cα trajectory
        │   │           ├── topCa.prmtop                # Matching Cα topology
        │   │           ├── md.log                   # Copy of the merged production log
        │   │           └── MERGE_DONE               # Verified merged outputs
        │   └── L119R/
        │       ├── prep/
        │       │   ├── model.pdb                    # Prepared solute coordinates
        │       │   ├── solvated.pdb                 # Protein, water and ions; includes ATP/Mg
        │       │   ├── system.prmtop                # Amber topology passed to equilibration
        │       │   ├── system.inpcrd                # Amber coordinates passed to equilibration
        │       │   └── PREP_DONE                    # Verified preparation outputs
        │       ├── eq/
        │       │   ├── solvated.pdb                 # Topology coordinates used by OpenMM
        │       │   ├── system.xml                   # Serialized OpenMM system
        │       │   ├── em.state.xml                 # Minimized state
        │       │   ├── nvt.state.xml                # State after NVT equilibration
        │       │   ├── npt.state.xml                # Equilibrated state passed to production
        │       │   ├── eq.pdb                       # Final equilibrated coordinates
        │       │   └── stages/
        │       │       ├── state_01.xml             # First staged restart
        │       │       └── state_24.xml             # Last of 24 staged restart states
        │       └── md/
        │           ├── r1/
        │           │   ├── system.xml               # Copied OpenMM system
        │           │   ├── solvated.pdb             # Copied topology coordinates
        │           │   ├── npt.state.xml            # Starting state copied from equilibration
        │           │   ├── prod.chk                 # Checkpoint used to resume production
        │           │   ├── chunks/
        │           │   │   ├── 0001.dcd             # First committed trajectory chunk
        │           │   │   └── 0001.log             # Matching committed thermodynamic log
        │           │   ├── traj.dcd     # Autoimaged trajectory without bulk solvent
        │           │   ├── top.prmtop          # Matching topology; use with the trajectory above
        │           │   ├── trajCa.dcd              # Cα trajectory
        │           │   ├── topCa.prmtop                # Matching Cα topology
        │           │   ├── md.log                   # Copy of the merged production log
        │           │   └── MERGE_DONE               # Verified merged outputs
        │           ├── r2/
        │           │   ├── system.xml               # Copied OpenMM system
        │           │   ├── solvated.pdb             # Copied topology coordinates
        │           │   ├── npt.state.xml            # Starting state copied from equilibration
        │           │   ├── prod.chk                 # Checkpoint used to resume production
        │           │   ├── chunks/
        │           │   │   ├── 0001.dcd             # First committed trajectory chunk
        │           │   │   └── 0001.log             # Matching committed thermodynamic log
        │           │   ├── traj.dcd     # Autoimaged trajectory without bulk solvent
        │           │   ├── top.prmtop          # Matching topology; use with the trajectory above
        │           │   ├── trajCa.dcd              # Cα trajectory
        │           │   ├── topCa.prmtop                # Matching Cα topology
        │           │   ├── md.log                   # Copy of the merged production log
        │           │   └── MERGE_DONE               # Verified merged outputs
        │           └── r3/
        │               ├── system.xml               # Copied OpenMM system
        │               ├── solvated.pdb             # Copied topology coordinates
        │               ├── npt.state.xml            # Starting state copied from equilibration
        │               ├── prod.chk                 # Checkpoint used to resume production
        │               ├── chunks/
        │               │   ├── 0001.dcd             # First committed trajectory chunk
        │               │   └── 0001.log             # Matching committed thermodynamic log
        │               ├── traj.dcd     # Autoimaged trajectory without bulk solvent
        │               ├── top.prmtop          # Matching topology; use with the trajectory above
        │               ├── trajCa.dcd              # Cα trajectory
        │               ├── topCa.prmtop                # Matching Cα topology
        │               ├── md.log                   # Copy of the merged production log
        │               └── MERGE_DONE               # Verified merged outputs
        └── analysis/
            ├── WT/
            │   ├── eq/
            │   │   ├── qc.txt                       # Equilibration diagnostics
            │   │   └── qc.png                       # Thermodynamic and structural QC figure
            │   ├── r1/
            │   │   ├── rmsd.csv                     # RMSD values
            │   │   ├── rg.csv                       # Radius of gyration values
            │   │   ├── rmsf.csv                     # Per-residue fluctuations
            │   │   ├── figures/
            │   │   │   ├── rmsd_vs_time.png         # RMSD plot
            │   │   │   ├── rg_vs_time.png           # Radius of gyration plot
            │   │   │   └── rmsf_per_residue.png     # Fluctuation plot
            │   │   └── ANALYSIS_DONE                # Verified replica tables and figures
            │   ├── r2/
            │   │   ├── rmsd.csv                     # RMSD values
            │   │   ├── rg.csv                       # Radius of gyration values
            │   │   ├── rmsf.csv                     # Per-residue fluctuations
            │   │   ├── figures/
            │   │   │   ├── rmsd_vs_time.png         # RMSD plot
            │   │   │   ├── rg_vs_time.png           # Radius of gyration plot
            │   │   │   └── rmsf_per_residue.png     # Fluctuation plot
            │   │   └── ANALYSIS_DONE                # Verified replica tables and figures
            │   ├── r3/
            │   │   ├── rmsd.csv                     # RMSD values
            │   │   ├── rg.csv                       # Radius of gyration values
            │   │   ├── rmsf.csv                     # Per-residue fluctuations
            │   │   ├── figures/
            │   │   │   ├── rmsd_vs_time.png         # RMSD plot
            │   │   │   ├── rg_vs_time.png           # Radius of gyration plot
            │   │   │   └── rmsf_per_residue.png     # Fluctuation plot
            │   │   └── ANALYSIS_DONE                # Verified replica tables and figures
            │   └── aggregate/
            │       ├── rmsd_avg.csv                 # Replica mean and standard deviation
            │       ├── rg_avg.csv                   # Replica mean and standard deviation
            │       ├── rmsf_avg.csv                 # Replica mean and standard deviation
            │       ├── figures/
            │       │   ├── rmsd_aggregate.png       # Replica RMSD summary
            │       │   ├── rg_aggregate.png         # Replica compactness summary
            │       │   └── rmsf_aggregate.png       # Replica fluctuation summary
            │       └── AGGREGATE_DONE               # Verified multi-replica summary
            ├── L119R/
            │   ├── eq/
            │   │   ├── qc.txt                       # Equilibration diagnostics
            │   │   └── qc.png                       # Thermodynamic and structural QC figure
            │   ├── r1/
            │   │   ├── rmsd.csv                     # RMSD values
            │   │   ├── rg.csv                       # Radius of gyration values
            │   │   ├── rmsf.csv                     # Per-residue fluctuations
            │   │   ├── figures/
            │   │   │   ├── rmsd_vs_time.png         # RMSD plot
            │   │   │   ├── rg_vs_time.png           # Radius of gyration plot
            │   │   │   └── rmsf_per_residue.png     # Fluctuation plot
            │   │   └── ANALYSIS_DONE                # Verified replica tables and figures
            │   ├── r2/
            │   │   ├── rmsd.csv                     # RMSD values
            │   │   ├── rg.csv                       # Radius of gyration values
            │   │   ├── rmsf.csv                     # Per-residue fluctuations
            │   │   ├── figures/
            │   │   │   ├── rmsd_vs_time.png         # RMSD plot
            │   │   │   ├── rg_vs_time.png           # Radius of gyration plot
            │   │   │   └── rmsf_per_residue.png     # Fluctuation plot
            │   │   └── ANALYSIS_DONE                # Verified replica tables and figures
            │   ├── r3/
            │   │   ├── rmsd.csv                     # RMSD values
            │   │   ├── rg.csv                       # Radius of gyration values
            │   │   ├── rmsf.csv                     # Per-residue fluctuations
            │   │   ├── figures/
            │   │   │   ├── rmsd_vs_time.png         # RMSD plot
            │   │   │   ├── rg_vs_time.png           # Radius of gyration plot
            │   │   │   └── rmsf_per_residue.png     # Fluctuation plot
            │   │   └── ANALYSIS_DONE                # Verified replica tables and figures
            │   └── aggregate/
            │       ├── rmsd_avg.csv                 # Replica mean and standard deviation
            │       ├── rg_avg.csv                   # Replica mean and standard deviation
            │       ├── rmsf_avg.csv                 # Replica mean and standard deviation
            │       ├── figures/
            │       │   ├── rmsd_aggregate.png       # Replica RMSD summary
            │       │   ├── rg_aggregate.png         # Replica compactness summary
            │       │   └── rmsf_aggregate.png       # Replica fluctuation summary
            │       └── AGGREGATE_DONE               # Verified multi-replica summary
            └── compare/
                ├── rmsd.png                         # Variant RMSD comparison
                ├── rg.png                           # Variant compactness comparison
                ├── rmsf.png                         # Variant fluctuation comparison
                ├── figures/
                │   ├── rmsd_avg.png                 # Copy of RMSD comparison
                │   ├── rg_avg.png                   # Copy of compactness comparison
                │   └── rmsf_avg.png                 # Copy of fluctuation comparison
                └── summary.txt                      # Comparison statistics
```

</details>

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
