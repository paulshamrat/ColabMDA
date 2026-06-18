# PyMOL Visualization Guide & Snapshot Composition

This guide details how to set up the local shared PyMOL visualization environment and run the automated publication-ready snapshot grid generation commands in `ColabMDA`.

---

## 1. Conda Environment Setup (`pymol_env`)

To perform high-fidelity image rendering, PyMOL scripting, and snapshot composition locally, you must set up a shared conda environment (`pymol_env`) containing the correct version of PyMOL open-source along with analysis/image processing libraries.

```bash
# 1. Create the environment with PyMOL and core packages from conda-forge
conda create -y -n pymol_env -c conda-forge python=3.11 pymol-open-source=3.1.0 pillow=12.1.1 numpy=2.4.2 biopython=1.86

# 2. Install analysis and formatting dependencies via pip
conda run -n pymol_env pip install mdanalysis==2.10.0 scipy==1.17.1 svgwrite==1.4.3 tqdm==4.67.3

# 3. Install ColabMDA in editable mode inside the environment
conda run -n pymol_env pip install -e .
```

---

## 2. Automated Snapshot Grid Generation

The `snapshots` command renders a configurable grid of trajectory frames so you can compare structural changes across systems, variants, or replicas using the same camera, alignment, and styling rules.

### Running the Command
To generate a snapshot grid, create a custom JSON configuration file mapping your local trajectory files and run:

```bash
colabmda openmm snapshots -c scratch/snapshots_config.json -o scratch/master_3x8_snapshots.png --temp-dir scratch/master_grid_3x8
```

### Configuration Format (`snapshots_config.json`)
The JSON file maps system names, topologies, trajectories, and key frame states/times:

```json
{
    "camera_view": [
        -0.597478449, -0.617595792, -0.511463583,
        -0.502994895,  0.785388291, -0.360776752,
         0.624511778,  0.041705273, -0.779902637,
         0.000050262, -0.000109358, -182.850265503,
        31.930021286, 32.076828003,  34.52369679,
       141.553222656, 224.141143799, -20.000000000
    ],
    "align_ref_pdb": "analysis/clean_trajectories/g12c_protein.pdb",
    "stable_core_sel": "name CA and (resi 1-24 or resi 41-54 or resi 81-166)",
    "systems": {
        "WT": {
            "pdb": "analysis/clean_trajectories/wt_protein.pdb",
            "dcd": "analysis/clean_trajectories/wt_protein.dcd",
            "states": [1, 51, 91, 100, 101, 106, 111, 121],
            "times": ["t = 0.00 ns", "t = 0.50 ns", "t = 0.90 ns", "t = 0.99 ns", "t = 1.00 ns", "t = 1.05 ns", "t = 1.10 ns", "t = 1.20 ns"],
            "mut_residue": 12
        },
        "G12C": {
            "pdb": "analysis/clean_trajectories/g12c_protein.pdb",
            "dcd": "analysis/clean_trajectories/g12c_protein.dcd",
            "states": [1, 51, 91, 100, 101, 106, 111, 121],
            "times": ["t = 0.00 ns", "t = 0.50 ns", "t = 0.90 ns", "t = 0.99 ns", "t = 1.00 ns", "t = 1.05 ns", "t = 1.10 ns", "t = 1.20 ns"],
            "mut_residue": 12
        },
        "G12D": {
            "pdb": "analysis/clean_trajectories/g12d_protein.pdb",
            "dcd": "analysis/clean_trajectories/g12d_protein.dcd",
            "states": [1, 51, 91, 100, 101, 106, 111, 121],
            "times": ["t = 0.00 ns", "t = 0.50 ns", "t = 0.90 ns", "t = 0.99 ns", "t = 1.00 ns", "t = 1.05 ns", "t = 1.10 ns", "t = 1.20 ns"],
            "mut_residue": 12
        }
    }
}
```

---

## 3. High-Quality Rendering Style & Key Commands

To reproduce consistent publication-quality visuals, the command runs headless PyMOL (`pymol -qc`) and executes the following visualization script logic:

### 3.1. Structure Alignment & Fitting
```text
# 1. Fit the trajectory frame-by-frame to the first frame using Cα atoms of the stable core
# This locks the core G-domain in place to make loop fluctuations easy to track.
intra_fit prot and name CA and (resi 1-24 or resi 41-54 or resi 81-166), 1

# 2. Align the overall structure to the common reference structure
align prot and stable_core, ref_align and stable_core, mobile_state=1, target_state=1
```

### 3.2. Color Coding and Transparency
```text
# 1. Base cartoon setup (semi-transparent gray for non-mutant and non-loop regions)
color gray90, prot
set cartoon_transparency, 0.6, prot and not (resid 57-75 or resid 12)

# 2. Highlight Switch II loop (residues 57-75) in solid opaque orange
color orange, prot and resid 57-75
set cartoon_transparency, 0.0, resid 57-75

# 3. Highlight mutated Cysteine/Aspartate 12 in solid yellow sticks
color yellow, prot and resid 12
show sticks, prot and resi 12 and not name N+C+O+H
set stick_radius, 0.25
```

### 3.3. Ray Tracing and Shadow Styles
```text
set ray_shadows, 1
set ray_trace_mode, 1     # Adds a clean silhouette outline around cartoons/sticks
set ray_trace_gain, 0.5   # Softens the outlines
set depth_cue, 1          # Enables depth shading (fog) to emphasize foreground loops
set specular, 0.3         # Controls surface shininess
set ray_opaque_background, 1
```

### 3.4. Canonical Camera Angle
A consistent camera matrix is applied via `set_view` before saving the image to ensure all snapshots are in the exact same spatial frame:
```text
set_view (-0.597478449, -0.617595792, -0.511463583, -0.502994895, 0.785388291, -0.360776752, 0.624511778, 0.041705273, -0.779902637, 0.000050262, -0.000109358, -182.850265503, 31.930021286, 32.076828003, 34.52369679, 141.553222656, 224.141143799, -20.000000000)
```

---

## 4. Programmatic Post-Processing
After rendering the raw PNG panels, `ColabMDA` uses `Pillow` (PIL) to:
1. Crop outer white space margins dynamically to auto-zoom on the pocket.
2. Arrange the panels into a 3 rows × 8 columns grid composite.
3. Add row labels ("A", "B", "C") and column time labels ("t = 0.00 ns", etc.).
4. Add a vertical gray dashed line to the G12C row indicating the transition threshold between the closed state (0.00–0.99 ns) and the open state (1.00–1.20 ns).
