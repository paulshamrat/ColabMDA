"""Downstream figures and structure views; simulation code stays separate."""

import os
import shutil
import subprocess
import sys
from pathlib import Path


def _headless_pyplot():
    """Override inherited notebook backends before Matplotlib validates them."""
    os.environ["MPLBACKEND"] = "Agg"
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def _publication_style():
    plt = _headless_pyplot()
    from cycler import cycler

    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.titlesize": 9,
            "axes.labelsize": 8.5,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "legend.fontsize": 7.5,
            "axes.prop_cycle": cycler(
                color=["#0072B2", "#D55E00", "#CC79A7", "#000000", "#E69F00", "#56B4E9"]
            )
            + cycler(linestyle=["-", "--", "-.", ":", "-", "--"]),
        }
    )


def openmm_compare(series_list, outdir):
    plt = _headless_pyplot()
    import pandas as pd

    outpath = Path(outdir) / "figures"
    outpath.mkdir(parents=True, exist_ok=True)

    _publication_style()

    metrics = {
        "rmsd.csv": ("Time (ps)", "RMSD (Å)", "System Stability (RMSD)"),
        "rg.csv": ("Time (ps)", "Radius of Gyration (Å)", "Compactness (Rg)"),
        "rmsf.csv": ("Residue Index", "RMSF (Å)", "Flexibility (RMSF)"),
    }

    def aggregate_system(dirs, metric_file):
        dfs = []
        for d in dirs:
            p = Path(d) / metric_file
            if p.exists():
                dfs.append(pd.read_csv(p))
        if not dfs:
            return None, None, None
        combined = pd.concat(dfs)
        col_name = dfs[0].columns[1]
        grouped = combined.groupby(combined.iloc[:, 0])
        mean = grouped[col_name].mean()
        std = grouped[col_name].std()
        return mean.index, mean.values, std.values

    for filename, (xlabel, ylabel, title) in metrics.items():
        plt.figure(figsize=(7, 3))
        found_any = False
        for item in series_list:
            if "=" not in item:
                continue
            label, dirs_str = item.split("=", 1)
            dirs = [d.strip() for d in dirs_str.split(",")]
            x, mean, std = aggregate_system(dirs, filename)
            if x is None:
                continue
            found_any = True
            p = plt.plot(x, mean, label=f"{label} (avg)", lw=2)
            color = p[0].get_color()
            plt.fill_between(x, mean - std, mean + std, color=color, alpha=0.2)

        if found_any:
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.title(title)
            plt.legend()
            plt.tight_layout()
            plt.savefig(outpath / filename.replace(".csv", "_avg.png"), dpi=300)
        plt.close()
    print(f"Aggregate plots saved in {outpath}")


def openmm_view(pdb_dir: str | None, topology: str | None, trajectory: str | None, resi: int = 12):
    from pathlib import Path

    # 1. Resolve simulation directory
    sim_dir = Path(pdb_dir or os.getcwd()).resolve()
    if not sim_dir.is_dir():
        sys.exit(f"Error: directory '{sim_dir}' does not exist.")

    # 2. Resolve topology and trajectory filenames
    topo_name = topology
    traj_name = trajectory

    if not topo_name:
        for candidate in ["prod_full.pdb", "kras_protein.pdb", "solvated.pdb"]:
            if (sim_dir / candidate).is_file():
                topo_name = candidate
                break
        if not topo_name:
            sys.exit(f"Error: No topology PDB file found in '{sim_dir}'. Run merge first.")

    if not traj_name:
        for candidate in ["prod_full.dcd", "kras_protein.dcd"]:
            if (sim_dir / candidate).is_file():
                traj_name = candidate
                break
        if not traj_name:
            sys.exit(f"Error: No trajectory DCD file found in '{sim_dir}'. Run merge first.")

    topo_path = sim_dir / topo_name
    traj_path = sim_dir / traj_name

    if not topo_path.is_file():
        sys.exit(f"Error: topology file not found: {topo_path}")
    if not traj_path.is_file():
        sys.exit(f"Error: trajectory file not found: {traj_path}")

    # Find crystal reference PDB (e.g., 4ldj_wt.pdb) to align coordinate space
    ref_name = f"{sim_dir.parent.name}.pdb"
    ref_path = sim_dir / ref_name
    if not ref_path.is_file():
        for p in sim_dir.glob("*.pdb"):
            if (
                p.name not in [topo_name, "solvated.pdb", "equilibrated.pdb"]
                and "cleaned" not in p.name
                and "protein" not in p.name
            ):
                ref_name = p.name
                ref_path = p
                break

    align_cmds = ""
    if ref_path.is_file():
        align_cmds = f"""
# 2b. Align trajectory to the crystal reference to ensure identical viewing orientation
load {ref_name}, crystal_ref
align kras and name CA, crystal_ref and name CA, mobile_state=1, target_state=1
delete crystal_ref
"""

    # 3. Create visualize.pml file
    pml_content = f"""# PyMOL Script generated by ColabMDA
# 1. Clear out all the old overlapping objects from memory
reinitialize
bg_color white

# 2. Load the trajectory files
load {topo_name}, kras
load_traj {traj_name}, kras
{align_cmds}
# 3. Hide solvent water and ions
hide everything, all
hide nonbonded, all
hide nb_spheres, all

# 4. Freeze the backbone tumbling rotation frame-by-frame
intra_fit kras and name CA

# 5. Generate your crisp secondary structure ribbon
dss kras
cartoon automatic, kras
show cartoon, kras

# 6. Apply your smooth cyan color scheme
color cyan, kras and name C*
util.cnc("kras")

# 7. Highlight the mutated residue side chain sticks
show sticks, kras and resi {resi} and not name N+C+O+H
color yellow, kras and resi {resi} and name SG
set stick_radius, 0.25

# 8. Focus camera right onto the protein
zoom kras and polymer, buffer=4
"""
    pml_file = sim_dir / "visualize.pml"
    try:
        pml_file.write_text(pml_content)
        print(f"→ Generated PyMOL script: {pml_file}")
    except Exception as e:
        sys.exit(f"Error writing PyMOL script: {e}")

    # 4. Run PyMOL
    if not shutil.which("pymol"):
        print("[WARN] PyMOL executable ('pymol') not found in PATH.")
        print("Please install PyMOL or run manually using the generated visualize.pml file:")
        print(f"  cd {sim_dir} && pymol visualize.pml")
        return

    print(f"→ Launching PyMOL: pymol visualize.pml in {sim_dir}...")
    try:
        subprocess.run(["pymol", "visualize.pml"], cwd=str(sim_dir), check=True)
    except Exception as e:
        print(f"Error running PyMOL: {e}")


def openmm_snapshots(
    config_path: str | None = None,
    output: str = "data/analysis/snapshots/master_3x8_snapshots.png",
    temp_dir: str = "scratch/master_grid_3x8",
    font_path: str | None = None,
):
    import json
    from pathlib import Path

    # 1. Load configuration
    default_config = {
        "camera_view": [
            -0.597478449,
            -0.617595792,
            -0.511463583,
            -0.502994895,
            0.785388291,
            -0.360776752,
            0.624511778,
            0.041705273,
            -0.779902637,
            0.000050262,
            -0.000109358,
            -182.850265503,
            31.930021286,
            32.076828003,
            34.52369679,
            141.553222656,
            224.141143799,
            -20.000000000,
        ],
        "align_ref_pdb": "analysis/clean_trajectories/g12c_protein.pdb",
        "stable_core_sel": "name CA and (resi 1-24 or resi 41-54 or resi 81-166)",
        "systems": {
            "WT": {
                "pdb": "analysis/clean_trajectories/wt_protein.pdb",
                "dcd": "analysis/clean_trajectories/wt_protein.dcd",
                "states": [1, 51, 91, 100, 101, 106, 111, 121],
                "times": [
                    "t = 0.00 ns",
                    "t = 0.50 ns",
                    "t = 0.90 ns",
                    "t = 0.99 ns",
                    "t = 1.00 ns",
                    "t = 1.05 ns",
                    "t = 1.10 ns",
                    "t = 1.20 ns",
                ],
                "mut_residue": 12,
            },
            "G12C": {
                "pdb": "analysis/clean_trajectories/g12c_protein.pdb",
                "dcd": "analysis/clean_trajectories/g12c_protein.dcd",
                "states": [1, 51, 91, 100, 101, 106, 111, 121],
                "times": [
                    "t = 0.00 ns",
                    "t = 0.50 ns",
                    "t = 0.90 ns",
                    "t = 0.99 ns",
                    "t = 1.00 ns",
                    "t = 1.05 ns",
                    "t = 1.10 ns",
                    "t = 1.20 ns",
                ],
                "mut_residue": 12,
            },
            "G12D": {
                "pdb": "analysis/clean_trajectories/g12d_protein.pdb",
                "dcd": "analysis/clean_trajectories/g12d_protein.dcd",
                "states": [1, 51, 91, 100, 101, 106, 111, 121],
                "times": [
                    "t = 0.00 ns",
                    "t = 0.50 ns",
                    "t = 0.90 ns",
                    "t = 0.99 ns",
                    "t = 1.00 ns",
                    "t = 1.05 ns",
                    "t = 1.10 ns",
                    "t = 1.20 ns",
                ],
                "mut_residue": 12,
            },
        },
    }

    if config_path:
        print(f"Loading custom configuration from {config_path}...")
        with open(config_path) as f:
            config = json.load(f)
    else:
        print("Using default KRAS 3x8 snapshot configuration...")
        config = default_config

    camera_view = config["camera_view"]
    align_ref_pdb = config.get("align_ref_pdb")
    stable_core_sel = config.get(
        "stable_core_sel", "name CA and (resi 1-24 or resi 41-54 or resi 81-166)"
    )
    systems = config["systems"]

    # 2. Check for PyMOL
    try:
        import pymol
        from pymol import cmd
    except ImportError:
        print("[ERROR] PyMOL python package is not installed in the active environment.")
        print("Please run this command using an environment with PyMOL (e.g. 'pymol-viz').")
        sys.exit(1)

    # Start PyMOL headlessly
    pymol.finish_launching(["pymol", "-qc"])

    # Create temp directory
    temp_path = Path(temp_dir)
    temp_path.mkdir(parents=True, exist_ok=True)

    rendered_files = {}

    for sys_name, sys_cfg in systems.items():
        pdb_path = sys_cfg["pdb"]
        dcd_path = sys_cfg["dcd"]
        states = sys_cfg["states"]
        mut_residue = sys_cfg.get("mut_residue", 12)

        if not os.path.exists(pdb_path):
            sys.exit(f"Error: PDB file not found: {pdb_path}")
        if not os.path.exists(dcd_path):
            sys.exit(f"Error: DCD file not found: {dcd_path}")

        print(f"Rendering snapshots for {sys_name}...")
        cmd.reinitialize()
        cmd.bg_color("white")

        # Load reference PDB for alignment if specified
        if align_ref_pdb and os.path.exists(align_ref_pdb):
            cmd.load(align_ref_pdb, "ref_align")

        # Load system
        cmd.load(pdb_path, "prot")
        cmd.load_traj(dcd_path, "prot")

        # Hide nonbonded
        cmd.hide("everything", "all")
        cmd.hide("nonbonded", "all")
        cmd.hide("nb_spheres", "all")

        # Step 1: Align all states of prot to state 1 of prot using stable core
        cmd.intra_fit(f"prot and {stable_core_sel}", 1)

        # Step 2: Align state 1 of prot to ref_align using stable core and apply to all states
        if align_ref_pdb and os.path.exists(align_ref_pdb):
            cmd.align(
                f"prot and {stable_core_sel}",
                f"ref_align and {stable_core_sel}",
                mobile_state=1,
                target_state=1,
            )
            cmd.delete("ref_align")

        cmd.dss("prot")
        cmd.cartoon("automatic", "prot")
        cmd.cartoon("loop", "prot and resid 57-75")
        cmd.show("cartoon", "prot")

        # Colors
        cmd.color("gray90", "prot")
        cmd.color("orange", "prot and resid 57-75")  # Switch II
        cmd.color("yellow", f"prot and resid {mut_residue}")

        # Mutation sticks
        cmd.show("sticks", f"prot and resi {mut_residue} and not name N+C+O+H")
        cmd.set("stick_radius", 0.25)

        # Transparency
        cmd.set("cartoon_transparency", 0.6, f"prot and not (resid 57-75 or resid {mut_residue})")
        cmd.set("cartoon_transparency", 0.0, f"resid 57-75 or resid {mut_residue}")

        # Ray tracing options
        cmd.set("ray_shadows", 1)
        cmd.set("ray_trace_mode", 1)
        cmd.set("ray_trace_gain", 0.5)
        cmd.set("depth_cue", 1)
        cmd.set("specular", 0.3)
        cmd.set("ray_opaque_background", 1)

        cmd.set_view(camera_view)

        rendered_files[sys_name] = []
        for state in states:
            p_out = temp_path / f"{sys_name}_state_{state}.png"
            rendered_files[sys_name].append(str(p_out))
            cmd.set("state", state)
            cmd.png(str(p_out), width=800, height=800, ray=1)

    pymol.cmd.quit()

    print("All snapshots rendered! Compiling image grid...")

    import PIL.ImageOps
    from PIL import Image, ImageDraw, ImageFont

    all_images = []
    for sys_name in systems.keys():
        for path in rendered_files[sys_name]:
            all_images.append(Image.open(path))

    # Bounding boxes for cropping
    bboxes = []
    for img in all_images:
        gray = img.convert("L")
        inverted = PIL.ImageOps.invert(gray)
        bbox = inverted.getbbox()
        if bbox:
            bboxes.append(bbox)

    orig_width, orig_height = all_images[0].size
    if bboxes:
        left_crop = min(b[0] for b in bboxes)
        top_crop = min(b[1] for b in bboxes)
        right_crop = max(b[2] for b in bboxes)
        bottom_crop = max(b[3] for b in bboxes)

        padding = 20
        left_crop = max(0, left_crop - padding)
        top_crop = max(0, top_crop - padding)
        right_crop = min(orig_width, right_crop + padding)
        bottom_crop = min(orig_height, bottom_crop + padding)

        cropped_images = [
            img.crop((left_crop, top_crop, right_crop, bottom_crop)) for img in all_images
        ]
    else:
        cropped_images = all_images

    # Grid configuration (N columns, M rows total)
    num_rows = len(systems)
    num_cols = len(list(systems.values())[0]["states"])

    panel_w, panel_h = cropped_images[0].size
    spacing_x = 12
    spacing_y = 15
    label_y_area = 120

    total_width = panel_w * num_cols + spacing_x * (num_cols - 1)
    total_height = (panel_h + label_y_area) * num_rows + spacing_y * (num_rows - 1)

    combined = Image.new("RGB", (total_width, total_height), "white")
    draw = ImageDraw.Draw(combined)

    # Font resolution
    font_path_to_use = font_path
    if not font_path_to_use:
        for path in [
            "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
            "/usr/share/fonts/truetype/freefont/FreeSansBold.ttf",
            "/usr/share/fonts/liberation/LiberationSans-Bold.ttf",
        ]:
            if os.path.exists(path):
                font_path_to_use = path
                break

    try:
        if font_path_to_use:
            font_panel_letter = ImageFont.truetype(font_path_to_use, 64)
            # Find corresponding regular font for time labels
            regular_font_path = font_path_to_use.replace("-Bold", "").replace("Bold", "")
            if not os.path.exists(regular_font_path):
                regular_font_path = font_path_to_use
            font_time = ImageFont.truetype(regular_font_path, 72)
        else:
            font_panel_letter = ImageFont.load_default()
            font_time = ImageFont.load_default()
    except OSError:
        font_panel_letter = ImageFont.load_default()
        font_time = ImageFont.load_default()

    img_idx = 0
    for r, (sys_name, sys_cfg) in enumerate(systems.items()):
        row_y = r * (panel_h + label_y_area + spacing_y)
        for c in range(num_cols):
            img = cropped_images[img_idx]
            x = c * (panel_w + spacing_x)
            combined.paste(img, (x, row_y))

            # Panel letter on the first panel of each row
            if c == 0:
                letter = chr(65 + r)  # A, B, C...
                draw.text((x + 20, row_y + 20), letter, fill="black", font=font_panel_letter)

            # Time label
            time_str = sys_cfg["times"][c]
            try:
                time_w = draw.textlength(time_str, font=font_time)
            except AttributeError:
                time_w, _ = draw.textsize(time_str, font=font_time)

            time_x = x + (panel_w - time_w) // 2
            time_y = row_y + panel_h + 20
            draw.text((time_x, time_y), time_str, fill="black", font=font_time)

            img_idx += 1

    # Draw faint dashed separator line only for G12C row (row r = 1) between columns 3 (0.99ns) and 4 (1.00ns)
    # If the system name G12C exists as the second row (standard KRAS layout)
    sys_list = list(systems.keys())
    if "G12C" in sys_list and sys_list.index("G12C") == 1 and num_cols >= 8:
        x_divider = 4 * panel_w + 3 * spacing_x + spacing_x // 2
        y_start = 1 * (panel_h + label_y_area + spacing_y)
        y_end = y_start + panel_h
        dash_len = 15
        gap_len = 10
        curr_y = y_start
        while curr_y < y_end:
            draw.line(
                [(x_divider, curr_y), (x_divider, min(curr_y + dash_len, y_end))],
                fill=(180, 180, 180),
                width=4,
            )
            curr_y += dash_len + gap_len

    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    combined.save(str(output_path), dpi=(300, 300))
    print(f"Master snapshots grid successfully saved to {output_path}!")

    # Clean up temp files
    for sys_name in systems.keys():
        for path in rendered_files[sys_name]:
            if os.path.exists(path):
                os.remove(path)
    if temp_path.exists():
        try:
            temp_path.rmdir()
        except OSError:
            pass


def plot_analysis(outdir):
    """Render saved measurements without recalculating trajectory metrics."""
    plt = _headless_pyplot()
    import pandas as pd

    root = Path(outdir)
    _publication_style()
    names = {
        "rmsd": ("RMSD (Å)", "rmsd_vs_time.png"),
        "rg": ("Rg (Å)", "rg_vs_time.png"),
        "rmsf": ("RMSF (Å)", "rmsf_per_residue.png"),
    }
    for key, (ylabel, filename) in names.items():
        for table in root.rglob(f"{key}.csv"):
            data = pd.read_csv(table)
            fig, ax = plt.subplots(figsize=(7, 3), layout="constrained")
            ax.plot(data.iloc[:, 0], data.iloc[:, 1], color="#0072B2")
            ax.set(xlabel=data.columns[0], ylabel=ylabel)
            target = table.parent / "figures"
            target.mkdir(exist_ok=True)
            fig.savefig(target / filename, dpi=300)
            plt.close(fig)
        table = root / "aggregate" / f"{key}_avg.csv"
        if table.exists():
            data = pd.read_csv(table)
            fig, ax = plt.subplots(figsize=(7, 3), layout="constrained")
            ax.plot(data.iloc[:, 0], data["Mean"], color="#0072B2")
            ax.fill_between(
                data.iloc[:, 0],
                data["Mean"] - data["Std"],
                data["Mean"] + data["Std"],
                color="#0072B2",
                alpha=0.2,
            )
            ax.set(xlabel=data.columns[0], ylabel=ylabel)
            target = table.parent / "figures"
            target.mkdir(exist_ok=True)
            fig.savefig(target / f"{key}_aggregate.png", dpi=300)
            plt.close(fig)


def plot_equilibration(workdir):
    """Plot saved thermodynamic logs after the numerical QC stage."""
    import json

    plt = _headless_pyplot()
    import pandas as pd

    _publication_style()
    root = Path(workdir)
    if not (root / "equilibration_qc.json").exists():
        return
    qc = json.loads((root / "equilibration_qc.json").read_text())
    fig, axes = plt.subplots(2, 2, figsize=(7, 4.5), layout="constrained")
    for ax, (filename, column) in zip(
        axes.flat,
        [
            ("nvt.log", "Temperature (K)"),
            ("npt.log", "Density (g/mL)"),
            ("nvt.log", "Potential Energy (kJ/mole)"),
            ("npt.log", "Potential Energy (kJ/mole)"),
        ],
    ):
        data = pd.read_csv(root / filename)
        ax.plot(data["Time (ps)"], data[column], color="#0072B2")
        ax.set(xlabel="Time (ps)", ylabel=column)
    fig.suptitle(
        "Equilibration QC: " + ("passed" if qc["passed"] else "review required"), fontsize=9
    )
    target = root / "figures"
    target.mkdir(exist_ok=True)
    fig.savefig(target / "equilibration_qc.png", dpi=300)
    plt.close(fig)
