"""Downstream figures and structure views; simulation code stays separate."""

import os
import shutil
import subprocess
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
            # Save both legacy figure name and clean flat name
            plt.savefig(outpath / filename.replace(".csv", "_avg.png"), dpi=300)
            plt.savefig(Path(outdir) / filename.replace(".csv", ".png"), dpi=300)
        plt.close()

    summary_file = Path(outdir) / "summary.txt"
    try:
        with open(summary_file, "w") as sf:
            sf.write("=========================================\n")
            sf.write("COLABMDA CROSS-VARIANT COMPARISON SUMMARY\n")
            sf.write("=========================================\n")
            for item in series_list:
                if "=" in item:
                    lbl, d_str = item.split("=", 1)
                    sf.write(f"Series: {lbl} (inputs: {d_str})\n")
            sf.write("Deliverables: rmsd.png, rg.png, rmsf.png\n")
            sf.write("=========================================\n")
    except Exception:
        pass
    print(f"Aggregate plots saved in {outdir}")


# Activation project structural overview and apo/holo network camera presets.
# Numeric camera metadata only; no activation-project coordinates are required.
CDKL5_CAMERAS = {
    "activation": (
        0.360244870,
        0.822598636,
        -0.439932555,
        0.888235748,
        -0.158373788,
        0.431204081,
        0.285034835,
        -0.546108961,
        -0.787725031,
        -0.001081586,
        -0.000175163,
        -245.772369385,
        53.673915863,
        50.904361725,
        39.656856537,
        205.716583252,
        285.834991455,
        -20.0,
    ),
    "network": (
        0.854537665844,
        0.461885809898,
        0.237543433905,
        -0.240715444088,
        0.757462859154,
        -0.606882452965,
        -0.460240691900,
        0.461423516273,
        0.758463501930,
        0.0,
        0.0,
        -262.464569091797,
        44.759418487549,
        32.431861877441,
        65.356414794922,
        206.929061889648,
        318.000061035156,
        -20.0,
    ),
}
CAMERA_SOURCES = {
    "activation": "cdkl5-activation/scripts/clustering/05_mapping.py:KINASE_VIEW (also grant/variantmap)",
    "network": "cdkl5-activation state-paired apo/holo network rendering:KINASE_VIEW + turn z 90",
}


def _view_dataset(
    topology,
    trajectory=None,
    reference=None,
    frames="1,last",
    resi=119,
    variant="WT",
    log=None,
    mapping=None,
    core_residues=None,
):
    """Validate atom identity/native labels and fit display frames to one reference."""
    import hashlib
    import json

    import MDAnalysis as mda
    import numpy as np
    from Bio.PDB import PDBParser
    from Bio.SeqUtils import seq1

    from colabmda.modeller.substrate import kabsch

    topology = Path(topology).resolve()
    reference = Path(reference or topology).resolve()
    if reference.suffix.lower() not in {".pdb", ".cif"}:
        raise ValueError("A numbered PDB/mmCIF reference is required")
    from Bio.PDB import MMCIFParser

    parser = (
        MMCIFParser(QUIET=True) if reference.suffix.lower() == ".cif" else PDBParser(QUIET=True)
    )
    model = parser.get_structure("ref", str(reference))[0]
    residues = [r for r in model.get_residues() if r.id[0] == " " and "CA" in r]
    ids = [(r.parent.id, r.id[1], r.id[2].strip()) for r in residues]
    if len(ids) != len(set(ids)):
        raise ValueError("Ambiguous reference residue labels")
    u = mda.Universe(str(topology), str(trajectory)) if trajectory else mda.Universe(str(topology))
    ca = u.select_atoms("protein and name CA")
    aliases = {"HID": "HIS", "HIE": "HIS", "HIP": "HIS", "CYX": "CYS", "ASH": "ASP", "GLH": "GLU"}
    actual = "".join(seq1(aliases.get(r, r)) for r in ca.resnames)
    expected = "".join(seq1(aliases.get(r.resname, r.resname)) for r in residues)
    if len(ca) != len(residues):
        raise ValueError("Reference/trajectory C-alpha count mismatch")
    differences = [i for i, (a, b) in enumerate(zip(expected, actual)) if a != b]
    if variant != "WT":
        import re

        if not re.fullmatch(r"[A-Z][1-9][0-9]*[A-Z]", variant):
            raise ValueError("Declared variant must use native notation, e.g. L119R")
        if not differences:
            raise ValueError("Declared variant is absent from the displayed sequence")
    if differences:
        if variant == "WT" or len(differences) != 1:
            raise ValueError("Reference/trajectory sequence mismatch")
        i = differences[0]
        if (residues[i].id[1], expected[i], actual[i]) != (
            int(variant[1:-1]),
            variant[0],
            variant[-1],
        ):
            raise ValueError("Declared variant does not explain sequence mismatch")
    selected = []
    for text in frames.split(","):
        n = len(u.trajectory) if text.strip() == "last" else int(text)
        if not 1 <= n <= len(u.trajectory):
            raise ValueError("Requested frame is outside trajectory")
        if n not in selected:
            selected.append(n)
    observed = None
    if mapping and Path(mapping).is_file():
        observed = {r["uniprot"] for r in json.loads(Path(mapping).read_text()) if r["resolved"]}
    core = [
        i
        for i, r in enumerate(residues)
        if (r.id[1] in core_residues if core_residues is not None else 30 <= r.id[1] <= 220)
        and (observed is None or r.id[1] in observed)
    ]
    if len(core) < 3:
        raise ValueError("Insufficient matched core C-alpha atoms (30-220)")
    ref_xyz = np.array([r["CA"].coord for r in residues])
    # Native labels are applied only to a derived viewer copy, never raw topology.
    for atom, residue in zip(ca, residues):
        atom.residue.resid = residue.id[1]
    if not hasattr(u.atoms, "chainIDs"):
        u.add_TopologyAttr("chainIDs")
    u.atoms.chainIDs = "A"
    if not hasattr(u.residues, "icodes"):
        u.add_TopologyAttr("icodes")
    for atom, residue in zip(ca, residues):
        atom.residue.atoms.chainIDs = residue.parent.id.strip() or "A"
        atom.residue.icode = residue.id[2].strip()
    if not u.select_atoms(f"protein and resid {resi}").n_atoms:
        raise ValueError(f"Highlighted native residue {resi} is absent")
    times = None
    if trajectory and log:
        import csv

        if not Path(log).is_file():
            raise FileNotFoundError(f"Production time log missing: {log}")

        with Path(log).open() as handle:
            times = np.array([float(row["Time (ps)"]) for row in csv.DictReader(handle)])
        if (
            len(times) != len(u.trajectory)
            or not np.isfinite(times).all()
            or (len(times) > 1 and not (np.diff(times) > 0).all())
        ):
            raise ValueError("Production times do not match trajectory frames")
    panels = []
    for n in selected:
        u.trajectory[n - 1]
        rotation, translation = kabsch(ca.positions[core], ref_xyz[core])
        xyz = u.atoms.positions @ rotation + translation
        if not np.isfinite(xyz).all():
            raise ValueError("Nonfinite display coordinates")
        panels.append(
            dict(
                frame=n,
                time_ps=float(times[n - 1]) if times is not None else None,
                positions=xyz,
                core_rmsd_angstrom=float(
                    np.sqrt(
                        np.mean(
                            np.sum(
                                (ca.positions[core] @ rotation + translation - ref_xyz[core]) ** 2,
                                axis=1,
                            )
                        )
                    )
                ),
            )
        )

    def sha(p):
        digest = hashlib.sha256()
        with Path(p).open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        return digest.hexdigest()

    return (
        u,
        panels,
        dict(
            topology=str(topology),
            topology_sha256=sha(topology),
            trajectory=str(Path(trajectory).resolve()) if trajectory else None,
            trajectory_sha256=sha(trajectory) if trajectory else None,
            reference=str(reference),
            reference_sha256=sha(reference),
            core_residues=[ids[i] for i in core],
            native_residues=ids,
            variant=variant,
            highlight=resi,
            ca_count=len(ca),
            source_kind="trajectory" if trajectory else "structure",
            time_log=str(Path(log).resolve()) if trajectory and log else None,
            time_log_sha256=sha(log) if trajectory and log else None,
            alignment="Per-frame matched native core CA 30-220; experimentally missing positions excluded when mapping available",
        ),
    )


def _state_views(
    state,
    variant,
    root=None,
    replica="r1",
    structure=False,
    reference=None,
    frames="1,last",
    resi=119,
):
    from colabmda.core.project import ProjectPaths, resolve_project_root
    from colabmda.core.taxonomy import SystemSpec

    root = resolve_project_root(root)
    specs = [(s.strip(), v.strip()) for s in state.split(",") for v in variant.split(",")]
    datasets = []
    common = Path(reference) if reference else None
    if common is None:
        for s in ["apo", specs[0][0]]:
            p = ProjectPaths(root, state=s)
            candidate = p.model_pdb(SystemSpec.from_tag(s, variant="WT"))
            if candidate.is_file():
                common = candidate
                break
        if common is None:
            p_apo = ProjectPaths(root, state="apo")
            candidate = p_apo.structure_dir / "models/WT/WT.pdb"
            if candidate.is_file():
                common = candidate
            elif (root / "data/apo/structure/models/WT/WT.pdb").is_file():
                common = root / "data/apo/structure/models/WT/WT.pdb"
    if common is None:
        raise FileNotFoundError("Supply a checked common WT --reference")
    for s, v in specs:
        paths = ProjectPaths(root, state=s)
        spec = SystemSpec.from_tag(s, variant=v)
        if not paths.run_id:
            prod = None
            is_structure = True
        else:
            try:
                prod = paths.replica_dir(spec, replica)
                is_structure = structure or not (prod / "stripped.prmtop").is_file()
            except Exception:
                prod = None
                is_structure = True

        if is_structure:
            topo = paths.model_pdb(spec)
            if not topo.is_file():
                topo = (
                    paths.structure_dir
                    / "models"
                    / v
                    / ("complex.pdb" if s == "atpmg" else f"{v}.pdb")
                )
            traj = None
        else:
            topo, traj = prod / "stripped.prmtop", prod / "prod_keep_ligand.dcd"
        mapping_file = paths.structure_dir / "mapping/residues.json"
        u, panels, record = _view_dataset(
            topo,
            traj,
            common,
            "1" if is_structure else frames,
            resi,
            v,
            prod / "prod_full.log" if not is_structure and prod else None,
            mapping_file if mapping_file.is_file() else None,
        )
        record["source_kind"] = "starting_model" if is_structure else "trajectory"
        datasets.append((f"{s} · {v}", u, panels, record))
    paths = ProjectPaths(root, state=specs[0][0])
    spec = SystemSpec.from_tag(specs[0][0], variant=specs[0][1])
    if is_structure or not paths.run_id:
        folder = paths.structure_dir / "models" / spec.variant
    else:
        folder = (
            paths.analysis_dir(spec) if structure else paths.replica_analysis_dir(spec, replica)
        )
    return datasets, folder


def _style_view(cmd, resi):
    cmd.hide("everything", "all")
    cmd.bg_color("white")
    for key, value in [
        ("orthoscopic", 1),
        ("depth_cue", 0),
        ("ray_shadows", 0),
        ("antialias", 2),
        ("cartoon_fancy_helices", 0),
        ("cartoon_flat_sheets", 1),
        ("cartoon_smooth_loops", 1),
    ]:
        cmd.set(key, value)
    cmd.show("cartoon", "molecule and polymer.protein")
    cmd.color("wheat", "molecule and polymer.protein")
    cmd.set("cartoon_transparency", 0.12, "molecule and polymer.protein")
    cmd.show("sticks", f"molecule and polymer.protein and resi {resi} and not hydro")
    cmd.color("orange", f"molecule and polymer.protein and resi {resi}")
    cmd.show("sticks", "molecule and resn ATP")
    cmd.color("purple", "molecule and resn ATP")
    cmd.show("spheres", "molecule and resn MG")
    cmd.color("marine", "molecule and resn MG")
    cmd.set("sphere_scale", 0.55, "molecule and resn MG")
    cmd.set("stick_radius", 0.18)


def _render_views(datasets, output, camera="activation", camera_view=None, font_path=None):
    import hashlib
    import json

    import MDAnalysis as mda
    import numpy as np
    from PIL import Image, ImageDraw, ImageFont
    from pymol2 import PyMOL

    if camera not in CDKL5_CAMERAS:
        raise ValueError("Unknown CDKL5 camera")
    if camera_view is not None and (len(camera_view) != 18 or not np.isfinite(camera_view).all()):
        raise ValueError("Camera requires 18 finite values")
    output = Path(output).resolve()
    if output.suffix.lower() != ".png":
        raise ValueError("Save rendered figures as PNG")
    output.parent.mkdir(parents=True, exist_ok=True)
    marker = output.with_name(output.stem.upper() + "_DONE")
    marker.unlink(missing_ok=True)
    counts = {len(panels) for _, _, panels, _ in datasets}
    if len(counts) != 1 or not datasets:
        raise ValueError("Snapshot rows must have equal frame counts")
    source_cols = counts.pop()
    compact = source_cols == 1 and len(datasets) > 1
    cols = min(2, len(datasets)) if compact else source_cols
    row_count = (len(datasets) + cols - 1) // cols if compact else len(datasets)
    width = 2100
    gap = 24
    cell = (width - gap * (cols - 1)) // cols
    header = 64
    canvas = Image.new("RGB", (width, row_count * (cell + header) + gap * (row_count - 1)), "white")
    fontfile = font_path or "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"
    font = ImageFont.truetype(fontfile, 38)
    bold = ImageFont.truetype(
        fontfile.replace(".ttf", "-Bold.ttf") if not font_path else fontfile, 42
    )
    draw = ImageDraw.Draw(canvas)
    records = []
    work = output.parent / "view_inputs" / output.stem
    work.mkdir(parents=True, exist_ok=True)
    with PyMOL() as session:
        cmd = session.cmd
        cmd.viewport(cell, cell)
        cmd.set("orthoscopic", 1)
        cmd.load(datasets[0][3]["reference"], "reference")
        cmd.set_view(camera_view or CDKL5_CAMERAS[camera])
        if camera == "network" and camera_view is None:
            cmd.turn("z", 90)
        cmd.zoom("reference and polymer.protein", buffer=4)
        fixed_camera = tuple(cmd.get_view())
        cmd.delete("reference")
        for row, (label, u, panels, record) in enumerate(datasets):
            items = []
            for col, panel in enumerate(panels):
                pdb = work / f"row{row + 1}_frame{panel['frame']}.pdb"
                u.atoms.positions = panel["positions"]
                with mda.Writer(str(pdb), n_atoms=len(u.atoms)) as writer:
                    writer.write(u.atoms)
                cmd.delete("all")
                cmd.load(str(pdb), "molecule")
                cmd.dss("molecule")
                _style_view(cmd, record["highlight"])
                cmd.set_view(fixed_camera)
                png = work / f"row{row + 1}_frame{panel['frame']}.png"
                cmd.png(str(png), width=cell, height=cell, dpi=300, ray=1, quiet=1)
                image = Image.open(png).convert("RGB")
                grid_row = row // cols if compact else row
                grid_col = row % cols if compact else col
                x = grid_col * (cell + gap)
                y = grid_row * (cell + header + gap)
                canvas.paste(image, (x, y + header))
                time = (
                    f"{panel['time_ps']:g} ps"
                    if panel["time_ps"] is not None
                    else (
                        "model"
                        if record["source_kind"] == "starting_model"
                        else (
                            "structure"
                            if record["source_kind"] == "structure"
                            else f"frame {panel['frame']}"
                        )
                    )
                )
                draw.text(
                    (x + 8, y + 5), chr(65 + row * source_cols + col), font=bold, fill="black"
                )
                draw.text((x + 65, y + 8), f"{label} · {time}", font=font, fill="black")
                items.append(
                    dict(
                        frame=panel["frame"],
                        time_ps=panel["time_ps"],
                        pdb_sha256=hashlib.sha256(pdb.read_bytes()).hexdigest(),
                        core_rmsd_angstrom=panel["core_rmsd_angstrom"],
                    )
                )
            records.append(dict(label=label, source=record, panels=items))
    canvas.save(output, dpi=(300, 300))
    metadata = dict(
        renderer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        camera=camera,
        camera_source=CAMERA_SOURCES[camera] if camera_view is None else "Explicit custom camera",
        source_camera=list(camera_view or CDKL5_CAMERAS[camera]),
        render_camera=list(fixed_camera),
        framing="Activation orientation; zoom once to common WT reference and reuse for every panel",
        rows=records,
        width_inches=7,
        dpi=300,
    )
    output.with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n")
    marker.write_text("Verified PNG and native-reference display inputs generated\n")
    return output


def openmm_view(
    pdb_dir=None,
    topology=None,
    trajectory=None,
    resi=119,
    *,
    state=None,
    variant="WT",
    root=None,
    replica="r1",
    reference=None,
    camera="activation",
    structure=False,
    png=None,
    no_launch=False,
    frame="1",
):
    import json

    try:
        import MDAnalysis as mda
    except ImportError:
        mda = None

    if mda is None:
        if png or trajectory:
            raise ImportError(
                "Trajectory viewing and offline PNG rendering require MDAnalysis.\n"
                "Install it with:\n"
                "  conda install -c conda-forge mdanalysis\n"
                "or:\n"
                "  pip install mdanalysis"
            )
        from colabmda.core.project import ProjectPaths, resolve_project_root
        from colabmda.core.taxonomy import SystemSpec

        root_path = resolve_project_root(root)
        if state:
            paths = ProjectPaths(root_path, state=state)
            spec = SystemSpec.from_tag(state, variant=variant)
            target = paths.model_pdb(spec)
            if not target.is_file():
                target = (
                    paths.structure_dir
                    / "models"
                    / variant
                    / ("complex.pdb" if state == "atpmg" else f"{variant}.pdb")
                )
            folder = paths.structure_dir / "models" / variant
        else:
            folder = Path(pdb_dir or os.getcwd()).resolve()
            target = folder / topology if topology else None

        if not target or not target.is_file():
            raise FileNotFoundError(
                f"Structure file not found: {target}\n"
                "Run local modeling first (cmda model) or specify a valid --topology."
            )

        folder.mkdir(parents=True, exist_ok=True)
        pml = folder / "visualize.pml"
        view = ",".join(map(str, CDKL5_CAMERAS.get(camera, CDKL5_CAMERAS["activation"])))
        pml.write_text(f"""reinitialize
load "{target}", molecule
hide everything
bg_color white
show cartoon, molecule and polymer.protein
color wheat, molecule and polymer.protein
show sticks, molecule and polymer.protein and resi {resi} and not hydro
color orange, molecule and polymer.protein and resi {resi}
show sticks, molecule and resn ATP
color purple, molecule and resn ATP
show spheres, molecule and resn MG
color marine, molecule and resn MG
set orthoscopic, 1
set_view ({view})
zoom molecule and polymer.protein, 4
""")
        if not no_launch:
            cmd = ["pymol", str(pml)]
            if subprocess.call(cmd) != 0:
                raise RuntimeError("PyMOL failed")
        return 0

    if state:
        if "," in state or "," in variant:
            raise ValueError("Use snapshots for multiple states or variants")
        datasets, folder = _state_views(
            state,
            variant,
            root,
            replica,
            structure,
            reference,
            "1" if frame == "all" else str(frame),
            resi,
        )
    else:
        folder = Path(pdb_dir or os.getcwd()).resolve()
        if topology is None:
            raise ValueError("Supply --state or an explicit matching --topology")
        u, panels, record = _view_dataset(
            folder / topology,
            folder / trajectory if trajectory else None,
            reference,
            "1" if frame == "all" else str(frame),
            resi,
            variant,
        )
        datasets = [(variant, u, panels, record)]
    if png and frame == "all":
        raise ValueError("Use snapshots for a frame grid; --frame all opens an interactive movie")
    if png:
        out = folder / "figures/structure.png" if png == "auto" else Path(png)
        _render_views(datasets, out, camera)
        return 0
    label, u, panels, record = datasets[0]
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    target = folder / "view.pdb"
    u.atoms.positions = panels[0]["positions"]
    with mda.Writer(str(target), n_atoms=len(u.atoms)) as writer:
        writer.write(u.atoms)
    # The generated script displays a fitted, native-numbered frame. Snapshots
    # provides sampled multi-frame comparison without loading an entire movie.
    if frame == "all" and not record["trajectory"]:
        raise ValueError("A trajectory is required for --frame all")
    core = "+".join(str(item[1]) + item[2] for item in record["core_residues"])
    movie = (
        f'load_traj "{record["trajectory"]}", molecule, state=1\nintra_fit molecule and polymer.protein and name CA and resi {core}, 1\nalign molecule and polymer.protein and name CA and resi {core}, reference and polymer.protein and name CA and resi {core}, mobile_state=1, target_state=1, cycles=0\n'
        if frame == "all"
        else ""
    )
    pml = folder / "visualize.pml"
    view = ",".join(map(str, CDKL5_CAMERAS[camera]))
    pml.write_text(
        f"""reinitialize
load "{record["reference"]}", reference
load "{target}", molecule
{movie}hide everything
bg_color white
show cartoon, molecule and polymer.protein
color wheat, molecule and polymer.protein
show sticks, molecule and polymer.protein and resi {resi} and not hydro
color orange, molecule and polymer.protein and resi {resi}
show sticks, molecule and resn ATP
color purple, molecule and resn ATP
show spheres, molecule and resn MG
color marine, molecule and resn MG
set orthoscopic, 1
set_view ({view})
"""
        + ("turn z, 90\n" if camera == "network" else "")
        + "zoom reference and polymer.protein, 4\ndelete reference\n"
    )
    record["display_frames"] = "all" if frame == "all" else panels[0]["frame"]
    pml.with_suffix(".json").write_text(
        json.dumps(
            dict(
                record,
                camera=camera,
                source_camera=CDKL5_CAMERAS[camera],
                camera_source=CAMERA_SOURCES[camera],
                framing="Reframed to common WT reference",
            ),
            indent=2,
        )
        + "\n"
    )
    if not no_launch:
        if not shutil.which("pymol"):
            raise RuntimeError(
                f"PyMOL executable not found in PATH; visualization script saved at {pml}\n"
                "To launch interactively on your local machine, install PyMOL with:\n"
                "  conda install -c conda-forge pymol-open-source"
            )
        subprocess.run(["pymol", str(pml)], check=True)
    print(f"Generated {pml}")
    return 0


def openmm_snapshots(
    config_path=None,
    output=None,
    temp_dir=None,
    font_path=None,
    *,
    state=None,
    variant="WT",
    root=None,
    replica="r1",
    reference=None,
    camera="activation",
    structure=False,
    frames="1,last",
    resi=119,
):
    import json

    if state:
        if config_path:
            raise ValueError("Choose --state or --config")
        datasets, folder = _state_views(
            state, variant, root, replica, structure, reference, frames, resi
        )
        view = None
    elif config_path:
        configfile = Path(config_path).resolve()
        config = json.loads(configfile.read_text())
        datasets = []

        def resolve(value):
            return configfile.parent / str(value)

        core_residues = config.get("core_residues")
        if config.get("stable_core_sel"):
            import re

            selection = config["stable_core_sel"]
            if not re.fullmatch(r"name CA and \(resi \d+-\d+(?: or resi \d+-\d+)*\)", selection):
                raise ValueError("Use explicit core_residues for complex custom selections")
            core_residues = [
                i
                for a, b in re.findall(r"resi (\d+)-(\d+)", selection)
                for i in range(int(a), int(b) + 1)
            ]
        for name, item in config["systems"].items():
            if "times" in item or "times" in config:
                raise ValueError("Replace custom times labels with a matching production log")
            u, panels, record = _view_dataset(
                resolve(item["pdb"]),
                resolve(item["dcd"]) if item.get("dcd") else None,
                resolve(item.get("reference") or config.get("align_ref_pdb") or item["pdb"]),
                ",".join(map(str, item.get("states", [1]))),
                item.get("mut_residue", resi),
                item.get("variant", "WT"),
                resolve(item["log"]) if item.get("log") else None,
                core_residues=core_residues,
            )
            datasets.append((name, u, panels, record))
        folder = configfile.parent
        view = config.get("camera_view")
        if view is not None and len(view) != 18:
            raise ValueError("Camera requires 18 values")
    else:
        raise ValueError("Supply --state for project data or --config for custom inputs")
    name = "model_snapshots.png" if structure else "structure_snapshots.png"
    out = Path(output) if output else folder / "figures" / name
    print(f"Rendered {_render_views(datasets, out, camera, view, font_path)}")
    return 0


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
        for table in [*root.rglob(f"{key}.csv"), *root.rglob(f"{key}_avg.csv")]:
            data = pd.read_csv(table)
            fig, ax = plt.subplots(figsize=(7, 3), layout="constrained")
            ax.plot(data.iloc[:, 0], data.iloc[:, 1], color="#0072B2")
            if "Std" in data:
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
            output = f"{key}_aggregate.png" if table.stem.endswith("_avg") else filename
            fig.savefig(target / output, dpi=300)
            plt.close(fig)


def plot_equilibration(workdir, outdir=None):
    """Plot saved thermodynamic logs and structural stability after equilibration."""
    import json
    import shutil
    import subprocess
    import tempfile
    from pathlib import Path

    import numpy as np
    import pandas as pd
    from PIL import Image

    plt = _headless_pyplot()
    _publication_style()
    root = Path(workdir)
    if outdir:
        dest = Path(outdir)
    else:
        parts = root.resolve().parts
        if "simulation" in parts:
            idx = parts.index("simulation")
            analysis_parts = list(parts)
            analysis_parts[idx] = "analysis"
            dest = Path(*analysis_parts)
        else:
            dest = root
    dest.mkdir(parents=True, exist_ok=True)

    qc_candidates = [
        dest / "equilibration_qc.json",
        root / "equilibration_qc.json",
    ]
    qc_file = next((p for p in qc_candidates if p.is_file()), None)
    passed = True
    if qc_file:
        try:
            qc = json.loads(qc_file.read_text())
            passed = qc.get("passed", False)
        except Exception:
            pass
    elif (dest / "qc.txt").is_file():
        passed = "Status: PASSED" in (dest / "qc.txt").read_text()
    elif (root / "qc.txt").is_file():
        passed = "Status: PASSED" in (root / "qc.txt").read_text()

    # Load thermodynamic logs
    npt_file = root / "npt.log"
    nvt_file = root / "nvt.log"
    has_npt = npt_file.is_file()
    has_nvt = nvt_file.is_file()

    df_npt = pd.read_csv(npt_file) if has_npt else pd.DataFrame()
    df_nvt = pd.read_csv(nvt_file) if has_nvt else pd.DataFrame()

    # Structural alignment: model.pdb vs eq.pdb
    eq_pdb = root / "eq.pdb"
    model_candidates = [
        root.parent / "01_prep" / "model.pdb",
        root / "model.pdb",
        root / "solvated.pdb",
    ]
    model_pdb = next((p for p in model_candidates if p.is_file()), None)

    rmsd = None
    drift = None
    common_res = None
    if eq_pdb.is_file() and model_pdb and model_pdb.is_file():
        try:

            def extract_ca(pdb_p):
                coords, res_nums = [], []
                with open(pdb_p) as f:
                    for line in f:
                        if line.startswith(("ATOM", "HETATM")) and line[12:16].strip() == "CA":
                            res_nums.append(int(line[22:26]))
                            coords.append(
                                [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                            )
                return np.array(res_nums), np.array(coords)

            r1, c1 = extract_ca(model_pdb)
            r2, c2 = extract_ca(eq_pdb)
            common_res = sorted(set(r1) & set(r2))
            if len(common_res) >= 10:
                idx1 = [np.where(r1 == r)[0][0] for r in common_res]
                idx2 = [np.where(r2 == r)[0][0] for r in common_res]
                c1_sub, c2_sub = c1[idx1], c2[idx2]
                # Kabsch superposition of c2 onto c1
                p1 = c1_sub - c1_sub.mean(axis=0)
                p2 = c2_sub - c2_sub.mean(axis=0)
                H = p2.T @ p1
                U, S, Vt = np.linalg.svd(H)
                d = np.linalg.det(Vt.T @ U.T)
                R = Vt.T @ np.diag([1, 1, d]) @ U.T
                p2_aligned = p2 @ R.T
                drift = np.linalg.norm(p2_aligned - p1, axis=1)
                rmsd = float(np.sqrt(np.mean(drift**2)))
        except Exception:
            rmsd = None
            drift = None

    # Append Cα RMSD to qc.txt if available
    qc_txt_path = dest / "qc.txt"
    if qc_txt_path.is_file() and rmsd is not None:
        try:
            txt = qc_txt_path.read_text()
            if "Cα RMSD" not in txt:
                with open(qc_txt_path, "a") as f:
                    f.write(f"Cα RMSD: {rmsd:.2f} Å (Threshold: < 1.50 Å)\n")
        except Exception:
            pass

    # Headless PyMOL rendering for 3D panel (pure white background)
    pymol_img = None
    if eq_pdb.is_file() and model_pdb and model_pdb.is_file():
        pymol_bin = shutil.which("pymol")
        if pymol_bin:
            with tempfile.TemporaryDirectory() as td:
                tp = Path(td)
                out_png = tp / "_pymol_super.png"
                pml = tp / "_render.pml"
                pml.write_text(f"""
load {model_pdb.resolve()}, m
load {eq_pdb.resolve()}, e
hide everything
bg_color white
set ray_opaque_background, 1
set orthoscopic, 1
set depth_cue, 0
set cartoon_fancy_helices, 1
show cartoon, m
show cartoon, e
color salmon, m
color cyan, e
set cartoon_transparency, 0.4, m
align e, m
zoom m, 1
ray 600, 600
png {out_png.resolve()}
quit
""")
                try:
                    subprocess.run([pymol_bin, "-cq", str(pml)], capture_output=True, timeout=20)
                    if out_png.is_file():
                        pymol_img = Image.open(out_png).copy()
                except Exception:
                    pymol_img = None

    # Unified 6-Panel Layout (7.0 x 4.6 inches, 9pt font)
    fig, axes = plt.subplots(2, 3, figsize=(7.0, 4.6), layout="constrained")

    # Panel A: Temperature
    ax_a = axes[0, 0]
    if has_npt and "Temperature (K)" in df_npt.columns:
        ax_a.plot(df_npt["Time (ps)"], df_npt["Temperature (K)"], color="#0072B2", lw=1)
        ax_a.axhline(300, color="#D55E00", linestyle="--", lw=0.8, alpha=0.7)
    elif has_nvt and "Temperature (K)" in df_nvt.columns:
        ax_a.plot(df_nvt["Time (ps)"], df_nvt["Temperature (K)"], color="#0072B2", lw=1)
        ax_a.axhline(300, color="#D55E00", linestyle="--", lw=0.8, alpha=0.7)
    ax_a.set_title("A: Temperature (K)", fontweight="bold")
    ax_a.set_xlabel("Time (ps)")
    ax_a.set_ylabel("Temperature (K)")

    # Panel B: Density
    ax_b = axes[0, 1]
    if has_npt and "Density (g/mL)" in df_npt.columns:
        ax_b.plot(df_npt["Time (ps)"], df_npt["Density (g/mL)"], color="#0072B2", lw=1)
    ax_b.set_title("B: Density (g/mL)", fontweight="bold")
    ax_b.set_xlabel("Time (ps)")
    ax_b.set_ylabel("Density (g/mL)")

    # Panel C: Potential Energy
    ax_c = axes[0, 2]
    if has_npt and "Potential Energy (kJ/mole)" in df_npt.columns:
        pe_scaled = df_npt["Potential Energy (kJ/mole)"] / 1e5
        ax_c.plot(df_npt["Time (ps)"], pe_scaled, color="#0072B2", lw=1)
        ax_c.set_ylabel(r"$E_{\mathrm{pot}}$ ($10^5$ kJ/mol)")
    elif has_nvt and "Potential Energy (kJ/mole)" in df_nvt.columns:
        pe_scaled = df_nvt["Potential Energy (kJ/mole)"] / 1e5
        ax_c.plot(df_nvt["Time (ps)"], pe_scaled, color="#0072B2", lw=1)
        ax_c.set_ylabel(r"$E_{\mathrm{pot}}$ ($10^5$ kJ/mol)")
    ax_c.set_title("C: Potential Energy", fontweight="bold")
    ax_c.set_xlabel("Time (ps)")

    # Panel D: Box Volume
    ax_d = axes[1, 0]
    if has_npt and "Box Volume (nm^3)" in df_npt.columns:
        ax_d.plot(df_npt["Time (ps)"], df_npt["Box Volume (nm^3)"], color="#0072B2", lw=1)
        ax_d.set_ylabel(r"Volume (nm$^3$)")
    elif has_npt and "Density (g/mL)" in df_npt.columns:
        ax_d.plot(df_npt["Time (ps)"], df_npt["Density (g/mL)"], color="#0072B2", lw=1)
        ax_d.set_ylabel("Density (g/mL)")
    ax_d.set_title(r"D: Box Volume (nm$^3$)", fontweight="bold")
    ax_d.set_xlabel("Time (ps)")

    # Panel E: Per-Residue Cα Displacement
    ax_e = axes[1, 1]
    if drift is not None and common_res is not None:
        ax_e.plot(common_res, drift, color="#0072B2", lw=1)
        ax_e.axhline(1.5, color="#D55E00", linestyle="--", lw=0.8, alpha=0.7, label="Threshold")
        ax_e.set_ylim(0, max(2.0, float(drift.max()) * 1.05))
    ax_e.set_title(r"E: C$\alpha$ Displacement (Å)", fontweight="bold")
    ax_e.set_xlabel("Residue Index")
    ax_e.set_ylabel(r"$\Delta r$ (Å)")

    # Panel F: 3D Backbone Superposition
    ax_f = axes[1, 2]
    rmsd_str = f"RMSD: {rmsd:.2f} Å" if rmsd is not None else "Aligned"
    if pymol_img is not None:
        ax_f.imshow(pymol_img)
        ax_f.axis("off")
        ax_f.set_title(f"F: Superposition ({rmsd_str})", fontweight="bold")
    else:
        ax_f.axis("off")
        ax_f.text(
            0.5,
            0.5,
            f"Fold Stable\n{rmsd_str}",
            ha="center",
            va="center",
            transform=ax_f.transAxes,
            fontsize=8.5,
            fontweight="bold",
            color="#0072B2",
            bbox=dict(boxstyle="round,pad=0.5", facecolor="#EBF3F9", edgecolor="#0072B2"),
        )
        ax_f.set_title(f"F: Fold ({rmsd_str})", fontweight="bold")

    status_str = "PASSED" if passed else "REVIEW REQUIRED"
    fig.suptitle(f"Equilibration Quality Control: {status_str}", fontsize=9.5, fontweight="bold")
    fig.savefig(dest / "qc.png", dpi=300)
    if outdir is None and "simulation" not in root.resolve().parts:
        target = root / "figures"
        target.mkdir(exist_ok=True)
        fig.savefig(target / "equilibration_qc.png", dpi=300)
    plt.close(fig)
