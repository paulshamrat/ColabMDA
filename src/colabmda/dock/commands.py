#!/usr/bin/env python3
"""Small, explicit docking workflow helpers.

The docking layer intentionally wraps standard command-line tools instead of
trying to hide the chemistry.  Inputs and outputs stay visible on disk:

1. build or provide a ligand library (SDF/SMILES)
2. prepare receptor and ligands as PDBQT
3. dock with AutoDock Vina
4. export ranked poses back to SDF for MD setup
"""

from __future__ import annotations

import csv
import json
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

DEMO_LIBRARY = [
    ("aspirin", "CC(=O)Oc1ccccc1C(=O)O", "small real-molecule smoke library"),
    ("caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "small real-molecule smoke library"),
    ("acetaminophen", "CC(=O)Nc1ccc(O)cc1", "small real-molecule smoke library"),
    ("ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O", "small real-molecule smoke library"),
]


@dataclass
class DockingHit:
    ligand_id: str
    affinity_kcal_mol: float
    output_pdbqt: Path
    log_file: Path


def _require_executable(name: str) -> str:
    exe = shutil.which(name)
    if not exe:
        raise SystemExit(
            f"ERROR: required command not found: {name}\n"
            "Install the docking tools first, for example with the Colab bootstrap:\n"
            "  WITH_DOCKING=1 bash ./bootstrap_colab_openmm_gpu.sh\n"
        )
    return exe


def _safe_id(text: str) -> str:
    clean = re.sub(r"[^A-Za-z0-9_.-]+", "_", text.strip())
    return clean.strip("._-") or "ligand"


def _read_smiles_table(path: Path) -> list[tuple[str, str, str]]:
    rows: list[tuple[str, str, str]] = []
    with path.open(newline="") as handle:
        sniff = handle.read(2048)
        handle.seek(0)
        has_header = csv.Sniffer().has_header(sniff) if sniff.strip() else False
        if has_header:
            reader = csv.DictReader(handle)
            for idx, row in enumerate(reader, start=1):
                smiles = row.get("smiles") or row.get("SMILES") or row.get("canonical_smiles")
                name = row.get("id") or row.get("name") or row.get("compound_id") or f"ligand_{idx}"
                source = row.get("source") or ""
                if smiles:
                    rows.append((_safe_id(name), smiles.strip(), source))
        else:
            reader = csv.reader(handle)
            for idx, row in enumerate(reader, start=1):
                if not row:
                    continue
                smiles = row[0].strip()
                name = row[1].strip() if len(row) > 1 and row[1].strip() else f"ligand_{idx}"
                source = row[2].strip() if len(row) > 2 else ""
                rows.append((_safe_id(name), smiles, source))
    return rows


def dock_library(
    outdir: str,
    preset: str = "demo",
    smiles_csv: str | None = None,
    max_mols: int | None = None,
    seed: int = 2026,
) -> None:
    """Create a small 3D SDF library from demo molecules or a SMILES table."""

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except Exception as exc:
        raise SystemExit(
            "ERROR: RDKit is required to build ligand libraries.\n"
            "Install docking extras with:\n"
            "  WITH_DOCKING=1 bash ./bootstrap_colab_openmm_gpu.sh"
        ) from exc

    out = Path(outdir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    library_sdf = out / "library.sdf"
    metadata_json = out / "library_metadata.json"
    manifest_csv = out / "library_manifest.csv"

    if smiles_csv:
        source = str(Path(smiles_csv).resolve())
        entries = _read_smiles_table(Path(smiles_csv))
    elif preset == "demo":
        source = "built-in demo library; not a screening database"
        entries = DEMO_LIBRARY.copy()
    else:
        raise SystemExit("ERROR: currently supported --preset values: demo")

    if max_mols is not None:
        entries = entries[:max_mols]

    written = []
    writer = Chem.SDWriter(str(library_sdf))
    for idx, (ligand_id, smiles, source_note) in enumerate(entries, start=1):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"[WARN] Skipping invalid SMILES for {ligand_id}: {smiles}")
            continue
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = int(seed + idx)
        if AllChem.EmbedMolecule(mol, params) != 0:
            print(f"[WARN] Could not generate 3D conformer for {ligand_id}")
            continue
        AllChem.UFFOptimizeMolecule(mol, maxIters=300)
        mol.SetProp("_Name", ligand_id)
        mol.SetProp("SMILES", smiles)
        if source_note:
            mol.SetProp("SOURCE", source_note)
        writer.write(mol)
        written.append((ligand_id, smiles, source_note))
    writer.close()

    with manifest_csv.open("w", newline="") as handle:
        writer_csv = csv.writer(handle)
        writer_csv.writerow(["ligand_id", "smiles", "source"])
        writer_csv.writerows(written)

    metadata = {
        "source": source,
        "preset": preset if not smiles_csv else None,
        "n_requested": len(entries),
        "n_written": len(written),
        "sdf": str(library_sdf),
        "manifest": str(manifest_csv),
        "note": (
            "The demo preset is only for smoke tests. For real virtual screening, "
            "download molecules from a named database/query, save SMILES or SDF, "
            "and keep that file with the campaign."
        ),
    }
    metadata_json.write_text(json.dumps(metadata, indent=2) + "\n")

    print("[OK] Ligand library written")
    print(f"  SDF      : {library_sdf}")
    print(f"  Manifest : {manifest_csv}")
    print(f"  Molecules: {len(written)}")


def dock_prep_receptor(
    receptor: str,
    outdir: str,
    name: str,
    center: list[float],
    size: list[float],
) -> None:
    """Prepare receptor PDBQT and box files with Meeko."""

    mk_prepare_receptor = _require_executable("mk_prepare_receptor.py")
    receptor_path = Path(receptor).resolve()
    if not receptor_path.is_file():
        raise SystemExit(f"ERROR: receptor file not found: {receptor_path}")

    out = Path(outdir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    prefix = out / name
    cmd = [
        mk_prepare_receptor,
        "-i",
        str(receptor_path),
        "-o",
        str(prefix),
        "-p",
        "-v",
        "--box_center",
        *(str(v) for v in center),
        "--box_size",
        *(str(v) for v in size),
    ]
    print("[RUN]", " ".join(cmd))
    subprocess.check_call(cmd)

    config = out / "vina_box.txt"
    config.write_text(
        "\n".join(
            [
                f"center_x = {center[0]}",
                f"center_y = {center[1]}",
                f"center_z = {center[2]}",
                f"size_x = {size[0]}",
                f"size_y = {size[1]}",
                f"size_z = {size[2]}",
                "",
            ]
        )
    )
    metadata = {
        "receptor": str(receptor_path),
        "receptor_pdbqt": str(prefix.with_suffix(".pdbqt")),
        "box_center_angstrom": center,
        "box_size_angstrom": size,
        "vina_config": str(config),
    }
    (out / "receptor_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print("[OK] Receptor prepared")
    print(f"  PDBQT : {prefix.with_suffix('.pdbqt')}")
    print(f"  Box   : {config}")


def _pdb_atom_coords(path: Path, resname: str | None = None, chain: str | None = None):
    res_filter = resname.upper() if resname else None
    chain_filter = chain.strip() if chain else None
    with path.open() as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            residue_names = {
                line[17:20].strip().upper(),
                line[16:20].strip().upper(),
                line[17:21].strip().upper(),
            }
            if res_filter and res_filter not in residue_names:
                parts = line.split()
                if res_filter not in {part.upper() for part in parts[2:5]}:
                    continue
            if chain_filter and line[21].strip() != chain_filter:
                continue
            try:
                yield (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue


def dock_box_from_pdb(
    pdb: str,
    out: str,
    resname: str | None = None,
    chain: str | None = None,
    padding: float = 8.0,
    min_size: float = 18.0,
) -> None:
    """Write a Vina box config from a PDB atom selection."""

    pdb_path = Path(pdb).resolve()
    if not pdb_path.is_file():
        raise SystemExit(f"ERROR: PDB file not found: {pdb_path}")

    coords = list(_pdb_atom_coords(pdb_path, resname=resname, chain=chain))
    if not coords:
        selection = []
        if resname:
            selection.append(f"resname={resname}")
        if chain:
            selection.append(f"chain={chain}")
        selection_text = ", ".join(selection) if selection else "all atoms"
        raise SystemExit(f"ERROR: no atoms matched selection ({selection_text}) in {pdb_path}")

    xs, ys, zs = zip(*coords)
    center = [sum(axis) / len(axis) for axis in (xs, ys, zs)]
    spans = [max(axis) - min(axis) for axis in (xs, ys, zs)]
    size = [max(float(min_size), span + 2.0 * float(padding)) for span in spans]

    out_path = Path(out).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(
        "\n".join(
            [
                f"center_x = {center[0]:.3f}",
                f"center_y = {center[1]:.3f}",
                f"center_z = {center[2]:.3f}",
                f"size_x = {size[0]:.3f}",
                f"size_y = {size[1]:.3f}",
                f"size_z = {size[2]:.3f}",
                "",
            ]
        )
    )

    metadata = {
        "pdb": str(pdb_path),
        "resname": resname,
        "chain": chain,
        "n_atoms": len(coords),
        "padding_angstrom": padding,
        "min_size_angstrom": min_size,
        "center_angstrom": center,
        "size_angstrom": size,
        "config": str(out_path),
    }
    out_path.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")

    print("[OK] Docking box written")
    print(f"  Config : {out_path}")
    print(f"  Center : {center[0]:.3f} {center[1]:.3f} {center[2]:.3f}")
    print(f"  Size   : {size[0]:.3f} {size[1]:.3f} {size[2]:.3f}")


def _iter_sdf_records(sdf: Path):
    text = sdf.read_text(errors="ignore")
    for idx, block in enumerate(text.split("$$$$"), start=1):
        block = block.strip()
        if not block:
            continue
        lines = block.splitlines()
        name = _safe_id(lines[0] if lines else f"ligand_{idx}")
        yield name, block + "\n$$$$\n"


def dock_prep_ligands(library_sdf: str, outdir: str) -> None:
    """Prepare one ligand PDBQT per SDF record."""

    mk_prepare_ligand = _require_executable("mk_prepare_ligand.py")
    library = Path(library_sdf).resolve()
    if not library.is_file():
        raise SystemExit(f"ERROR: SDF library not found: {library}")
    out = Path(outdir).resolve()
    sdf_dir = out / "sdf"
    pdbqt_dir = out / "pdbqt"
    sdf_dir.mkdir(parents=True, exist_ok=True)
    pdbqt_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for idx, (ligand_id, sdf_block) in enumerate(_iter_sdf_records(library), start=1):
        unique_id = _safe_id(ligand_id)
        if any(row[0] == unique_id for row in rows):
            unique_id = f"{unique_id}_{idx}"
        sdf_path = sdf_dir / f"{unique_id}.sdf"
        pdbqt_path = pdbqt_dir / f"{unique_id}.pdbqt"
        sdf_path.write_text(sdf_block)
        cmd = [mk_prepare_ligand, "-i", str(sdf_path), "-o", str(pdbqt_path)]
        print("[RUN]", " ".join(cmd))
        subprocess.check_call(cmd)
        rows.append((unique_id, str(sdf_path), str(pdbqt_path)))

    manifest = out / "prepared_ligands.csv"
    with manifest.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["ligand_id", "sdf", "pdbqt"])
        writer.writerows(rows)

    print("[OK] Ligands prepared")
    print(f"  PDBQT dir: {pdbqt_dir}")
    print(f"  Manifest : {manifest}")
    print(f"  Count    : {len(rows)}")


def _parse_best_affinity(log_text: str) -> float | None:
    for line in log_text.splitlines():
        match = re.match(r"\s*1\s+(-?\d+(?:\.\d+)?)\s+", line)
        if match:
            return float(match.group(1))
    return None


def dock_run(
    receptor_pdbqt: str,
    ligands_dir: str,
    config: str,
    outdir: str,
    exhaustiveness: int = 8,
    num_modes: int = 9,
    cpu: int = 0,
    scoring: str = "vina",
) -> None:
    """Run AutoDock Vina over prepared PDBQT ligands."""

    vina = _require_executable("vina")
    receptor = Path(receptor_pdbqt).resolve()
    ligand_dir = Path(ligands_dir).resolve()
    config_path = Path(config).resolve()
    if not receptor.is_file():
        raise SystemExit(f"ERROR: receptor PDBQT not found: {receptor}")
    if not ligand_dir.is_dir():
        raise SystemExit(f"ERROR: ligand PDBQT directory not found: {ligand_dir}")
    if not config_path.is_file():
        raise SystemExit(f"ERROR: Vina config file not found: {config_path}")

    out = Path(outdir).resolve()
    poses_dir = out / "poses_pdbqt"
    logs_dir = out / "logs"
    poses_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    hits: list[DockingHit] = []
    ligands = sorted(ligand_dir.glob("*.pdbqt"))
    if not ligands:
        raise SystemExit(f"ERROR: no .pdbqt ligands found in {ligand_dir}")

    for ligand in ligands:
        ligand_id = ligand.stem
        pose = poses_dir / f"{ligand_id}_out.pdbqt"
        log_file = logs_dir / f"{ligand_id}.log"
        cmd = [
            vina,
            "--receptor",
            str(receptor),
            "--ligand",
            str(ligand),
            "--config",
            str(config_path),
            "--exhaustiveness",
            str(exhaustiveness),
            "--num_modes",
            str(num_modes),
            "--cpu",
            str(cpu),
            "--scoring",
            scoring,
            "--out",
            str(pose),
        ]
        print("[RUN]", " ".join(cmd))
        completed = subprocess.run(cmd, text=True, capture_output=True, check=False)
        log_text = completed.stdout + completed.stderr
        log_file.write_text(log_text)
        if completed.returncode != 0:
            print(f"[WARN] Docking failed for {ligand_id}; see {log_file}")
            continue
        affinity = _parse_best_affinity(log_text)
        if affinity is None:
            print(f"[WARN] Could not parse affinity for {ligand_id}; see {log_file}")
            continue
        hits.append(DockingHit(ligand_id, affinity, pose, log_file))

    hits = sorted(hits, key=lambda hit: hit.affinity_kcal_mol)
    ranked_csv = out / "ranked_hits.csv"
    with ranked_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["rank", "ligand_id", "affinity_kcal_mol", "pose_pdbqt", "log"])
        for rank, hit in enumerate(hits, start=1):
            writer.writerow(
                [
                    rank,
                    hit.ligand_id,
                    hit.affinity_kcal_mol,
                    str(hit.output_pdbqt),
                    str(hit.log_file),
                ]
            )

    metadata = {
        "engine": "AutoDock Vina",
        "receptor_pdbqt": str(receptor),
        "ligands_dir": str(ligand_dir),
        "config": str(config_path),
        "exhaustiveness": exhaustiveness,
        "num_modes": num_modes,
        "cpu": cpu,
        "scoring": scoring,
        "ranked_hits": str(ranked_csv),
    }
    (out / "docking_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")

    print("[OK] Docking complete")
    print(f"  Ranked hits: {ranked_csv}")
    if hits:
        print(f"  Best hit   : {hits[0].ligand_id} ({hits[0].affinity_kcal_mol:.2f} kcal/mol)")


def _read_ranked_hits(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _keep_first_sdf_record(path: Path) -> None:
    """Keep only the first SDF molecule record in a pose file."""

    text = path.read_text(errors="ignore")
    records = [record.strip() for record in text.split("$$$$") if record.strip()]
    if len(records) <= 1:
        return
    path.write_text(records[0] + "\n$$$$\n")


def dock_export_top(results: str, outdir: str, top_n: int = 3) -> None:
    """Export top Vina PDBQT poses to SDF so they can be used as MD ligands."""

    mk_export = _require_executable("mk_export.py")
    results_path = Path(results).resolve()
    ranked = results_path / "ranked_hits.csv" if results_path.is_dir() else results_path
    if not ranked.is_file():
        raise SystemExit(f"ERROR: ranked hits CSV not found: {ranked}")

    out = Path(outdir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    exported_rows = []
    for row in _read_ranked_hits(ranked)[:top_n]:
        ligand_id = row["ligand_id"]
        pose_pdbqt = Path(row["pose_pdbqt"]).resolve()
        pose_sdf = out / f"{int(row['rank']):02d}_{ligand_id}_pose.sdf"
        cmd = [mk_export, str(pose_pdbqt), "-s", str(pose_sdf)]
        print("[RUN]", " ".join(cmd))
        subprocess.check_call(cmd)
        _keep_first_sdf_record(pose_sdf)
        exported_rows.append(
            [row["rank"], ligand_id, row["affinity_kcal_mol"], str(pose_sdf), str(pose_pdbqt)]
        )

    selected_csv = out / "selected_poses.csv"
    with selected_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["rank", "ligand_id", "affinity_kcal_mol", "pose_sdf", "pose_pdbqt"])
        writer.writerows(exported_rows)

    print("[OK] Top docked poses exported for MD")
    print(f"  Selected poses: {selected_csv}")
    if exported_rows:
        print("  Example MD handoff:")
        print(
            "  colabmda openmm run --name receptor_"
            f"{exported_rows[0][1]} --replica r1 --ligand {exported_rows[0][3]} --total-ns 10.0"
        )


def dock_report(results: str, outdir: str | None = None, top_n: int = 20) -> None:
    """Create lightweight plots and a Markdown summary for docking results."""

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        raise SystemExit(
            "ERROR: matplotlib is required for docking reports.\n"
            "Install analysis extras or use the bootstrap environment."
        ) from exc

    results_path = Path(results).resolve()
    ranked = results_path / "ranked_hits.csv" if results_path.is_dir() else results_path
    if not ranked.is_file():
        raise SystemExit(f"ERROR: ranked hits CSV not found: {ranked}")

    rows = _read_ranked_hits(ranked)
    if not rows:
        raise SystemExit(f"ERROR: ranked hits CSV has no hits: {ranked}")

    report_dir = Path(outdir).resolve() if outdir else ranked.parent / "report"
    report_dir.mkdir(parents=True, exist_ok=True)

    subset = rows[:top_n]
    names = [row["ligand_id"] for row in subset]
    scores = [float(row["affinity_kcal_mol"]) for row in subset]

    fig_height = max(3.5, 0.35 * len(subset) + 1.5)
    fig, ax = plt.subplots(figsize=(8, fig_height))
    y = list(range(len(subset)))
    ax.barh(y, scores, color="#2b83ba")
    ax.set_yticks(y)
    ax.set_yticklabels(names)
    ax.invert_yaxis()
    ax.set_xlabel("Vina score (kcal/mol)")
    ax.set_title(f"Top {len(subset)} docked ligands")
    ax.axvline(0, color="black", linewidth=0.8)
    fig.tight_layout()
    score_plot = report_dir / "docking_scores.png"
    fig.savefig(score_plot, dpi=200)
    plt.close(fig)

    score_csv = report_dir / "top_hits.csv"
    with score_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(subset)

    best = subset[0]
    markdown = report_dir / "docking_report.md"
    markdown.write_text(
        "\n".join(
            [
                "# Docking Report",
                "",
                f"- Ranked hits: `{ranked}`",
                f"- Hits plotted: {len(subset)}",
                f"- Best hit: `{best['ligand_id']}` ({float(best['affinity_kcal_mol']):.2f} kcal/mol)",
                "",
                "![Docking scores](docking_scores.png)",
                "",
                "## Top hits",
                "",
                "| Rank | Ligand | Vina score (kcal/mol) |",
                "|---:|---|---:|",
                *[
                    f"| {row['rank']} | {row['ligand_id']} | {float(row['affinity_kcal_mol']):.2f} |"
                    for row in subset
                ],
                "",
                "Scores are useful for prioritizing poses within the same campaign, not as experimental binding free energies.",
                "",
            ]
        )
    )

    print("[OK] Docking report written")
    print(f"  Plot   : {score_plot}")
    print(f"  Summary: {markdown}")
    print(f"  CSV    : {score_csv}")


def dock_merge_results(campaign_dir: str, outdir: str | None = None) -> None:
    """Merge ranked_hits.csv files from many docking shards."""

    campaign = Path(campaign_dir).resolve()
    if not campaign.is_dir():
        raise SystemExit(f"ERROR: campaign directory not found: {campaign}")

    ranked_files = sorted(campaign.glob("shards/*/results/ranked_hits.csv"))
    if not ranked_files:
        ranked_files = sorted(campaign.glob("*/results/ranked_hits.csv"))
    if not ranked_files:
        raise SystemExit(
            f"ERROR: no shard ranked_hits.csv files found below {campaign}\n"
            "Expected paths like: data/docking/<campaign>/shards/000001/results/ranked_hits.csv"
        )

    merged = []
    for ranked in ranked_files:
        shard = ranked.parent.parent.name
        for row in _read_ranked_hits(ranked):
            merged.append(
                {
                    "shard": shard,
                    "ligand_id": row["ligand_id"],
                    "affinity_kcal_mol": float(row["affinity_kcal_mol"]),
                    "pose_pdbqt": row["pose_pdbqt"],
                    "log": row["log"],
                }
            )

    merged = sorted(merged, key=lambda row: row["affinity_kcal_mol"])
    output_dir = Path(outdir).resolve() if outdir else campaign / "merged"
    output_dir.mkdir(parents=True, exist_ok=True)
    merged_csv = output_dir / "ranked_hits.csv"
    with merged_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["rank", "shard", "ligand_id", "affinity_kcal_mol", "pose_pdbqt", "log"],
        )
        writer.writeheader()
        for rank, row in enumerate(merged, start=1):
            writer.writerow({"rank": rank, **row})

    metadata = {
        "campaign_dir": str(campaign),
        "n_shards": len(ranked_files),
        "n_hits": len(merged),
        "ranked_hits": str(merged_csv),
        "shard_files": [str(path) for path in ranked_files],
    }
    (output_dir / "merge_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")

    print("[OK] Shard results merged")
    print(f"  Shards : {len(ranked_files)}")
    print(f"  Hits   : {len(merged)}")
    print(f"  Output : {merged_csv}")
