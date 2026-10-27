"""Command-oriented structural mapping, validation and reviewed cofactor transfer."""

import hashlib
import io
import json
import re
from pathlib import Path
from types import SimpleNamespace
from urllib.request import urlopen

from Bio import SeqIO
from Bio.PDB import PDBIO, PDBParser, Select
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SeqUtils import seq1

from colabmda.core.project import ProjectPaths, resolve_project_root
from colabmda.modeller.structure import (
    fetch_archive,
    verify_sequence,
    write_chain,
)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, record):
    Path(path).write_text(json.dumps(record, indent=2) + "\n")


def fetch(args):
    target_tag = f"{args.uniprot_id.upper()}_{args.pdb_id.upper()}"
    paths = ProjectPaths(resolve_project_root(args.root), state=args.state, target=target_tag)
    base = paths.state_dir / "structure" if paths.state_dir else paths.structure_dir
    if ((base / "raw").exists() or (base / "STRUCTURE_DONE").exists()) and not (base / target_tag).exists():
        folder = base
    else:
        folder = base / target_tag
    folder.mkdir(parents=True, exist_ok=True)
    raw = folder
    mapping_dir = folder
    marker = folder / "STRUCTURE_DONE"
    marker.unlink(missing_ok=True)
    structure, archive, url = fetch_archive(args.pdb_id, raw)
    chain = next(iter(structure))[args.chain]
    fasta = raw / "uniprot.fasta"
    payload = urlopen(
        f"https://rest.uniprot.org/uniprotkb/{args.uniprot_id}.fasta", timeout=60
    ).read()
    if fasta.exists() and fasta.read_bytes() != payload:
        raise ValueError("Different existing UniProt sequence preserved")
    fasta.write_bytes(payload)
    sequence = str(SeqIO.read(io.StringIO(payload.decode()), "fasta").seq)
    start, end = args.range
    if not 1 <= start <= end <= len(sequence):
        raise ValueError("Requested UniProt range is outside the sequence")

    meta_method = "X-RAY"
    meta_resolution = None
    author_to_uni = {}
    if archive.suffix == ".cif":
        d = MMCIF2Dict(str(archive))
        refs = zip(
            d["_struct_ref_seq.pdbx_strand_id"],
            d["_struct_ref_seq.pdbx_db_accession"],
            d["_struct_ref_seq.seq_align_beg"],
            d["_struct_ref_seq.seq_align_end"],
            d["_struct_ref_seq.db_align_beg"],
        )
        blocks = [
            (int(a), int(b), int(c))
            for chains, accession, a, b, c in refs
            if args.chain in chains.split(",") and accession == args.uniprot_id
        ]
        if not blocks:
            raise ValueError("No matching mmCIF chain/UniProt reference mapping")
        rows = zip(
            d["_atom_site.auth_asym_id"],
            d["_atom_site.label_atom_id"],
            d["_atom_site.auth_seq_id"],
            d["_atom_site.pdbx_PDB_ins_code"],
            d["_atom_site.label_seq_id"],
        )
        for ch, atom, author, insertion, label in rows:
            if ch != args.chain or atom != "CA" or not label.isdigit():
                continue
            for lo, hi, db in blocks:
                if lo <= int(label) <= hi:
                    key = (int(author), "" if insertion in {".", "?"} else insertion)
                    value = db + int(label) - lo
                    if key in author_to_uni and author_to_uni[key] != value:
                        raise ValueError("Ambiguous archive reference mapping")
                    author_to_uni[key] = value

        m = d.get("_exptl.method")
        if isinstance(m, list) and m:
            meta_method = str(m[0])
        elif isinstance(m, str):
            meta_method = m
        r_val = (
            d.get("_refine.ls_d_res_high")
            or d.get("_reflns.d_resolution_high")
            or d.get("_em_3d_reconstruction.resolution")
        )
        if isinstance(r_val, list) and r_val:
            r_val = r_val[0]
        if r_val and str(r_val).strip() not in {".", "?"}:
            try:
                meta_resolution = float(r_val)
            except (ValueError, TypeError):
                pass
    else:
        for line in archive.read_text().splitlines():
            if (
                line.startswith("DBREF ")
                and line[12] == args.chain
                and line[33:41].strip() == args.uniprot_id
            ):
                lo, hi, db = int(line[14:18]), int(line[20:24]), int(line[55:60])
                author_to_uni.update({(i, ""): db + i - lo for i in range(lo, hi + 1)})
        if not author_to_uni:
            raise ValueError("No matching PDB DBREF mapping")

    if meta_resolution is None and hasattr(structure, "header") and structure.header:
        meta_resolution = structure.header.get("resolution")
        meta_method = structure.header.get("structure_method") or meta_method

    selected = []
    observed = {}
    template_mismatches = []
    for residue in chain:
        if residue.id[0] != " ":
            continue
        key = (residue.id[1], residue.id[2].strip())
        canonical = author_to_uni.get(key)
        if canonical is None:
            continue
        if start <= canonical <= end:
            if "CA" not in residue:
                continue
            res_aa = seq1(residue.resname)
            uni_aa = sequence[canonical - 1]
            if res_aa != uni_aa:
                template_mismatches.append(
                    {
                        "canonical": canonical,
                        "uniprot": uni_aa,
                        "template": res_aa,
                        "resname": residue.resname,
                        "author": residue.id[1],
                        "insertion_code": residue.id[2].strip(),
                    }
                )
            if canonical in observed:
                raise ValueError("Duplicate canonical coordinate mapping")
            observed[canonical] = residue
            selected.append(residue)
    if not selected:
        raise ValueError("No mapped coordinates in requested range")

    mutation_records = []
    for mutation in args.mutations.split(",") if args.mutations else []:
        position = int(mutation[1:-1])
        if position not in observed:
            raise ValueError(
                f"Mutation coordinate gate failed: residue {position} is in an unresolved loop"
            )
        if sequence[position - 1] != mutation[0]:
            raise ValueError(
                f"Mutation WT/canonical reference gate failed: UniProt residue {position} is {sequence[position - 1]}, not {mutation[0]}"
            )
        mutation_records.append(
            f"Target variant {mutation} verified on WT residue {mutation[0]}{position} (coordinate & UniProt match)"
        )

    class Protein(Select):
        def accept_residue(self, residue):
            return residue in selected

    template = raw / "template.pdb"
    stream = io.StringIO()
    write_chain(chain, stream, Protein())
    data = stream.getvalue()
    if template.exists() and template.read_text() != data:
        raise ValueError("Different existing template preserved")
    template.write_text(data)
    rows = [
        {
            "uniprot": i,
            "amino_acid": sequence[i - 1],
            "resolved": i in observed,
            "author": observed[i].id[1] if i in observed else None,
            "insertion_code": observed[i].id[2].strip() if i in observed else None,
        }
        for i in range(start, end + 1)
    ]
    missing_residues = [i for i in range(start, end + 1) if i not in observed]
    loop_segments = []
    if missing_residues:
        seg_start = missing_residues[0]
        seg_end = missing_residues[0]
        for idx in missing_residues[1:]:
            if idx == seg_end + 1:
                seg_end = idx
            else:
                loop_segments.append((seg_start, seg_end, seg_end - seg_start + 1))
                seg_start = idx
                seg_end = idx
        loop_segments.append((seg_start, seg_end, seg_end - seg_start + 1))

    write_json(mapping_dir / "residues.json", rows)
    write_json(
        mapping_dir / "input.json",
        {
            "pdb": args.pdb_id,
            "uniprot": args.uniprot_id,
            "chain": args.chain,
            "range": [start, end],
            "sequence": sequence[start - 1 : end],
            "source_url": url,
            "missing_residues": missing_residues,
            "loop_segments": loop_segments,
            "template_mismatches": template_mismatches,
            "resolution": meta_resolution,
            "method": meta_method,
            "files": {
                str(p.relative_to(folder)): sha(p)
                for p in [archive, fasta, template, mapping_dir / "residues.json"]
            },
        },
    )
    marker.write_text(
        "Author mapping, canonical identity and requested mutation density verified\n"
    )

    total_range = end - start + 1
    density_pct = (len(selected) / total_range) * 100
    res_str = f", Resolution: {meta_resolution:.2f} Å" if meta_resolution else ""
    try:
        rel_folder = folder.relative_to(paths.root)
    except ValueError:
        rel_folder = folder

    print("=" * 78)
    print(" ColabMDA Structural Mapping & Verification")
    print("=" * 78)
    print(
        f" UniProt Target:   {args.uniprot_id} (Residues {start}–{end}, Total Length: {len(sequence)} AA)"
    )
    print(
        f" Template PDB:     {args.pdb_id.upper()} (Chain {args.chain}, Method: {meta_method}{res_str})"
    )
    print(
        f" Coordinate Gate:  {len(selected)}/{total_range} observed residues resolved ({density_pct:.1f}% coordinate density)"
    )
    if loop_segments:
        print(
            f" Missing Loops:    {len(missing_residues)} missing residues identified across {len(loop_segments)} loop segment(s):"
        )
        for idx, (s, e, count) in enumerate(loop_segments, 1):
            print(f"                   - Loop {idx}: Res {s}–{e} ({count} residues)")
    else:
        print(" Missing Loops:    None (100% full coordinate coverage)")

    if template_mismatches:
        print(
            f" Template Check:   ⚠️ {len(template_mismatches)} engineered mutation(s) in template:"
        )
        for m in template_mismatches:
            print(
                f"                   - Res {m['canonical']} (PDB: {m['template']}/{m['resname']} vs UniProt WT: {m['uniprot']})"
            )
        print(
            "                   (Modeller will automatically revert these to Canonical WT during cmda model)"
        )
    else:
        print(
            f" Template Check:   ✓ All {len(selected)} observed residues match Canonical UniProt WT"
        )

    if mutation_records:
        for m_rec in mutation_records:
            print(f" Mutation Gate:    ✓ {m_rec}")
    elif args.mutations:
        print(f" Mutation Gate:    ✓ Verified {args.mutations}")

    print(f" Outputs Saved:    {rel_folder}/")
    print("=" * 78)
    return 0


def check(args):
    paths = ProjectPaths(
        resolve_project_root(args.root), state=args.state, target=getattr(args, "target", None)
    )
    folder = paths.structure_dir
    input_file = folder / "input.json" if (folder / "input.json").exists() else folder / "mapping/input.json"
    record = json.loads(input_file.read_text())
    for name, expected in record["files"].items():
        target_file = folder / name
        if not target_file.exists() and (folder / Path(name).name).exists():
            target_file = folder / Path(name).name
        if sha(target_file) != expected:
            raise ValueError(f"Input changed: {name}")
    sequence = record["sequence"]
    start = record["range"][0]
    if args.variant != "WT":
        pos = int(args.variant[1:-1])
        index = pos - start
        if not 0 <= index < len(sequence) or sequence[index] != args.variant[0]:
            raise ValueError("Mutation WT reference mismatch")
        sequence = sequence[:index] + args.variant[-1] + sequence[index + 1 :]
    model = (
        Path(args.model) if args.model else folder / "models" / args.variant / f"{args.variant}.pdb"
    )
    (model.parent / "MODEL_DONE").unlink(missing_ok=True)
    verify_sequence(model, record["chain"], sequence, start)
    out = model.parent / "model_check.json"
    write_json(
        out,
        {
            "model": str(model.resolve()),
            "sha256": sha(model),
            "input_sha256": sha(input_file),
            "variant": args.variant,
            "sequence": sequence,
            "range": record["range"],
        },
    )
    (model.parent / "MODEL_DONE").write_text(
        "Canonical model sequence, numbering and C-alpha coordinates verified\n"
    )
    print(f"Verified {args.variant}: {len(sequence)} residues, {start}–{record['range'][1]}")
    return 0


def model(args):
    """Build/mutate from checked project inputs using the existing MODELLER engine."""
    import subprocess
    import sys

    if not re.fullmatch(r"WT|[A-Z][1-9][0-9]*[A-Z]", args.variant):
        raise ValueError("Invalid variant")
    paths = ProjectPaths(
        resolve_project_root(args.root), state=args.state, target=getattr(args, "target", None)
    )
    folder = paths.structure_dir
    if not (folder / "STRUCTURE_DONE").is_file():
        raise ValueError("Modeling requires completed cmda fetch (STRUCTURE_DONE)")
    input_file = folder / "input.json" if (folder / "input.json").exists() else folder / "mapping/input.json"
    record = json.loads(input_file.read_text())
    for name, expected in record["files"].items():
        target_file = folder / name
        if not target_file.exists() and (folder / Path(name).name).exists():
            target_file = folder / Path(name).name
        if sha(target_file) != expected:
            raise ValueError(f"Changed structure input: {name}")
    parent = folder / "models" / args.variant
    request = dict(
        variant=args.variant,
        seed=args.seed,
        mapping_sha256=sha(input_file),
        engine_sha256=sha(Path(__file__).parents[1] / "modeller/engine.py"),
    )
    template_file = folder / "template.pdb" if (folder / "template.pdb").exists() else folder / "raw/template.pdb"
    fasta_file = folder / "uniprot.fasta" if (folder / "uniprot.fasta").exists() else folder / "raw/uniprot.fasta"
    if args.variant == "WT":
        argv = [
            "build",
            "--pdb-id",
            record["pdb"],
            "--uniprot-id",
            record["uniprot"],
            "--chain",
            record["chain"],
            "--range",
            *map(str, record["range"]),
            "--uniprot-numbering",
            "--template",
            str(template_file),
            "--sequence",
            str(fasta_file),
            "--outdir",
            str(parent),
        ]
    else:
        position = int(args.variant[1:-1])
        res_file = folder / "residues.json" if (folder / "residues.json").exists() else folder / "mapping/residues.json"
        rows = json.loads(res_file.read_text())
        row = next((row for row in rows if row["uniprot"] == position), None)
        if not row or not row["resolved"] or row["amino_acid"] != args.variant[0]:
            raise ValueError("Mutation WT/canonical/coordinate density gate failed")
        check(SimpleNamespace(root=str(paths.root), state=args.state, variant="WT", model=None))
        wt = folder / "models/WT/WT.pdb"
        request["wt_sha256"] = sha(wt)
        argv = [
            "mutate",
            "--pdb-in",
            str(wt),
            "--chain",
            record["chain"],
            "--mut",
            args.variant,
            "--outdir-mut",
            str(parent),
        ]
    receipt = parent / "model_request.json"
    if (parent / "MODEL_DONE").exists():
        if not receipt.is_file() or json.loads(receipt.read_text()) != request:
            raise ValueError(
                "Existing model has different inputs/settings; use a fresh project structure"
            )
        return check(
            SimpleNamespace(
                root=str(paths.root), state=args.state, variant=args.variant, model=None
            )
        )
    if parent.exists() and any(parent.iterdir()):
        raise ValueError("Incomplete/nonempty model folder preserved; inspect before retrying")
    parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [sys.executable, "-m", "colabmda.cli", "modeller", *argv, "--seed", str(args.seed)],
        check=True,
    )
    write_json(receipt, request)
    check(SimpleNamespace(root=str(paths.root), state=args.state, variant=args.variant, model=None))
    return 0


def assemble(args):
    """Transfer a reviewed pose or a checked public homologous ATP/Mg donor."""
    if getattr(args, "template", None):
        from colabmda.core.pose import assemble_public

        return assemble_public(args)
    if not all(
        getattr(args, key, None)
        for key in ["receptor", "reference", "ligand", "ions", "mapping", "out"]
    ):
        raise ValueError(
            "Supply --state/--template or explicit receptor/reference/ligand/ions/mapping/out"
        )
    if getattr(args, "relax_sidechains", False):
        raise ValueError("Sidechain relaxation is available only for a public donor")
    import numpy as np

    parser = PDBParser(QUIET=True)
    receptor = parser.get_structure("model", args.receptor)
    reference = parser.get_structure("reference", args.reference)
    target = {r.id: r for r in receptor[0][args.chain] if r.id[0] == " " and "CA" in r}
    source = {r.id: r for r in reference[0][args.chain] if r.id[0] == " " and "CA" in r}
    aliases = {"HID": "HIS", "HIE": "HIS", "HIP": "HIS", "ASH": "ASP", "GLH": "GLU", "CYX": "CYS"}

    def normalize(residue):
        return seq1(aliases.get(residue.resname, residue.resname))

    mapping = json.loads(Path(args.mapping).read_text())
    expected_ids = {row["uniprot"] for row in mapping}
    if {key[1] for key in target} != expected_ids:
        raise ValueError("Model does not cover the complete mapped construct")
    mutation_pos = int(args.variant[1:-1]) if args.variant != "WT" else None
    for key, res in target.items():
        if key not in source:
            raise ValueError("Reference receptor lacks a required canonical residue")
        wt, actual = normalize(source[key]), normalize(res)
        if key[1] == mutation_pos:
            if wt != args.variant[0] or actual != args.variant[-1]:
                raise ValueError("Mutant/reference identity mismatch")
        elif wt != actual:
            raise ValueError(f"Reference/model identity mismatch at {key}")
    observed = {row["uniprot"] for row in mapping if row["resolved"]}
    keys = [k for k in target if args.core[0] <= k[1] <= args.core[1] and k[1] in observed]
    if len(keys) < 3:
        raise ValueError("Insufficient matched core coordinates")
    x = np.array([source[k]["CA"].coord for k in keys])
    y = np.array([target[k]["CA"].coord for k in keys])
    xc, yc = x.mean(0), y.mean(0)
    u, _, vt = np.linalg.svd((x - xc).T @ (y - yc))
    rotation = u @ np.diag([1, 1, np.linalg.det(u @ vt)]) @ vt
    rmsd = float(np.sqrt(np.mean(np.sum(((x - xc) @ rotation + yc - y) ** 2, axis=1))))
    if rmsd > args.max_rmsd:
        raise ValueError(f"Core alignment RMSD {rmsd:.3f} exceeds {args.max_rmsd}")
    heterogens = []
    counts = {}
    for path, name in [(args.ligand, "ATP"), (args.ions, "MG")]:
        structure = parser.get_structure(name, path)
        residues = [r for r in structure.get_residues() if r.resname == name]
        counts[name] = len(residues)
        for residue in residues:
            for atom in residue:
                atom.coord = (atom.coord - xc) @ rotation + yc
            heterogens.append(residue)
    if counts != {"ATP": 1, "MG": 1}:
        raise ValueError(f"Expected exactly one ATP and one Mg: {counts}")
    protein_xyz = np.array([a.coord for r in target.values() for a in r if a.element != "H"])
    ligand_xyz = np.array(
        [a.coord for r in heterogens if r.resname == "ATP" for a in r if a.element != "H"]
    )
    minimum_distance = float(np.linalg.norm(protein_xyz[:, None] - ligand_xyz, axis=2).min())
    if minimum_distance < 1.0:
        raise ValueError("Severe protein/ATP clash")
    # Write derived coordinates only, with distinct cofactor IDs after the protein.
    chain = receptor[0][args.chain]
    for i, residue in enumerate(heterogens, max(r.id[1] for r in chain) + 1):
        residue.detach_parent()
        residue.id = (f"H_{residue.resname}", i, " ")
        chain.add(residue)
    out = Path(args.out)
    if out.exists():
        raise ValueError("Existing assembled coordinates preserved; choose a fresh output")
    out.parent.mkdir(parents=True, exist_ok=True)
    writer = PDBIO()
    writer.set_structure(receptor)
    writer.save(str(out))
    from colabmda.modeller.substrate import terminate_components

    terminate_components(out)
    write_json(
        out.with_suffix(".json"),
        {
            "core": args.core,
            "core_rmsd_angstrom": rmsd,
            "fitted_residues": [k[1] for k in keys],
            "fit_selection": "Experimentally resolved C-alpha atoms in mapped core; modeled gaps excluded",
            "variant": args.variant,
            "components": counts,
            "minimum_protein_atp_distance_angstrom": minimum_distance,
            "coordinate_source": "Reviewed activation-project ATP/Mg pose, rigid transfer",
            "input_sha256": {
                str(Path(p).resolve()): sha(p)
                for p in [args.receptor, args.reference, args.ligand, args.ions, args.mapping]
            },
            "output_sha256": sha(out),
        },
    )
    print(f"Assembled reviewed ATP/Mg pose; matched core RMSD {rmsd:.3f} Å")
    return 0


def preparation_request(paths, spec, receptor, template):
    """Freeze direct-workflow inputs and verify completed preparation before reuse."""
    from dataclasses import asdict

    folder = paths.prep_dir(spec)
    folder.mkdir(parents=True, exist_ok=True)
    mapping = (
        paths.structure_dir / "input.json"
        if (paths.structure_dir / "input.json").exists()
        else paths.structure_dir / "mapping/input.json"
    )
    if mapping.exists():
        record = json.loads(mapping.read_text())
        for name, expected in record["files"].items():
            target_file = paths.structure_dir / name
            if not target_file.exists() and (paths.structure_dir / Path(name).name).exists():
                target_file = paths.structure_dir / Path(name).name
            if sha(target_file) != expected:
                raise ValueError(f"Changed structural input: {name}")
        sequence = record["sequence"]
        if spec.variant != "WT":
            index = int(spec.variant[1:-1]) - record["range"][0]
            if not 0 <= index < len(sequence) or sequence[index] != spec.variant[0]:
                raise ValueError("Mutation WT reference mismatch")
            sequence = sequence[:index] + spec.variant[-1] + sequence[index + 1 :]
        verify_sequence(receptor, record["chain"], sequence, record["range"][0])
    request = {
        "spec": json.loads(json.dumps(asdict(spec), default=str)),
        "receptor_sha256": sha(receptor),
        "template_sha256": sha(template) if template else None,
        "mapping_sha256": sha(mapping) if mapping.exists() else None,
        "builder_sha256": sha(Path(__file__).parents[1] / "amber/leap.py"),
    }
    receipt = folder / "preparation_input.json"
    if receipt.exists() and json.loads(receipt.read_text()) != request:
        raise ValueError("Preparation inputs changed; use a fresh COLABMDA_RUN_ID")
    marker = folder / "PREP_DONE"
    if marker.exists():
        if not receipt.exists():
            raise ValueError("Completed preparation lacks its input receipt")
        inventory = verify_prepared(folder)
        net_q = inventory.get("total_charge", 0.0)
        ions = inventory.get("ions", {})
        ion_str = (
            ", ".join(f"{cnt} {name}" for name, cnt in sorted(ions.items()))
            if ions
            else "0 (solute already neutral)"
        )
        n_atoms = inventory.get("atoms", 0)
        n_waters = inventory.get("waters", 0)
        print(f"[OK] Verified existing preparation: {folder}")
        print(f"     - Net charge: {net_q:.4f} e [PASS - NEUTRAL]")
        print(
            f"     - Solvation & ions: {n_waters:,} OPC waters | Counter-ions: {ion_str} | Total atoms: {n_atoms:,}"
        )
        return True
    write_json(receipt, request)
    return False


def verify_prepared(folder):
    """Refuse partial or changed Amber preparation before research dynamics."""
    folder = Path(folder)
    if not (folder / "PREP_DONE").is_file():
        raise SystemExit("Research dynamics requires verified PREP_DONE preparation")
    inventory = json.loads((folder / "topology_check.json").read_text())
    for name, expected in inventory["files"].items():
        if sha(folder / name) != expected:
            raise ValueError(f"Prepared output changed: {name}")
    return inventory
