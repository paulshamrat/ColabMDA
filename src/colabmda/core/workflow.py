"""Command-oriented structural mapping, validation and reviewed cofactor transfer."""

import hashlib
import io
import json
from pathlib import Path
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
    paths = ProjectPaths(resolve_project_root(args.root), state=args.state)
    folder = paths.structure_dir
    raw = folder / "raw"
    mapping_dir = folder / "mapping"
    raw.mkdir(parents=True, exist_ok=True)
    mapping_dir.mkdir(parents=True, exist_ok=True)
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
    selected = []
    observed = {}
    for residue in chain:
        if residue.id[0] != " ":
            continue
        key = (residue.id[1], residue.id[2].strip())
        canonical = author_to_uni.get(key)
        if canonical is None:
            continue
        if start <= canonical <= end:
            if "CA" not in residue or seq1(residue.resname) != sequence[canonical - 1]:
                raise ValueError(f"Coordinate/UniProt identity or density mismatch at {key}")
            if canonical in observed:
                raise ValueError("Duplicate canonical coordinate mapping")
            observed[canonical] = residue
            selected.append(residue)
    if not selected:
        raise ValueError("No mapped coordinates in requested range")
    for mutation in args.mutations.split(",") if args.mutations else []:
        position = int(mutation[1:-1])
        if position not in observed or sequence[position - 1] != mutation[0]:
            raise ValueError(f"Mutation WT/canonical/coordinate gate failed: {mutation}")

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
            "missing_residues": [i for i in range(start, end + 1) if i not in observed],
            "files": {
                str(p.relative_to(folder)): sha(p)
                for p in [archive, fasta, template, mapping_dir / "residues.json"]
            },
        },
    )
    marker.write_text(
        "Author mapping, canonical identity and requested mutation density verified\n"
    )
    print(
        f"Mapped {len(selected)} observed residues; {end - start + 1 - len(selected)} missing residues recorded"
    )
    return 0


def check(args):
    paths = ProjectPaths(resolve_project_root(args.root), state=args.state)
    folder = paths.structure_dir
    record = json.loads((folder / "mapping/input.json").read_text())
    for name, expected in record["files"].items():
        if sha(folder / name) != expected:
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
            "input_sha256": sha(folder / "mapping/input.json"),
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


def assemble(args):
    """Transfer a reviewed ATP/Mg pose by matching canonical CDKL5 core residues."""
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
    mapping = paths.structure_dir / "mapping/input.json"
    if mapping.exists():
        record = json.loads(mapping.read_text())
        for name, expected in record["files"].items():
            if sha(paths.structure_dir / name) != expected:
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
        verify_prepared(folder)
        print(f"Verified existing preparation: {folder}")
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
