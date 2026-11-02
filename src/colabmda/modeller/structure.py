"""Archive ingestion and strict coordinate checks, independent of MODELLER."""

import io
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

from Bio.PDB import PDBIO, MMCIFParser, PDBParser
from Bio.SeqUtils import seq1


def fetch_archive(pdb_id, folder):
    """Keep the original archive; prefer mmCIF and fall back to legacy PDB."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    for suffix in ("cif", "pdb"):
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.{suffix}"
        try:
            payload = urlopen(url, timeout=60).read()
        except (HTTPError, URLError, TimeoutError):
            continue
        parser = MMCIFParser(QUIET=True) if suffix == "cif" else PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_id, io.StringIO(payload.decode()))
        path = folder / f"{pdb_id.lower()}.{suffix}"
        if path.exists() and path.read_bytes() != payload:
            raise ValueError(f"Changed archive preserved: {path}; use a fresh input directory")
        path.write_bytes(payload)
        return structure, path, url
    raise RuntimeError(f"Neither mmCIF nor PDB retrieval succeeded for {pdb_id}")


def write_chain(chain, path, select=None):
    """Write a derived single-chain MODELLER input; never renumber the archive."""
    writer = PDBIO()
    writer.set_structure(chain)
    target = path if hasattr(path, "write") else str(path)
    if select is None:
        writer.save(target)
    else:
        writer.save(target, select)


def protein_residues(path, chain_id="A"):
    model = next(iter(PDBParser(QUIET=True).get_structure("input", str(path))))
    if chain_id not in model:
        raise ValueError(f"Requested author chain {chain_id} absent from {path}")
    residues = [r for r in model[chain_id] if r.id[0] == " "]
    if not residues or any("CA" not in r for r in residues):
        raise ValueError("Missing protein residues or C-alpha coordinates")
    return residues


def verify_sequence(path, chain_id, sequence, start):
    residues = protein_residues(path, chain_id)
    expected_ids = [(i, " ") for i in range(start, start + len(sequence))]
    if [(r.id[1], r.id[2]) for r in residues] != expected_ids:
        raise ValueError("Model residue count/numbering does not match the requested reference")
    if "".join(seq1(r.resname) for r in residues) != sequence:
        raise ValueError("Model amino acids do not match the requested reference")
    return residues
