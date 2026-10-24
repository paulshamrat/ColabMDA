"""Fetch Trp-cage (1L2Y) and select NMR model 1, author chain A, unchanged numbering."""

import hashlib
import io
import json
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

from Bio.PDB import PDBIO, MMCIFParser, PDBParser, Select
from Bio.SeqUtils import seq1

ROOT = Path(__file__).resolve().parents[2]
folder = ROOT / "data/01_structure/trpcage"
folder.mkdir(parents=True, exist_ok=True)
source = None
for suffix in ["cif", "pdb"]:
    url = f"https://files.rcsb.org/download/1L2Y.{suffix}"
    try:
        payload = urlopen(url, timeout=60).read()
    except (HTTPError, URLError, TimeoutError):
        continue
    parser = MMCIFParser(QUIET=True) if suffix == "cif" else PDBParser(QUIET=True)
    structure = parser.get_structure("1L2Y", io.StringIO(payload.decode()))
    source = url
    break
if source is None:
    raise RuntimeError("Neither mmCIF nor PDB retrieval succeeded")
model = next(iter(structure))
chain = model["A"]
residues = [r for r in chain if r.id[0] == " "]
assert len(residues) == 20
assert "".join(seq1(r.resname) for r in residues) == "NLYIQWLKDGGPSSGRPPPS"
assert [(r.id[1], r.id[2]) for r in residues] == [(i, " ") for i in range(1, 21)]
assert all("CA" in r for r in residues)


class Protein(Select):
    def accept_residue(self, residue):
        return residue.id[0] == " "


output = io.StringIO()
writer = PDBIO()
writer.set_structure(chain)
writer.save(output, Protein())
coordinates = output.getvalue().encode()
target = folder / "WT.pdb"
if target.exists() and target.read_bytes() != coordinates:
    raise RuntimeError("Different existing input preserved")
raw = folder / ("1L2Y." + suffix)
if raw.exists() and raw.read_bytes() != payload:
    raise RuntimeError("Different existing archive preserved")
raw.write_bytes(payload)
target.write_bytes(coordinates)
record = {
    "pdb": "1L2Y",
    "source_url": source,
    "archive_sha256": hashlib.sha256(payload).hexdigest(),
    "input_sha256": hashlib.sha256(coordinates).hexdigest(),
    "input": "data/01_structure/trpcage/WT.pdb",
    "model": 1,
    "author_chain": "A",
    "residues": 20,
    "sequence": "NLYIQWLKDGGPSSGRPPPS",
    "selection": "First NMR model, author chain A, standard protein residues; native numbering",
}
(folder / "input.json").write_text(json.dumps(record, indent=2) + "\n")
(folder / "STRUCTURE_DONE").write_text(
    "Deposited sequence, native numbering and coordinate presence verified\n"
)
print(json.dumps(record, indent=2))
