#!/usr/bin/env python3
"""
truncate_model.py

Truncate a PDB to a user‐specified residue range, preserving CRYST1 and other header lines.

Usage:
    python3 truncate_model.py \
      --input   model_with_cryst.pdb \
      --chain   A \
      --start   9 \
      --end     170 \
      --output  truncated.pdb
"""

import argparse
from Bio.PDB import PDBParser, PDBIO, Select

class RangeSelect(Select):
    def __init__(self, chain_id, start, end):
        self.chain = chain_id
        self.start = start
        self.end = end

    def accept_chain(self, chain):
        return chain.id == self.chain

    def accept_residue(self, residue):
        return self.start <= residue.id[1] <= self.end

    def accept_atom(self, atom):
        return True

def main():
    p = argparse.ArgumentParser(description="Truncate a PDB by residue range")
    p.add_argument("--input",  required=True, help="Input PDB (with CRYST1)")
    p.add_argument("--chain",  required=True, help="Chain ID to keep")
    p.add_argument("--start",  type=int, required=True, help="First residue to keep")
    p.add_argument("--end",    type=int, required=True, help="Last  residue to keep")
    p.add_argument("--output", required=True, help="Output truncated PDB")
    args = p.parse_args()

    # Preserve header lines (CRYST1, HEADER, REMARK, etc.)
    header = []
    with open(args.input) as f:
        for line in f:
            if line.startswith(("CRYST1", "HEADER", "TITLE", "REMARK")):
                header.append(line)
            else:
                break

    # Parse structure
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("model", args.input)

    # Write out truncated model
    io = PDBIO()
    io.set_structure(struct)
    with open(args.output, "w") as out:
        for h in header:
            out.write(h)
        io.save(out, select=RangeSelect(args.chain, args.start, args.end))

    print(f"Truncated model saved to {args.output}")

if __name__ == "__main__":
    main()
