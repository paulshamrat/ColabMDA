#!/usr/bin/env python3
"""
modeller4.py

Modeller homology modeling pipeline CLI with logging.
Usage:
    python3 modeller4.py <PDB_ID> <UNIPROT_ID>

This script will:
 1) Create folder named after the PDB ID (e.g. “4ldj/”)
 2) Download the PDB and UniProt FASTA into that folder
 3) Clean the PDB (remove HETATM, keep CRYST1)
 4) Extract the actual template sequence from the ATOM records
 5) Write an alignment.ali
 6) Run Modeller’s automodel to build one model (capturing its output)
 7) Insert the original CRYST1 record back into the final PDB
All console output (including Modeller’s) is saved to pipeline.log in that folder.
"""
import os
import sys
import subprocess
from datetime import datetime
from contextlib import redirect_stdout, redirect_stderr
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.Data.IUPACData import protein_letters_3to1 as P3
from Bio import SeqIO
from modeller import Environ, log
from modeller.automodel import automodel, assess

def timestamp():
    return datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC")

def run(cmd, logfh):
    line = f"[{timestamp()}] → {cmd}\n"
    print(line, end='')
    logfh.write(line)
    subprocess.run(cmd, shell=True, check=True, stdout=logfh, stderr=logfh)

if len(sys.argv) != 3:
    print("Usage: python3 modeller4.py <PDB_ID> <UNIPROT_ID>")
    sys.exit(1)

pdb_id = sys.argv[1].lower()
uni_id = sys.argv[2].upper()

# 1) Create output dir and open log
outdir = pdb_id
os.makedirs(outdir, exist_ok=True)
log_path = os.path.join(outdir, "pipeline.log")
logfh = open(log_path, "w")
logfh.write(f"[{timestamp()}] Pipeline started\n")

# Switch into it
os.chdir(outdir)

# 2) Download inputs
run(f"wget -q https://files.rcsb.org/download/{pdb_id}.pdb -O {pdb_id}_orig.pdb", logfh)
run(f"wget -q https://www.uniprot.org/uniprot/{uni_id}.fasta -O {uni_id}.fasta", logfh)

# 3) Clean PDB: remove HETATM, keep CRYST1
logfh.write(f"[{timestamp()}] Cleaning PDB\n")
cryst1 = ""
with open(f"{pdb_id}_orig.pdb") as fin:
    for L in fin:
        if L.startswith("CRYST1"):
            cryst1 = L
            break

parser = PDBParser(QUIET=True)
struct = parser.get_structure(pdb_id, f"{pdb_id}_orig.pdb")
class KeepProtein(Select):
    def accept_residue(self, r): return r.id[0] == ' '
io = PDBIO(); io.set_structure(struct)
io.save(f"{pdb_id}_clean.pdb", select=KeepProtein())

with open(f"{pdb_id}.pdb","w") as fout, open(f"{pdb_id}_clean.pdb") as fin:
    fout.write(cryst1)
    for L in fin:
        if not L.startswith("CRYST1"):
            fout.write(L)
open("cryst1.txt","w").write(cryst1)

# 4) Extract template sequence from ATOM
logfh.write(f"[{timestamp()}] Extracting template sequence\n")
struct2 = parser.get_structure(pdb_id, f"{pdb_id}.pdb")
chain = next(struct2[0].get_chains())
obs = {r.id[1]:P3.get(r.get_resname().capitalize(),"X") for r in chain if r.id[0]==' '}
start, end = min(obs), max(obs)
template_seq = "".join(obs.get(i,"-") for i in range(start, end+1))

# 5) Read full UniProt sequence
target_seq = str(SeqIO.read(f"{uni_id}.fasta","fasta").seq)
length = len(target_seq)
logfh.write(f"[{timestamp()}] UniProt length: {length}\n")
logfh.write(f"[{timestamp()}] Template range: {start}-{end}\n")

# 6) Write alignment.ali
logfh.write(f"[{timestamp()}] Writing alignment.ali\n")
with open("alignment.ali","w") as f:
    f.write(f">P1;{pdb_id}\n")
    f.write(f"structureX:{pdb_id}:{start}:A:{end}:A::::\n")
    f.write(template_seq + "*\n")
    f.write(">P1;target\n")
    f.write(f"sequence:target:1:A:{length}:A::::\n")
    f.write(target_seq + "*\n")

# 7) Run Modeller automodel
logfh.write(f"[{timestamp()}] Running Modeller automodel\n")
env = Environ()
env.io.atom_files_directory = ['.']
log.verbose()
a = automodel(env,
    alnfile='alignment.ali',
    knowns=pdb_id,
    sequence='target',
    assess_methods=(assess.DOPE, assess.GA341)
)
a.starting_model = a.ending_model = 1
with redirect_stdout(logfh), redirect_stderr(logfh):
    a.make()

# 8) Insert CRYST1 back into model
logfh.write(f"[{timestamp()}] Reinserting CRYST1\n")
for fn in os.listdir("."):
    if fn.startswith("target") and fn.endswith(".pdb"):
        outfn = fn.replace(".pdb","_with_cryst.pdb")
        with open(fn) as fin, open(outfn,"w") as fout:
            inserted = False
            for L in fin:
                if not inserted and (L.startswith("ATOM") or L.startswith("HETATM")):
                    fout.write(cryst1)
                    inserted = True
                if not L.startswith("CRYST1"):
                    fout.write(L)

logfh.write(f"[{timestamp()}] Pipeline complete\n")
logfh.close()

print(f"\nPipeline finished; see '{outdir}/pipeline.log' for full details.")  
