#!/usr/bin/env python3
"""
pipeline.py

Full Pipeline: Modeller Remodeling with UniProt Sequence & CRYST1 Preservation
"""

import os
import sys
import glob
import shutil
import subprocess
from datetime import datetime
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.Data.IUPACData import protein_letters_3to1 as P3
from Bio import SeqIO

# ───── Configuration ─────
PDB_ID     = "4bgq"
UNIPROT_ID = "O76039"
FULL_LEN   = 303
BASEDIR    = os.getcwd()
WORKDIR    = os.path.join(BASEDIR, "work")
LOGFILE    = os.path.join(BASEDIR, "pipeline.log")

# collect final model paths
FINAL_MODELS = []

def log(msg: str):
    ts = datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC")
    line = f"[{ts}] {msg}"
    print(line)
    with open(LOGFILE, "a") as f:
        f.write(line + "\n")


def run(cmd: str):
    log(f"→ {cmd}")
    res = subprocess.run(cmd, shell=True)
    if res.returncode != 0:
        log(f"✖ Command failed: {cmd}")
        sys.exit(1)


def step2_download():
    log("Step 2: Downloading PDB & UniProt")
    os.makedirs(WORKDIR, exist_ok=True)
    os.chdir(WORKDIR)
    run(f"wget -q https://files.rcsb.org/download/{PDB_ID}.pdb")
    run(f"wget -q https://www.uniprot.org/uniprot/{UNIPROT_ID}.fasta -O {UNIPROT_ID}.fasta")
    os.rename(f"{PDB_ID}.pdb", f"{PDB_ID}_orig.pdb")


def step3_clean_pdb():
    log("Step 3: Cleaning PDB (remove HETATM, keep CRYST1)")
    parser = PDBParser(QUIET=True)
    orig = f"{PDB_ID}_orig.pdb"
    clean = f"{PDB_ID}_clean.pdb"
    cryst1 = ''
    for L in open(orig):
        if L.startswith("CRYST1"):
            cryst1 = L
            break
    class StdSelect(Select):
        def accept_residue(self, r): return r.id[0] == ' '
    structure = parser.get_structure(PDB_ID, orig)
    io = PDBIO(); io.set_structure(structure)
    io.save(clean, select=StdSelect())
    log(f"Cleaned → {clean}")
    shutil.copy(clean, f"{PDB_ID}.pdb")
    open("cryst1.txt","w").write(cryst1)


def step4_seqres_vs_uniprot():
    log("Step 4: Reconciling SEQRES vs UniProt")
    orig = f"{PDB_ID}_orig.pdb"
    seqres = ''
    for L in open(orig):
        if L.startswith("SEQRES") and L[11]=='A':
            seqres += ''.join(P3.get(r.capitalize(),'X') for r in L.split()[4:])
    if len(seqres)==FULL_LEN+1 and seqres[0]!=seqres[1]:
        seqres = seqres[1:]
    uni = str(SeqIO.read(f"{UNIPROT_ID}.fasta",'fasta').seq)[:FULL_LEN]
    if seqres != uni:
        log(f"⚠️ mismatch; using UniProt 1–{FULL_LEN}")
        final = uni
    else:
        log("✅ SEQRES matches UniProt")
        final = seqres
    open("final_seq.txt","w").write(final)
    log(f"Final seq length: {len(final)}")


def step5_build_template():
    log("Step 5: Building template sequence from ATOM")
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(PDB_ID, f"{PDB_ID}.pdb")
    chain = next(structure[0].get_chains())
    obs = {r.id[1]:P3.get(r.get_resname().capitalize(),'X') for r in chain if r.id[0]==' '}
    if min(obs)==0:
        obs = {k+1:v for k,v in obs.items()}
    start,end = min(obs), max(obs)
    templ = ''.join(obs.get(i,'-') for i in range(1,FULL_LEN+1))
    open("template_seq.txt","w").write(templ)
    open("range.txt","w").write(f"{start} {end}")
    log(f"Observed residues {start}–{end}")


def step6_write_alignment():
    log("Step 6: Writing alignment.ali")
    start,end = map(int, open("range.txt").read().split())
    templ = open("template_seq.txt").read().strip()
    final = open("final_seq.txt").read().strip()
    with open("alignment.ali","w") as f:
        f.write(f">P1;{PDB_ID}\n")
        f.write(f"structureX:{PDB_ID}:{start}:A:{end}:A::::\n")
        f.write(templ+"*\n")
        f.write(">P1;target\n")
        f.write(f"sequence:target:1:A:{FULL_LEN}:A::::\n")
        f.write(final+"*\n")
    log("✅ alignment.ali written")


def step7_run_modeller():
    log("Step 7: Running Modeller")
    with open("run_modeller.py","w") as f:
        f.write(f"""from modeller import environ
from modeller.automodel import automodel, assess
env = environ()
env.io.hetatm=True
a = automodel(env,alnfile='alignment.ali',knowns='{PDB_ID}',sequence='target',assess_methods=(assess.DOPE,assess.GA341))
a.starting_model=1; a.ending_model=1; a.make()
""" )
    run("mod10.7 run_modeller.py")


def step8_insert_cryst1():
    log("Step 8: Inserting CRYST1")
    cryst = open("cryst1.txt").read()
    for fn in glob.glob("target*.pdb") + glob.glob(f"{PDB_ID}.B*.pdb"):
        out = fn.replace(".pdb","_with_cryst.pdb")
        with open(fn) as inp, open(out,"w") as outp:
            inserted=False
            for L in inp:
                if not inserted and (L.startswith("ATOM") or L.startswith("HETATM")):
                    outp.write(cryst)
                    inserted=True
                if not L.startswith("CRYST1"):
                    outp.write(L)
            if not inserted:
                outp.write(cryst)
        FINAL_MODELS.append(os.path.join(WORKDIR,out))
        log(f"Wrote {out}")


def main():
    open(LOGFILE,"w").write(f"Pipeline started at {datetime.utcnow()}\n")
    step2_download()
    step3_clean_pdb()
    step4_seqres_vs_uniprot()
    step5_build_template()
    step6_write_alignment()
    step7_run_modeller()
    step8_insert_cryst1()
    log("Pipeline complete")

    # Summary
    print("\n===== Workspace Summary =====")
    print(f"Workspace directory: {WORKDIR}\n")
    print("Files:")
    for fn in [f"{PDB_ID}_orig.pdb", f"{UNIPROT_ID}.fasta", f"{PDB_ID}.pdb", "final_seq.txt", "template_seq.txt", "alignment.ali"]:
        print(f"  {os.path.join(WORKDIR, fn)}")
    print("Final Models:")
    for m in FINAL_MODELS:
        print(f"  {m}")
    print(f"Log file: {LOGFILE}")

if __name__=="__main__":
    main()
