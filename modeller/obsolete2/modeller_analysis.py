#!/usr/bin/env python3
"""
modeller_analysis.py

Analyze and summarize your Modeller output, using the truncated model.

Usage:
    # Basic (auto-detects template and truncated model in <dir>/)
    python3 modeller_analysis.py --dir 4ldj

    # Or specify all files explicitly:
    python3 modeller_analysis.py \
        --dir      4ldj \
        --template 4ldj/4ldj.pdb \
        --model    4ldj/4ldj_truncated_9-170.pdb \
        --profile  4ldj/target.V99990001

This script will:
 1) Parse pipeline.log for molpdf, DOPE, GA341
 2) Read alignment.ali to compute sequence identity & coverage
 3) Superimpose CA atoms between template and model, computing RMSD
 4) Plot the per-residue DOPE profile and save as dope_profile.png
"""

import os
import argparse
from Bio.PDB import PDBParser, Superimposer
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt

def parse_pipeline_log(log_path):
    metrics = {}
    with open(log_path) as f:
        for line in f:
            if line.startswith("target.B"):
                parts = line.split()
                metrics['molpdf'] = float(parts[1])
                metrics['DOPE']   = float(parts[2])
                metrics['GA341']  = float(parts[3])
                break
    return metrics

def compute_sequence_identity(aln_path):
    seqs = []
    with open(aln_path) as f:
        for line in f:
            if not line.startswith('>'):
                seqs.append(line.strip().rstrip('*'))
    matches = sum(a == b for a, b in zip(seqs[0], seqs[1]))
    identity = matches / len(seqs[1]) * 100
    coverage = len(seqs[0].replace('-', '')) / len(seqs[1]) * 100
    return identity, coverage

def compute_rmsd(template_pdb, model_pdb):
    parser = PDBParser(QUIET=True)
    tpl = parser.get_structure('tpl', template_pdb)[0]
    mdl = parser.get_structure('mdl', model_pdb)[0]

    tpl_ca = {(c.id, r.id[1]): r['CA'] for c in tpl for r in c if 'CA' in r}
    mdl_ca = {(c.id, r.id[1]): r['CA'] for c in mdl for r in c if 'CA' in r}

    common = sorted(set(tpl_ca) & set(mdl_ca))
    if not common:
        raise ValueError("No overlapping CA residues to superimpose")

    fixed  = [tpl_ca[k] for k in common]
    moving = [mdl_ca[k] for k in common]

    sup = Superimposer()
    sup.set_atoms(fixed, moving)
    return sup.rms

def plot_dope(profile_path, out_png):
    df = pd.read_csv(profile_path, delim_whitespace=True, comment='#',
                     names=['residue','name','score'])
    plt.figure()
    plt.plot(df['residue'], df['score'], marker='.')
    plt.xlabel('Residue Number')
    plt.ylabel('DOPE Score')
    plt.title('Per-Residue DOPE Profile')
    plt.axhline(0, linestyle='--', color='gray')
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Summarize Modeller results")
    parser.add_argument('--dir',      required=True, help='Model directory')
    parser.add_argument('--template', help='Template PDB (default: <dir>/<dir>.pdb)')
    parser.add_argument('--model',    help='Model PDB (default: first <dir>_truncated_*.pdb)')
    parser.add_argument('--profile',  help='DOPE profile file (default: target.V99990001)')
    args = parser.parse_args()
    d = args.dir

    # Default paths
    log_path  = os.path.join(d, 'pipeline.log')
    aln_path  = os.path.join(d, 'alignment.ali')
    tpl_pdb   = args.template or os.path.join(d, f"{d}.pdb")
    model_pdb = args.model or os.path.join(
        d,
        next(fn for fn in os.listdir(d)
             if fn.startswith(f"{d}_truncated_") and fn.endswith(".pdb"))
    )
    profile   = args.profile or os.path.join(d, 'target.V99990001')

    # Compute metrics
    metrics  = parse_pipeline_log(log_path)
    identity, coverage = compute_sequence_identity(aln_path)
    rmsd     = compute_rmsd(tpl_pdb, model_pdb)

    # Print summary
    print("Global Quality Metrics:")
    print(f"  molpdf: {metrics['molpdf']}")
    print(f"   DOPE: {metrics['DOPE']}")
    print(f"  GA341: {metrics['GA341']:.4f}")
    print(f"Sequence identity: {identity:.2f}%")
    print(f"       coverage: {coverage:.2f}%")
    print(f"          RMSD: {rmsd:.2f} Å")

    # Plot DOPE profile
    out_png = os.path.join(d, 'dope_profile.png')
    plot_dope(profile, out_png)
    print(f"DOPE profile plot saved to {out_png}")

if __name__ == '__main__':
    main()
