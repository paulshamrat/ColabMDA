#!/usr/bin/env python3
"""
modeller6.py — Build & Mutate pipeline for Modeller

Two modes:
1) Build mode (default): PDB_ID + UniProt_ID → align2d → automodel → CRYST1 reinsertion
   (optional --range/--length to set UniProt region; optional --truncate; optional mutations)
2) Mutate-only: mutate an existing PDB with one or many mutations

Examples
========
# Help
  python3 modeller6.py --help

# Build model for chain A, UniProt 1–303
  python3 modeller6.py 4bgq O76039 --chain A --range 1 303

# Build + one mutation
  python3 modeller6.py 4bgq O76039 --chain A --range 1 303 --mut K76E

# Build + batch mutations saved under outdir
  python3 modeller6.py 4bgq O76039 --chain A --range 1 303 --list mutations.txt --outdir-mut 4bgq_fix/mutants

# Mutate-only (single)
  python3 modeller6.py --pdb-in 4bgq_fix/target.B99990001_with_cryst.pdb --chain A --mut K76E

# Mutate-only (batch)
  python3 modeller6.py --pdb-in 4bgq_fix/target.B99990001_with_cryst.pdb --chain A --list mutations.txt --outdir-mut 4bgq_fix/mutants --seed 123
"""

import argparse
import csv
import os
import re
import sys
import urllib.request
from contextlib import redirect_stderr, redirect_stdout
from datetime import datetime
from textwrap import dedent

from Bio import SeqIO
from Bio.Data.IUPACData import protein_letters_3to1 as P3
from Bio.PDB import PDBParser

# Modeller
try:
    from modeller import Alignment, Environ, Model, Selection
    from modeller import log as mlog
    from modeller.automodel import assess, automodel, autosched
    from modeller.optimizers import ConjugateGradients, MolecularDynamics

    _MODELLER_OK = True
except Exception as e:
    _MODELLER_OK = False
    _MODELLER_ERR = e


# ───────────── utilities ─────────────


def ts():
    return datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC")


def write(logfh, msg):
    line = f"[{ts()}] {msg}\n"
    logfh.write(line)
    logfh.flush()
    print(line, end="")


def download(url, path, logfh):
    write(logfh, f"Downloading {url} → {path}")
    urllib.request.urlretrieve(url, path)


def clean_pdb_keep_cryst1(orig_path, clean_path, logfh):
    """Write CRYST1 + ATOM lines to cleaned PDB; return CRYST1 line (if any)."""
    cryst1 = ""
    with open(orig_path) as fin, open(clean_path, "w") as fout:
        for L in fin:
            if L.startswith("CRYST1"):
                cryst1 = L
                fout.write(L)
            elif L.startswith("ATOM"):
                fout.write(L)
    write(logfh, "CRYST1 found and preserved." if cryst1 else "No CRYST1 found in original PDB.")
    return cryst1


def extract_seqres_chain(orig_pdb, chain_id, logfh):
    seq = []
    with open(orig_pdb) as f:
        for line in f:
            if line.startswith("SEQRES") and line[11] == (chain_id if chain_id else " "):
                for res3 in line.split()[4:]:
                    seq.append(P3.get(res3.capitalize(), "X"))
    write(logfh, f"SEQRES({chain_id}) length = {len(seq)}")
    return "".join(seq)


def build_template_from_atoms(clean_pdb, chain_id, logfh):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("tmpl", clean_pdb)
    model0 = list(structure)[0]

    chain = None
    for ch in model0:
        if ch.id.strip() == chain_id:
            chain = ch
            break
    if chain is None:
        chain = list(model0)[0]
        chain_id = chain.id
        write(logfh, f"Requested chain not found; using first chain '{chain_id}'.")

    pos2aa = {}
    for res in chain:
        if res.id[0] != " " or "CA" not in res:
            continue
        aa = P3.get(res.get_resname().strip().capitalize(), "X")
        pos2aa[res.id[1]] = aa

    if not pos2aa:
        raise RuntimeError("No standard residues with CA found in chosen chain.")

    start, end = min(pos2aa), max(pos2aa)
    template_seq = "".join(pos2aa.get(i, "-") for i in range(start, end + 1))
    write(
        logfh,
        f"Template observed chain {chain_id}: residues {start}-{end} (len={len(template_seq)})",
    )
    return chain_id, start, end, template_seq


def write_target_pir(target_seq, start_res=1, chain_id="A", out_path="target.ali"):
    end_res = start_res + len(target_seq) - 1
    with open(out_path, "w") as f:
        f.write(">P1;target\n")
        f.write(f"sequence:target:{start_res}:{chain_id}:{end_res}:{chain_id}::::\n")
        f.write(target_seq + "*\n")


def run_align2d_write_pir(pdb_id, out_pir, logfh):
    if not _MODELLER_OK:
        raise ImportError(f"Modeller not available: {_MODELLER_ERR}")
    env = Environ()
    aln = Alignment(env)
    mdl = Model(env, file=pdb_id, model_segment=("FIRST:@", "LAST:@"))  # keeps real numbering
    aln.append_model(mdl, atom_files=pdb_id, align_codes=pdb_id)
    aln.append(file="target.ali", align_codes="target")
    aln.align2d(overhang=50, gap_penalties_1d=(-450, -50))
    aln.write(file=out_pir, alignment_format="PIR")
    write(logfh, f"Wrote {out_pir} via align2d()")


def run_automodel(alnfile, knowns, sequence, logfh, verbose=False):
    if not _MODELLER_OK:
        raise ImportError(f"Modeller not available: {_MODELLER_ERR}")
    env = Environ()
    env.io.hetatm = True
    mlog.verbose() if verbose else mlog.minimal()
    a = automodel(
        env,
        alnfile=alnfile,
        knowns=knowns,
        sequence=sequence,
        assess_methods=(assess.DOPE, assess.GA341),
    )
    a.starting_model = 1
    a.ending_model = 1
    write(logfh, "Running Modeller automodel (1 model)...")
    with redirect_stdout(logfh), redirect_stderr(logfh):
        a.make()


def find_first_model_pdb(cwd):
    for fn in sorted(os.listdir(cwd)):
        if fn.startswith("target.") and fn.endswith(".pdb"):
            return fn
    return None


def reinsert_cryst1(src_pdb, out_pdb, cryst1_line):
    inserted = False
    with open(src_pdb) as fin, open(out_pdb, "w") as fout:
        for L in fin:
            if not inserted and (L.startswith("ATOM") or L.startswith("HETATM")):
                if cryst1_line:
                    fout.write(cryst1_line)
                inserted = True
                fout.write(L)
            elif L.startswith("CRYST1"):
                continue
            else:
                fout.write(L)
        if not inserted and cryst1_line:
            fout.write(cryst1_line)


def truncate_pdb_range(in_pdb, chain_id, start, end, out_pdb, cryst1_line=None):
    with open(in_pdb) as fin:
        lines = fin.readlines()
    with open(out_pdb, "w") as fout:
        inserted = False
        for L in lines:
            if L.startswith(("ATOM", "HETATM")):
                ch = L[21].strip() or " "
                try:
                    resseq = int(L[22:26])
                except Exception:
                    continue
                if ch == chain_id and (start <= resseq <= end):
                    if cryst1_line and not inserted:
                        fout.write(cryst1_line)
                        inserted = True
                    fout.write(L)
            else:
                if not L.startswith("CRYST1"):
                    fout.write(L)
        if cryst1_line and not inserted:
            fout.write(cryst1_line)


def print_alignment_summary(aln_path, logfh):
    """Parses a PIR alignment file and prints a summary to the screen."""
    try:
        with open(aln_path) as f:
            lines = f.readlines()

        entries = []
        current = None
        for L in lines:
            if L.startswith(">P1;"):
                if current:
                    entries.append(current)
                current = {"id": L[4:].strip(), "seq": ""}
            elif current and not L.startswith(("structure", "sequence")):
                current["seq"] += L.strip().replace("*", "")
        if current:
            entries.append(current)

        if len(entries) < 2:
            return

        t1, t2 = entries[0], entries[1]
        matches = sum(1 for a, b in zip(t1["seq"], t2["seq"]) if a == b and a != "-")
        total = sum(1 for a, b in zip(t1["seq"], t2["seq"]) if a != "-" or b != "-")
        identity = (matches / total * 100) if total > 0 else 0

        write(logfh, "--- Alignment Summary ---")
        write(logfh, f"Template: {t1['id']}")
        write(logfh, f"Target:   {t2['id']}")
        write(logfh, f"Identity: {identity:.1f}%")

        # Show full alignment in 60-residue blocks
        write(logfh, "\nFull Alignment (Template vs Target):")
        full1, full2 = t1["seq"], t2["seq"]
        for start in range(0, len(full1), 60):
            end = min(start + 60, len(full1))
            s1_chunk = full1[start:end]
            s2_chunk = full2[start:end]
            m_chunk = "".join(
                "*" if a == b and a != "-" else " " for a, b in zip(s1_chunk, s2_chunk)
            )

            write(logfh, f"{'Template:':<12} {s1_chunk}")
            write(logfh, f"{'Target:':<12} {s2_chunk}")
            write(logfh, f"{'Match:':<12} {m_chunk}")
            write(logfh, "")
    except Exception as e:
        write(logfh, f"Warning: Could not generate alignment summary: {e}")


def align_pdb_to_template_in_place(model_path, template_path, logfh):
    """Aligns model_path PDB to template_path PDB using Bio.PDB.Superimposer with sequence-based matching."""
    try:
        from Bio.PDB import PDBParser, Superimposer, PDBIO
        from Bio.Align import PairwiseAligner

        parser = PDBParser(QUIET=True)
        model_struct = parser.get_structure("model", model_path)
        template_struct = parser.get_structure("template", template_path)

        three_to_one = {
            "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
            "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
            "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
            "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
        }

        # Extract CA atoms and sequence
        model_atoms = []
        model_seq = []
        for model in model_struct:
            for chain in model:
                for residue in chain:
                    if residue.id[0] == " " and "CA" in residue:
                        resname = residue.get_resname().strip().upper()
                        model_seq.append(three_to_one.get(resname, "X"))
                        model_atoms.append(residue["CA"])

        template_atoms = []
        template_seq = []
        for model in template_struct:
            for chain in model:
                for residue in chain:
                    if residue.id[0] == " " and "CA" in residue:
                        resname = residue.get_resname().strip().upper()
                        template_seq.append(three_to_one.get(resname, "X"))
                        template_atoms.append(residue["CA"])

        aligner = PairwiseAligner()
        aligner.mode = 'global'
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -1

        alignments = aligner.align("".join(model_seq), "".join(template_seq))
        if not alignments:
            write(logfh, "Warning: Could not align sequences for coordinate alignment.")
            return

        best_aln = alignments[0]
        aligned_indices = best_aln.aligned

        matched_model_atoms = []
        matched_template_atoms = []

        for seg_mod, seg_tmpl in zip(aligned_indices[0], aligned_indices[1]):
            for m_idx, t_idx in zip(range(seg_mod[0], seg_mod[1]), range(seg_tmpl[0], seg_tmpl[1])):
                matched_model_atoms.append(model_atoms[m_idx])
                matched_template_atoms.append(template_atoms[t_idx])

        if not matched_model_atoms:
            write(logfh, "Warning: No matching residues found for coordinate alignment.")
            return

        superimposer = Superimposer()
        superimposer.set_atoms(matched_template_atoms, matched_model_atoms)
        superimposer.apply(model_struct.get_atoms())

        write(logfh, f"Aligned model to template structure. RMSD = {superimposer.rms:.4f} Angstroms")

        io = PDBIO()
        io.set_structure(model_struct)
        io.save(model_path)
    except Exception as e:
        write(logfh, f"Warning: Failed to align model to template coordinate frame: {e}")


def renumber_pdb(pdb_path, start_res, logfh):
    """
    Physically changes residue numbers in a PDB file to start from 'start_res'.
    """
    try:
        with open(pdb_path) as f:
            lines = f.readlines()

        current_min = None
        for L in lines:
            if L.startswith(("ATOM", "HETATM")):
                try:
                    resid = int(L[22:26])
                    if current_min is None or resid < current_min:
                        current_min = resid
                except Exception:
                    continue

        if current_min is None:
            return
        offset = start_res - current_min
        if offset == 0:
            return

        write(logfh, f"Applying numbering offset: {current_min} -> {start_res} (offset={offset})")
        with open(pdb_path, "w") as f:
            for L in lines:
                if L.startswith(("ATOM", "HETATM")):
                    try:
                        old_id = int(L[22:26])
                        new_id = old_id + offset
                        new_id_str = f"{new_id:4d}"
                        L = L[:22] + new_id_str + L[26:]
                    except Exception:
                        pass
                f.write(L)
    except Exception as e:
        write(logfh, f"Warning: Failed to renumber PDB: {e}")


def verify_model_sequence(pdb_path, chain_id, expected_seq, expected_start, logfh):
    """Sanity check: verifies model residue types and numbering against UniProt."""
    try:
        # Extract what we actually got in the PDB
        obs_chain, obs_start, obs_end, obs_seq = build_template_from_atoms(
            pdb_path, chain_id, logfh
        )

        expected_end = expected_start + len(expected_seq) - 1
        write(logfh, "--- Post-Build Sanity Check ---")
        write(logfh, f"Requested Range: {expected_start} to {expected_end}")
        write(logfh, f"Observed Range:  {obs_start} to {obs_end}")

        status = "OK"
        if obs_start != expected_start or obs_end != expected_end:
            write(
                logfh,
                f"WARNING: Range mismatch! Expected {expected_start}-{expected_end}, got {obs_start}-{obs_end}",
            )
            status = "WARN"

        # Compare sequences and show a snippet
        mismatches = 0
        comp_len = min(len(obs_seq), len(expected_seq))
        match_str = ""
        for i in range(comp_len):
            o, e = obs_seq[i], expected_seq[i]
            if o == e and o != "-":
                match_str += "*"
            elif o == "-" or e == "-":
                match_str += " "
            else:
                match_str += "!"
                mismatches += 1

        # Show full sequence comparison in 60-residue blocks
        write(logfh, "\nFinal PDB vs UniProt Sequence Comparison:")
        for start in range(0, len(expected_seq), 60):
            end = min(start + 60, len(expected_seq))
            pdb_idx_label = f"PDB Index {obs_start + start}:"
            ref_label = "UniProt Ref:"
            match_label = "Match Map:"

            # Align labels to 18 characters
            write(logfh, f"{pdb_idx_label:<18} {obs_seq[start:end]}")
            write(logfh, f"{ref_label:<18} {expected_seq[start:end]}")
            write(logfh, f"{match_label:<18} {match_str[start:end]}")
            write(logfh, "")

        write(logfh, "(* = match, ! = mismatch, space = gap)\n")

        if mismatches > 0:
            write(logfh, f"ERROR: Found {mismatches} amino acid mismatches!")
            status = "FAIL"

        if status == "OK":
            write(logfh, "✅ SUCCESS: Model numbering and sequence verified.")
        elif status == "WARN":
            write(logfh, "⚠️  COMPLETED: Model built but range is slightly different.")
        else:
            write(logfh, "❌ FAILED: Significant sequence errors detected.")

    except Exception as e:
        write(logfh, f"Warning: Sanity check encountered an error: {e}")


# ───────────── mutation helpers ─────────────


def one_to_three(aa):
    m = {
        "A": "ALA",
        "R": "ARG",
        "N": "ASN",
        "D": "ASP",
        "C": "CYS",
        "E": "GLU",
        "Q": "GLN",
        "G": "GLY",
        "H": "HIS",
        "I": "ILE",
        "L": "LEU",
        "K": "LYS",
        "M": "MET",
        "F": "PHE",
        "P": "PRO",
        "S": "SER",
        "T": "THR",
        "W": "TRP",
        "Y": "TYR",
        "V": "VAL",
    }
    aa = aa.upper()
    if aa not in m:
        raise ValueError(f"Unsupported residue: {aa}")
    return m[aa]


def parse_mut(s):
    s = s.strip()
    if not re.match(r"^[A-Za-z]\d+[A-Za-z]$", s):
        raise ValueError(f"Bad mutation format: {s} (expected F13S)")
    return s[0].upper(), int(re.findall(r"\d+", s)[0]), s[-1].upper()


def read_mut_list(path):
    muts = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            muts.append(line)
    return muts


def read_cryst1(pdb_path):
    with open(pdb_path) as f:
        for L in f:
            if L.startswith("CRYST1"):
                return L
    return ""


def observed_wt_letter(pdb_path, chain_id, pos):
    parser = PDBParser(QUIET=True)
    st = parser.get_structure("base", pdb_path)
    mdl = list(st)[0]
    chain = None
    for ch in mdl:
        if ch.id.strip() == chain_id:
            chain = ch
            break
    if chain is None:
        chain = list(mdl)[0]
    for res in chain:
        if res.id[0] != " ":
            continue
        if res.id[1] == pos:
            return P3.get(res.get_resname().strip().capitalize(), "X")
    return None


def make_restraints(mdl, aln):
    rsr = mdl.restraints
    rsr.clear()
    s = Selection(mdl)
    for typ in ("stereo", "phi-psi_binormal"):
        rsr.make(s, restraint_type=typ, aln=aln, spline_on_site=True)
    for typ in ("omega", "chi1", "chi2", "chi3", "chi4"):
        rsr.make(
            s,
            restraint_type=typ + "_dihedral",
            spline_range=4.0,
            spline_dx=0.3,
            spline_min_points=5,
            aln=aln,
            spline_on_site=True,
        )


def refine_md(atmsel):
    md = MolecularDynamics(cap_atom_shift=0.39, md_time_step=4.0, md_return="FINAL")
    init_vel = True
    for its, equil, temps in (
        (200, 20, (150.0, 250.0, 400.0, 700.0, 1000.0)),
        (200, 600, (1000.0, 800.0, 600.0, 500.0, 400.0, 300.0)),
    ):
        for T in temps:
            md.optimize(
                atmsel,
                init_velocities=init_vel,
                temperature=T,
                max_iterations=its,
                equilibrate=equil,
            )
            init_vel = False


def optimize_sel(atmsel, mdl):
    sched = autosched.loop.make_for_model(mdl)
    for step in sched:
        step.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)
    refine_md(atmsel)
    cg = ConjugateGradients()
    cg.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)


def mutate_once(pdb_in, chain, mut_str, seed, outdir_mut=None, verbose=False):
    """
    Perform one mutation on 'pdb_in' and return:
      (final_out_path, E_unopt, E_opt, observed_WT, status)

    Official sequence: mutate → generate_topology → transfer_xyz → build → (res_num_from)
                       → make_restraints → energies/opt → write → reinsert CRYST1

    All artifacts (core .pdb, *_with_cryst.pdb, .tmp) are written into outdir_mut if provided,
    otherwise next to the input PDB.
    """
    if not _MODELLER_OK:
        raise ImportError(f"Modeller not available: {_MODELLER_ERR}")

    wt, pos, to = parse_mut(mut_str)
    cryst1 = read_cryst1(pdb_in)

    # Check WT residue in the input structure (PDB numbering)
    obs = observed_wt_letter(pdb_in, chain, pos)
    if obs is None:
        return None, None, None, obs, f"FAIL: residue {chain}:{pos} not found"
    if obs != wt:
        return (
            None,
            None,
            None,
            obs,
            f"FAIL: WT mismatch at {chain}:{pos} (expected {wt}, found {obs})",
        )

    basedir, basefn = os.path.split(pdb_in)
    stem = os.path.splitext(basefn)[0]
    if stem.endswith("_with_cryst"):
        stem = stem[: -len("_with_cryst")]

    # Decide where to place ALL artifacts
    if outdir_mut:
        outdir_mut_abs = os.path.abspath(outdir_mut)
        os.makedirs(outdir_mut_abs, exist_ok=True)
        working_dir = outdir_mut_abs
    else:
        working_dir = basedir if basedir else os.getcwd()

    core_name = f"{stem}_{mut_str}.pdb"
    final_name = f"{stem}_{mut_str}_with_cryst.pdb"
    tmp_name = f"{stem}_{mut_str}.tmp"

    core_abs = os.path.join(working_dir, core_name)
    final_abs = os.path.join(working_dir, final_name)
    tmp_abs = os.path.join(working_dir, tmp_name)

    cwd = os.getcwd()
    if basedir:
        os.chdir(basedir)
    try:
        env = Environ(rand_seed=seed)
        env.io.hetatm = True
        env.edat.dynamic_sphere = False
        env.edat.dynamic_lennard = True
        env.edat.contact_shell = 4.0
        env.edat.update_dynamic = 0.39
        env.libs.topology.read(file="$(LIB)/top_heav.lib")
        env.libs.parameters.read(file="$(LIB)/par.lib")

        mlog.verbose() if verbose else mlog.minimal()

        # Load model and append once to alignment
        mdl1 = Model(env, file=os.path.basename(basefn))
        ali = Alignment(env)
        ali.append_model(
            mdl1, atom_files=os.path.basename(basefn), align_codes=os.path.basename(basefn)
        )

        # Mutate side chain at requested residue
        try:
            s = Selection(mdl1.chains[chain].residues[str(pos)])
        except Exception:
            s = Selection(mdl1.chains[chain].residues[pos])
        s.mutate(residue_type=one_to_three(to))

        # Append mutated model to alignment (Modeller trick)
        ali.append_model(mdl1, align_codes=os.path.basename(basefn))

        # Generate topology & build coordinates for mutant
        mdl1.clear_topology()
        mdl1.generate_topology(ali[-1])
        mdl1.transfer_xyz(ali)
        mdl1.build(initialize_xyz=False, build_method="INTERNAL_COORDINATES")

        # Keep original residue numbering (optional but recommended)
        mdl2 = Model(env, file=os.path.basename(basefn))
        mdl1.res_num_from(mdl2, ali)

        # Write/read tmp to refresh all internal tables — write to working_dir
        mdl1.write(file=tmp_abs)
        mdl1.read(file=tmp_abs)

        # Now make restraints and compute energies/optimize
        make_restraints(mdl1, ali)
        try:
            s = Selection(mdl1.chains[chain].residues[str(pos)])
        except Exception:
            s = Selection(mdl1.chains[chain].residues[pos])

        mdl1.restraints.unpick_all()
        mdl1.restraints.pick(s)

        mdl1.env.edat.nonbonded_sel_atoms = 1
        try:
            E_unopt = s.energy()
        except Exception:
            E_unopt = None

        s.randomize_xyz(deviation=4.0)
        mdl1.env.edat.nonbonded_sel_atoms = 2
        optimize_sel(s, mdl1)
        mdl1.env.edat.nonbonded_sel_atoms = 1
        optimize_sel(s, mdl1)

        try:
            E_opt = s.energy()
        except Exception:
            E_opt = None

        # Write mutant core and final WITH CRYST1 directly into working_dir
        mdl1.write(file=core_abs)
        reinsert_cryst1(core_abs, final_abs, cryst1)

        # Cleanup tmp (in working_dir)
        try:
            os.remove(tmp_abs)
        except Exception:
            pass

        return final_abs, E_unopt, E_opt, obs, "OK"

    finally:
        if basedir:
            os.chdir(cwd)


# ───────────── CLI ─────────────


def main():
    epilog = dedent("""
    Modes
    -----
    • Build mode (default): provide PDB ID and UniProt ID (positional args).
    • Mutate-only mode: provide --pdb-in to mutate an existing PDB; build steps are skipped.

    Notes
    -----
    • --chain applies to template extraction, truncation, and mutations (default A).
    • Mutations use PDB residue numbering of the input structure.
    • If you use --range, truncation is usually unnecessary because the model matches that subsequence.
    """)

    ap = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="Build: template PDB + UniProt → align2d → automodel → CRYST1 → (optional) truncate → (optional) mutate\n"
        "Mutate-only: mutate an existing PDB directly.",
        epilog=epilog,
    )

    # Positional for build mode
    ap.add_argument("pdb_id", nargs="?", help="Template PDB ID (e.g., 4bgq)")
    ap.add_argument("uni_id", nargs="?", help="UniProt ID (e.g., O76039)")

    # Mutate-only toggle
    ap.add_argument("--pdb-in", dest="pdb_in", help="Mutate-only mode: input PDB to mutate")

    # Build options
    ap.add_argument(
        "--outdir", default=None, help="Build output directory (default: <pdb_id>_fix/)"
    )
    ap.add_argument(
        "--chain", default="A", help="Chain ID for template/trim/mutations (default: A)"
    )
    ap.add_argument(
        "--length", type=int, default=None, help="Model UniProt 1..L (e.g., --length 303)"
    )
    ap.add_argument(
        "--range",
        nargs=2,
        type=int,
        metavar=("START", "END"),
        help="Model UniProt START..END (1-based inclusive), e.g., --range 1 303",
    )
    ap.add_argument(
        "--truncate",
        action="store_true",
        help="After modeling, also trim to --range on --chain (usually not needed if you used --range)",
    )
    ap.add_argument(
        "--uniprot-numbering",
        action="store_true",
        help="Use UniProt indices for PDB numbering (requires --range)",
    )

    # Mutation options
    ap.add_argument("--mut", help="Single mutation like K76E")
    ap.add_argument("--list", help="File with one mutation per line (e.g., F13S)")
    ap.add_argument("--outdir-mut", help="Directory to write mutant PDBs")
    ap.add_argument(
        "--seed", type=int, default=-49837, help="Modeller RNG seed (default deterministic)"
    )

    # Logging
    ap.add_argument(
        "--logfile", default=None, help="Log path (build: pipeline.log; mutate-only: mutate.log)"
    )
    ap.add_argument("--verbose", action="store_true", help="Verbose Modeller logs")

    args = ap.parse_args()

    if not _MODELLER_OK:
        sys.exit(f"ERROR: Modeller not available: {_MODELLER_ERR}")

    # ── Mutate-only mode ──
    if args.pdb_in:
        pdb_in = args.pdb_in
        if not os.path.isabs(pdb_in):
            script_dir = os.path.dirname(os.path.abspath(__file__))
            candidate = os.path.join(script_dir, pdb_in)
            pdb_in = candidate if os.path.exists(candidate) else os.path.abspath(args.pdb_in)
        if not os.path.exists(pdb_in):
            sys.exit(f"ERROR: --pdb-in file not found: {pdb_in}")

        basedir = os.path.dirname(os.path.abspath(pdb_in))
        log_path = args.logfile or os.path.join(basedir, "mutate.log")
        muts = []
        if args.mut:
            muts.append(args.mut)
        if args.list:
            muts.extend(read_mut_list(args.list))
        if not muts:
            sys.exit("ERROR: in mutate-only mode you must provide --mut or --list.")

        created = []
        with (
            open(log_path, "w") as logfh,
            open(os.path.join(basedir, "mutate_summary.csv"), "a", newline="") as cf,
        ):
            write(logfh, "Mutate-only mode started")
            write(logfh, f"Input PDB={pdb_in}, chain={args.chain}, seed={args.seed}")
            w = csv.writer(cf)
            if cf.tell() == 0:
                w.writerow(
                    [
                        "pdb_in",
                        "mutation",
                        "chain",
                        "seed",
                        "out_pdb",
                        "E_unopt",
                        "E_opt",
                        "observed_WT",
                        "status",
                    ]
                )

            outdir_mut = args.outdir_mut
            if outdir_mut:
                os.makedirs(outdir_mut, exist_ok=True)
                write(logfh, f"Mutants will be written to: {os.path.abspath(outdir_mut)}")

            for m in muts:
                try:
                    out_pdb, E_unopt, E_opt, obs, status = mutate_once(
                        pdb_in,
                        args.chain,
                        m,
                        args.seed,
                        outdir_mut=outdir_mut,
                        verbose=args.verbose,
                    )
                    w.writerow(
                        [
                            pdb_in,
                            m,
                            args.chain,
                            args.seed,
                            out_pdb or "",
                            E_unopt,
                            E_opt,
                            obs,
                            status,
                        ]
                    )
                    write(logfh, f"{m} → {status}; out={out_pdb}, E_unopt={E_unopt}, E_opt={E_opt}")
                    if out_pdb:
                        created.append(out_pdb)
                except Exception as e:
                    w.writerow([pdb_in, m, args.chain, args.seed, "", "", "", "", f"FAIL: {e}"])
                    write(logfh, f"{m} → FAIL: {e}")

            write(logfh, "Mutate-only complete")
            print(f"\n✅ Done. Log: {log_path}")
            if created:
                print("• Created mutants:")
                for p in created:
                    print("  -", p)
            print(f"• Summary: {os.path.join(basedir, 'mutate_summary.csv')}")
        return

    # ── Build mode ──
    if not (args.pdb_id and args.uni_id):
        ap.error(
            "Build mode requires positional arguments: pdb_id uni_id (or use --pdb-in for mutate-only mode)"
        )

    pdb_id = args.pdb_id.lower()
    uni_id = args.uni_id.upper()
    outdir = args.outdir or f"{pdb_id}_fix"
    os.makedirs(outdir, exist_ok=True)
    log_path = args.logfile or os.path.join(outdir, "pipeline.log")

    with open(log_path, "w") as logfh:
        write(logfh, "Build mode started")
        write(
            logfh,
            f"Inputs: pdb_id={pdb_id}, uni_id={uni_id}, chain={args.chain}, "
            f"length={args.length}, range={args.range}, truncate={args.truncate}, "
            f"mut={args.mut}, list={args.list}, seed={args.seed}, outdir_mut={args.outdir_mut}",
        )

        os.chdir(outdir)

        # 1) Download inputs
        pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        uni_url = f"https://rest.uniprot.org/uniprotkb/{uni_id}.fasta"
        download(pdb_url, f"{pdb_id}.pdb", logfh)
        download(uni_url, f"{uni_id}.fasta", logfh)
        os.system(f"cp {pdb_id}.pdb {pdb_id}_orig.pdb")

        # 2) Clean PDB + keep CRYST1 (Modeller reads <pdb_id>.pdb)
        cryst1_line = clean_pdb_keep_cryst1(f"{pdb_id}_orig.pdb", f"{pdb_id}_clean.pdb", logfh)
        os.system(f"cp {pdb_id}_clean.pdb {pdb_id}.pdb")

        # 3) Define target UniProt region
        seqres = extract_seqres_chain(f"{pdb_id}_orig.pdb", args.chain, logfh)
        uni_full = str(SeqIO.read(f"{uni_id}.fasta", "fasta").seq)
        if args.range:
            start, end = args.range
            if not (1 <= start <= end <= len(uni_full)):
                raise SystemExit(
                    f"--range {start}-{end} out of bounds for UniProt length {len(uni_full)}"
                )
            target_seq = uni_full[start - 1 : end]
            write(logfh, f"Target UniProt subsequence {start}-{end} (len={len(target_seq)})")
        elif args.length:
            L = args.length
            if L < 1 or L > len(uni_full):
                raise SystemExit(f"--length {L} out of bounds for UniProt length {len(uni_full)}")
            target_seq = uni_full[:L]
            write(logfh, f"Target UniProt 1..{L} (len={len(target_seq)})")
        else:
            target_seq = uni_full
            write(logfh, f"Target UniProt full length {len(target_seq)})")

        if seqres and (len(seqres) != len(target_seq) or seqres != target_seq):
            write(logfh, "SEQRES vs UniProt mismatch; using UniProt region for modeling.")
        elif seqres:
            write(logfh, "SEQRES matches UniProt region; using UniProt region anyway.")

        # 4) Log template coverage
        build_template_from_atoms(f"{pdb_id}.pdb", args.chain, logfh)

        # 5) Align & write PIR
        if args.uniprot_numbering and args.range:
            start_num = args.range[0]
        else:
            start_num = 1

        write_target_pir(
            target_seq, start_res=start_num, chain_id=args.chain, out_path="target.ali"
        )
        run_align2d_write_pir(pdb_id, out_pir="alignment.ali", logfh=logfh)
        print_alignment_summary("alignment.ali", logfh)

        # 6) Build model
        run_automodel(
            alnfile="alignment.ali",
            knowns=pdb_id,
            sequence="target",
            logfh=logfh,
            verbose=args.verbose,
        )

        # 7) Reinsert CRYST1 into result
        model_core = find_first_model_pdb(os.getcwd())
        if not model_core:
            raise SystemExit(
                "ERROR: Modeller did not produce a model PDB (target.*.pdb not found)."
            )
        base_with_cryst = model_core.replace(".pdb", "_with_cryst.pdb")
        reinsert_cryst1(model_core, base_with_cryst, cryst1_line)
        base_with_cryst_abs = os.path.abspath(base_with_cryst)
        write(logfh, f"Base model (CRYST1): {base_with_cryst_abs}")

        # Align model back to template coordinate space to ensure ligand/cofactor alignment
        align_pdb_to_template_in_place(base_with_cryst_abs, f"{pdb_id}_orig.pdb", logfh)

        # Fix Numbering
        renumber_pdb(base_with_cryst_abs, start_num, logfh)

        # Post-build Sanity Check
        verify_model_sequence(base_with_cryst_abs, args.chain, target_seq, start_num, logfh)

        # 8) Optional truncation (not needed if using --range, but available)
        trunc_for_mut = None
        if args.truncate and args.range:
            start, end = args.range
            trunc_core = model_core.replace(".pdb", f"_trunc_{args.chain}_{start}-{end}.pdb")
            truncate_pdb_range(
                base_with_cryst, args.chain, start, end, trunc_core, cryst1_line=cryst1_line
            )
            trunc_for_mut = os.path.abspath(trunc_core)
            write(logfh, f"Truncated model: {trunc_for_mut}")

        base_for_mut = trunc_for_mut or base_with_cryst_abs

        # 9) Optional mutations after build
        muts = []
        if args.mut:
            muts.append(args.mut)
        if args.list:
            muts.extend(read_mut_list(args.list))

        created = []
        if muts:
            outdir_mut = args.outdir_mut or os.getcwd()
            os.makedirs(outdir_mut, exist_ok=True)
            write(logfh, f"Mutants will be written to: {os.path.abspath(outdir_mut)}")
            csv_path = os.path.join(".", "mutate_summary.csv")
            write_header = not os.path.exists(csv_path)
            with open(csv_path, "a", newline="") as cf:
                w = csv.writer(cf)
                if write_header:
                    w.writerow(
                        [
                            "pdb_in",
                            "mutation",
                            "chain",
                            "seed",
                            "out_pdb",
                            "E_unopt",
                            "E_opt",
                            "observed_WT",
                            "status",
                        ]
                    )
                for m in muts:
                    try:
                        out_pdb, E_unopt, E_opt, obs, status = mutate_once(
                            base_for_mut,
                            args.chain,
                            m,
                            args.seed,
                            outdir_mut=outdir_mut,
                            verbose=args.verbose,
                        )
                        w.writerow(
                            [
                                base_for_mut,
                                m,
                                args.chain,
                                args.seed,
                                out_pdb or "",
                                E_unopt,
                                E_opt,
                                obs,
                                status,
                            ]
                        )
                        write(
                            logfh,
                            f"{m} → {status}; out={out_pdb}, E_unopt={E_unopt}, E_opt={E_opt}",
                        )
                        if out_pdb:
                            created.append(out_pdb)
                    except Exception as e:
                        w.writerow(
                            [base_for_mut, m, args.chain, args.seed, "", "", "", "", f"FAIL: {e}"]
                        )
                        write(logfh, f"{m} → FAIL: {e}")

        # Final summary
        print(f"\n✅ Build complete. Log: {log_path}")
        print("• Base model with CRYST1:", base_with_cryst_abs)
        if trunc_for_mut:
            print("• Truncated model      :", trunc_for_mut)
        if created:
            print("• Created mutants:")
            for p in created:
                print("  -", p)
            print("• mutate_summary.csv is in:", os.path.abspath("."))


if __name__ == "__main__":
    main()
