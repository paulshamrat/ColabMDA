"""Small natural-protein tutorials using the public modeling and simulation CLI."""

import hashlib
import io
import json
import os
import re
import subprocess
import sys
from importlib.metadata import version
from pathlib import Path
from urllib.request import urlopen

from Bio.PDB import PDBParser, Select
from Bio.SeqUtils import seq1

from colabmda.modeller.structure import (
    fetch_archive,
    verify_sequence,
    write_chain,
)

PROTEINS = {
    "fkbp": {
        "pdb": "1D6O",
        "uniprot": "P62942",
        "start": 2,
        "end": 108,
        "author_end": 107,
        "variants": ["WT", "F37V"],
        "state": "apo",
    },
    "hras": {
        "pdb": "4Q21",
        "uniprot": "P01112",
        "start": 1,
        "end": 166,
        "author_end": 166,
        "variants": ["WT", "G12V"],
        "state": "gdp_1mg",
    },
}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + "\n")


def cli(*args):
    subprocess.run([sys.executable, "-m", "colabmda.cli", *map(str, args)], check=True)


def fetch(root, name):
    spec = PROTEINS[name]
    root = Path(root).resolve()
    folder = root / "data/01_structure" / name
    folder.mkdir(parents=True, exist_ok=True)
    marker = folder / "STRUCTURE_DONE"
    marker.unlink(missing_ok=True)
    structure, archive, url = fetch_archive(spec["pdb"], folder)
    chain = next(iter(structure))["A"]
    fasta_url = f"https://rest.uniprot.org/uniprotkb/{spec['uniprot']}.fasta"
    fasta = urlopen(fasta_url, timeout=60).read()
    sequence = "".join(fasta.decode().splitlines()[1:])
    target = sequence[spec["start"] - 1 : spec["end"]]
    residues = [r for r in chain if r.id[0] == " " and 1 <= r.id[1] <= spec["author_end"]]
    if (
        len(residues) != len(target)
        or any("CA" not in r for r in residues)
        or "".join(seq1(r.resname) for r in residues) != target
        or [(r.id[1], r.id[2]) for r in residues] != [(i, " ") for i in range(1, len(target) + 1)]
    ):
        raise ValueError("Selected archive coordinates do not match the UniProt region")
    mapping = [
        {
            "chain": "A",
            "author_residue": r.id[1],
            "insertion_code": r.id[2].strip(),
            "uniprot_residue": i,
            "amino_acid": seq1(r.resname),
        }
        for i, r in enumerate(residues, spec["start"])
    ]
    for mutation in spec["variants"][1:]:
        position = int(mutation[1:-1])
        if target[position - spec["start"]] != mutation[0]:
            raise ValueError("Mutation WT identity disagrees with UniProt and coordinates")

    class Protein(Select):
        def accept_residue(self, residue):
            return residue in residues

    # Only heterogens from this author chain; never transfer from a different protein.
    cofactors = [r for r in chain if r.resname in {"GDP", "MG"}] if name == "hras" else []
    if name == "hras" and sorted(r.resname for r in cofactors) != ["GDP", "MG"]:
        raise ValueError("Expected exactly one experimental GDP and one Mg")
    # Keep experimental first-shell Mg waters, recording their original identities.
    mg = next((next(iter(r)) for r in cofactors if r.resname == "MG"), None)
    waters = [
        r for r in chain if r.resname == "HOH" and "O" in r and mg is not None and r["O"] - mg < 2.6
    ]

    class Complex(Select):
        def accept_residue(self, residue):
            return residue in [*residues, *cofactors, *waters]

    files = {"uniprot.fasta": fasta}
    for filename, selection in [("template.pdb", Protein()), ("complex.pdb", Complex())]:
        out = io.StringIO()
        write_chain(chain, out, selection)
        files[filename] = out.getvalue().encode()
    for filename, content in files.items():
        path = folder / filename
        if path.exists() and path.read_bytes() != content:
            raise ValueError(f"Changed input preserved: {path}")
        path.write_bytes(content)
    save_json(folder / "mapping.json", mapping)
    note = "Author and UniProt numbering agree for selected 1-166"
    if name == "fkbp":
        note = (
            "FKBP author 1-107 maps to current UniProt 2-108; archived reference offset is outdated"
        )
    record = {
        **spec,
        "name": name,
        "author_chain": "A",
        "model": 1,
        "archive_url": url,
        "uniprot_url": fasta_url,
        "sequence": target,
        "mapping_basis": "Exact selected-coordinate/UniProt-region identity",
        "numbering_note": note,
        "cofactors": [{"name": r.resname, "author_residue": r.id[1]} for r in cofactors],
        "mg_waters": [r.id[1] for r in waters],
        "files": {
            p.name: digest(p)
            for p in [archive, *[folder / n for n in files], folder / "mapping.json"]
        },
    }
    save_json(folder / "input.json", record)
    marker.write_text("Archive, canonical sequence, residue mapping and mutation sites verified\n")
    print(f"{name}: {len(target)} residues; mapping saved to {folder / 'mapping.json'}")


def verify_files(folder, files):
    for filename, expected in files.items():
        if digest(Path(folder) / filename) != expected:
            raise ValueError(f"Input changed: {filename}; use a fresh input/run")


def model(root, name):
    root = Path(root).resolve()
    folder = root / "data/01_structure" / name
    if not (folder / "STRUCTURE_DONE").exists():
        raise ValueError("Run structure fetch and mapping first")
    record = json.loads((folder / "input.json").read_text())
    verify_files(folder, record["files"])
    marker = folder / "MODEL_DONE"
    if marker.exists():
        existing = json.loads((folder / "models.json").read_text())
        verify_files(folder, existing["files"])
        if existing["input_sha256"] != digest(folder / "input.json"):
            raise ValueError("Model inputs changed")
        print(f"{name}: verified existing models")
        return
    spec = PROTEINS[name]
    wt = folder / "models/WT/WT.pdb"
    cli(
        "modeller",
        "build",
        "--pdb-id",
        spec["pdb"],
        "--uniprot-id",
        spec["uniprot"],
        "--chain",
        "A",
        "--range",
        spec["start"],
        spec["end"],
        "--uniprot-numbering",
        "--template",
        folder / "template.pdb",
        "--sequence",
        folder / "uniprot.fasta",
        "--outdir",
        wt.parent,
    )
    verify_sequence(wt, "A", record["sequence"], spec["start"])
    inputs = {}
    for variant in spec["variants"]:
        expected = record["sequence"]
        path = wt
        if variant != "WT":
            position = int(variant[1:-1])
            path = folder / "models" / variant / f"{variant}.pdb"
            cli(
                "modeller",
                "mutate",
                "--pdb-in",
                wt,
                "--chain",
                "A",
                "--mut",
                variant,
                "--outdir-mut",
                path.parent,
            )
            idx = position - spec["start"]
            expected = expected[:idx] + variant[-1] + expected[idx + 1 :]
        verify_sequence(path, "A", expected, spec["start"])
        # MODELLER aligns the WT to the original template; mutation preserves that frame.
        protein_lines = [
            line
            for line in path.read_text().splitlines(keepends=True)
            if line.startswith(("ATOM  ", "TER"))
        ]
        heterogens = [
            line
            for line in (folder / "complex.pdb").read_text().splitlines(keepends=True)
            if line.startswith("HETATM")
        ]
        out = folder / f"{variant}.pdb"
        out.write_text("".join([*protein_lines, *heterogens, "END\n"]))
        verify_sequence(out, "A", expected, spec["start"])
        if name == "hras":
            check_complex(out)
        inputs[variant] = {"path": out.name, "sequence": expected, "sha256": digest(out)}
    save_json(
        folder / "models.json",
        {
            "inputs": inputs,
            "input_sha256": digest(folder / "input.json"),
            "modeling": {
                "seed": -49837,
                "python": sys.version,
                "implementation": {
                    n: digest(Path(__file__).resolve().parents[1] / "modeller" / n)
                    for n in ["engine.py", "structure.py", "commands.py"]
                },
            },
            "files": {v["path"]: v["sha256"] for v in inputs.values()},
        },
    )
    marker.write_text(
        "WT and mutant sequences, numbering and coordinates verified; review model scores\n"
    )


def check_complex(path):
    import numpy as np

    structure = PDBParser(QUIET=True).get_structure("complex", str(path))
    residues = list(next(iter(structure)).get_residues())
    if sorted(r.resname for r in residues if r.resname in {"GDP", "MG"}) != ["GDP", "MG"]:
        raise ValueError("GDP/Mg composition changed")
    protein = np.array([a.coord for r in residues if r.id[0] == " " for a in r if a.element != "H"])
    ligand = np.array(
        [a.coord for r in residues if r.resname == "GDP" for a in r if a.element != "H"]
    )
    if np.linalg.norm(protein[:, None] - ligand, axis=2).min() < 1.0:
        raise ValueError("Severe protein/GDP coordinate clash")
    prep = Path(__file__).resolve().parents[1] / "amber/params/GDP-B3.prepi"
    expected = {
        line.split()[1]
        for line in prep.read_text().splitlines()
        if len(line.split()) == 11
        and line.split()[0].isdigit()
        and line.split()[1] != "DUMM"
        and not line.split()[1].startswith("H")
    }
    observed = {a.name for r in residues if r.resname == "GDP" for a in r if a.element != "H"}
    if observed != expected:
        raise ValueError(f"GDP atom names differ from parameter template: {observed ^ expected}")


def simulate(root, name, step, variant="WT"):
    """Run one public CLI stage; one immutable run ID is used per protein/variant."""
    root = Path(root).resolve()
    spec = PROTEINS[name]
    if variant not in spec["variants"]:
        raise ValueError(f"Choose a reviewed variant: {spec['variants']}")
    run_id = os.environ.get("COLABMDA_RUN_ID", "")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", run_id):
        raise ValueError("Set a fresh COLABMDA_RUN_ID, or the original ID to resume")
    os.environ.update(COLABMDA_PLATFORM="CUDA", COLABMDA_REQUIRE_GPU="1", MPLBACKEND="Agg")
    folder = root / "data/01_structure" / name
    if not (folder / "MODEL_DONE").exists():
        raise ValueError("Build and verify the models first")
    inputs = json.loads((folder / "input.json").read_text())
    models = json.loads((folder / "models.json").read_text())
    verify_files(folder, inputs["files"])
    verify_files(folder, models["files"])
    if models["input_sha256"] != digest(folder / "input.json"):
        raise ValueError("Model inputs changed")
    receipt = root / "data/00_setup" / f"{name}.json"
    setup = json.loads(receipt.read_text())
    verify_files(root, setup["files"])
    run = root / "data/02_simulation" / run_id
    system = run / spec["state"] / variant
    system.mkdir(parents=True, exist_ok=True)
    (root / "logs/02_simulation" / run_id).mkdir(parents=True, exist_ok=True)
    provenance = {
        "example": name,
        "variant": variant,
        "state": spec["state"],
        "setup": setup,
        "input_sha256": digest(folder / "input.json"),
        "models_sha256": digest(folder / "models.json"),
        "python": sys.version,
        "tool_version": version("colabmda"),
        "packages": {
            p: version(p) for p in ["openmm", "numpy", "biopython", "MDAnalysis", "mdtraj"]
        },
        "implementation": {
            str(p.relative_to(Path(__file__).resolve().parents[1])): digest(p)
            for p in [
                Path(__file__),
                Path(__file__).resolve().parents[1] / "openmm/md.py",
                Path(__file__).resolve().parents[1] / "cli.py",
            ]
        },
        "protocol": "quick",
        "nvt_ps": 100,
        "npt_ps": 100,
        "step_fs": 2,
        "frame_ps": 10,
        "checkpoint_ps": 100,
    }
    if (run / "provenance.json").exists() and json.loads(
        (run / "provenance.json").read_text()
    ) != provenance:
        raise ValueError("Run identity/settings changed; choose a new run ID")
    save_json(run / "provenance.json", provenance)
    args = ["--root", str(root), "--state", spec["state"], "--variant", variant]
    if step == "prepare":
        if (system / "01_prep/PREP_DONE").exists():
            recorded = json.loads((system / "01_prep/topology_check.json").read_text())
            verify_files(system / "01_prep", recorded["files"])
            print("Verified existing preparation; preserving coordinates/topology")
            return
        cli(
            "prep",
            *args,
            "--protein",
            name,
            "--receptor",
            folder / models["inputs"][variant]["path"],
        )
        for filename in [
            "model.pdb",
            "solvated.pdb",
            "system.prmtop",
            "system.inpcrd",
            "prep.json",
        ]:
            if not (system / "01_prep" / filename).stat().st_size:
                raise ValueError("Incomplete preparation")
        verify_sequence(
            system / "01_prep/model.pdb", "A", models["inputs"][variant]["sequence"], spec["start"]
        )
        if name == "hras":
            check_complex(system / "01_prep/model.pdb")
        check_topology(system / "01_prep", name)
        (system / "01_prep/PREP_DONE").write_text("Verified preparation outputs\n")
    elif step == "equilibrate":
        if not (system / "01_prep/PREP_DONE").exists():
            raise ValueError("Preparation incomplete")
        cli(
            "sim",
            *args,
            "--stage",
            "eq",
            "--protocol",
            "quick",
            "--equil-time",
            "100",
            "--software-test",
        )
        for filename in ["system.xml", "npt.state.xml", "equilibration_qc.json"]:
            if not (system / "02_eq" / filename).stat().st_size:
                raise ValueError("Incomplete equilibration")
        (system / "02_eq/EQ_DONE").write_text(
            "Quick equilibration outputs verified; review physical QC\n"
        )
    elif step == "production":
        if not (system / "02_eq/EQ_DONE").is_file():
            raise ValueError("Equilibration incomplete")
        ns = float(os.environ.get("COLABMDA_EXAMPLE_NS", "0.1"))
        if ns <= 0 or abs(ns * 10 - round(ns * 10)) > 1e-8:
            raise ValueError("Use a positive production duration divisible by 0.1 ns")
        requests = run / "production.json"
        history = json.loads(requests.read_text()) if requests.exists() else []
        if history and ns < history[-1]["requested_ns"]:
            raise ValueError("Cannot shorten a completed/requested run; choose a new run ID")
        history.append({"requested_ns": ns})
        save_json(requests, history)
        (run / "RUN_DONE").unlink(missing_ok=True)
        cli(
            "sim",
            *args,
            "--software-test",
            "--stage",
            "md",
            "--replica",
            "r1",
            "--total-ns",
            ns,
            "--checkpoint-ps",
            "100",
            "--traj-interval",
            "10",
        )
    elif step == "analysis":
        cli("merge", *args, "--replicas", "r1")
        cli("analyze", *args)
    elif step == "check":
        check_run(run, system, name, variant, models["inputs"][variant]["sequence"])
    else:
        raise ValueError(f"Unknown step: {step}")


def check_run(run, system, name, variant, sequence):
    import mdtraj as md
    import numpy as np
    import pandas as pd

    ns = json.loads((run / "production.json").read_text())[-1]["requested_ns"]
    prod = system / "03_md/r1"
    trajectory = md.load(str(prod / "prod_keep_ligand.dcd"), top=str(prod / "stripped.prmtop"))
    expected_times = np.arange(10, ns * 1000 + 0.001, 10)
    if trajectory.n_frames != len(expected_times) or not np.isfinite(trajectory.xyz).all():
        raise ValueError("Invalid trajectory frames/coordinates")
    if len(list(prod.glob("prod_*to*ps.dcd"))) != round(ns * 10):
        raise ValueError("Unexpected chunk count")
    log = pd.read_csv(prod / "prod_full.log")
    if len(log) != len(expected_times) or not np.allclose(log["Time (ps)"], expected_times):
        raise ValueError("Production times have gaps/overlaps")
    if not np.isfinite(log.select_dtypes(include="number").to_numpy()).all():
        raise ValueError("Nonfinite production log")
    analysis = system / "04_analysis"
    for metric in ["rmsd", "rg", "rmsf"]:
        table = pd.read_csv(analysis / f"{metric}.csv")
        if table.empty or not np.isfinite(table.to_numpy()).all():
            raise ValueError("Invalid numerical analysis")
    if pd.read_csv(analysis / "rmsf.csv")["Residue"].tolist() != list(
        range(PROTEINS[name]["start"], PROTEINS[name]["end"] + 1)
    ):
        raise ValueError("Biological residue numbering lost")
    ca = trajectory.topology.select("name CA and protein")
    if len(ca) != len(sequence):
        raise ValueError("Protein C-alpha count changed")
    if name == "hras":
        if (
            len(trajectory.topology.select("resname GDP")) != 40
            or len(trajectory.topology.select("resname MG")) != 1
        ):
            raise ValueError("GDP/Mg lost during preparation or merging")
    for filename in ["rmsd_vs_time.png", "rg_vs_time.png", "rmsf_per_residue.png"]:
        if not (analysis / "figures" / filename).stat().st_size:
            raise ValueError("Missing analysis figure")
    if not (system / "02_eq/figures/equilibration_qc.png").stat().st_size:
        raise ValueError("Missing equilibration figure")
    qc = json.loads((system / "02_eq/equilibration_qc.json").read_text())
    save_json(
        run / "summary.json",
        {
            "example": name,
            "variant": variant,
            "production_ns": ns,
            "frames": trajectory.n_frames,
            "residues": len(ca),
            "software_pass": True,
            "equilibration_qc_pass": qc["passed"],
            "analysis": json.loads((analysis / "analysis.json").read_text()),
        },
    )
    (run / "RUN_DONE").write_text("Software checks passed; physical equilibration QC is separate\n")
    print(f"RUN_DONE: {name} {variant}, {ns} ns; physical QC: {qc['passed']}")


def check_topology(folder, name=None, *, spec=None):
    """Check the actual Amber charges/composition without running dynamics."""
    from openmm import NonbondedForce, unit
    from openmm.app import AmberPrmtopFile

    folder = Path(folder)
    topology = AmberPrmtopFile(str(folder / "system.prmtop"))
    system = topology.createSystem()
    force = next(f for f in system.getForces() if isinstance(f, NonbondedForce))
    charges = [
        force.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)
        for i in range(system.getNumParticles())
    ]
    total = sum(charges)
    if abs(total) > 0.001:
        raise ValueError(f"Solvated system is not neutral: {total}")
    components = [
        {
            "residue": r.name,
            "atoms": len(list(r.atoms())),
            "charge": sum(charges[a.index] for a in r.atoms()),
        }
        for r in topology.topology.residues()
        if r.name in {"ATP", "ADP", "GTP", "GDP", "MG", "MN"}
    ]
    expected = {"GDP": (40, -3), "MG": (1, 2)} if name == "hras" else {}
    composition = {key: 1 for key in expected}
    if spec is not None:
        if spec.ligand:
            counts = {"ATP": (43, -4), "ADP": (39, -3), "GTP": (44, -4), "GDP": (40, -3)}
            if spec.ligand not in counts:
                raise ValueError("Topology validation requires a supported nucleotide")
            expected[spec.ligand] = counts[spec.ligand]
            composition[spec.ligand] = 1
        for ion, count in spec.ions.items():
            if ion != "MG":
                raise ValueError("This validation route currently supports Mg only")
            expected[ion] = (1, 2)
            composition[ion] = count
    actual = {key: sum(r["residue"] == key for r in components) for key in expected}
    if actual != composition or any(r["residue"] not in expected for r in components):
        raise ValueError("Wrong prepared nucleotide/ion composition")
    for residue in components:
        count, charge = expected[residue["residue"]]
        if residue["atoms"] != count or abs(residue["charge"] - charge) > 0.0001:
            raise ValueError("Wrong nucleotide/ion topology atoms or charge")
    # Defaults are OPC with its complete 12-6 ion set. Loading the 12-6-4 LJ
    # values without adding the separate C4 term is a different, incomplete model.
    expected_ions = {
        "Na+": (1.467, 0.02960343),
        "Cl-": (2.360, 0.67878870),
        "MG": (1.330, 0.00716930),
    }
    ion_lj = {}
    for residue in topology.topology.residues():
        if residue.name not in expected_ions:
            continue
        atom = next(residue.atoms())
        _, sigma, epsilon = force.getParticleParameters(atom.index)
        rmin_half = sigma.value_in_unit(unit.angstrom) * 2 ** (1 / 6) / 2
        eps = epsilon.value_in_unit(unit.kilocalorie_per_mole)
        expected_rmin, expected_epsilon = expected_ions[residue.name]
        if abs(rmin_half - expected_rmin) > 1e-5 or abs(eps - expected_epsilon) > 1e-7:
            raise ValueError(f"Incomplete or mismatched OPC ion model for {residue.name}")
        ion_lj[residue.name] = {"rmin_half_angstrom": rmin_half, "epsilon_kcal_mol": eps}
    save_json(
        folder / "topology_check.json",
        {
            "atoms": system.getNumParticles(),
            "total_charge": total,
            "components": components,
            "ion_model": "OPC Li-Merz 12-6",
            "ion_lj": ion_lj,
            "files": {
                n: digest(folder / n)
                for n in [
                    "model.pdb",
                    "system.prmtop",
                    "system.inpcrd",
                    "solvated.pdb",
                    "prep.json",
                ]
            },
        },
    )
