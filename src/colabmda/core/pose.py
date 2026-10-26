"""Public homologous ATP/Mg transfer with fixed-backbone pocket relaxation."""

import copy
import json
import re
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB import PDBIO, PDBParser
from Bio.SeqUtils import seq1

from colabmda.core.project import ProjectPaths, resolve_project_root
from colabmda.core.taxonomy import SystemSpec
from colabmda.modeller.structure import fetch_archive, verify_sequence
from colabmda.modeller.substrate import kabsch, terminate_components


def fit_donor(target, donor, observed, core):
    """Align homologous sequence positions, iteratively rejecting CA outliers."""
    fixed = [
        r
        for r in target
        if r.id[0] == " " and core[0] <= r.id[1] <= core[1] and r.id[1] in observed
    ]
    moving = [r for r in donor if r.id[0] == " " and "CA" in r]
    if len(fixed) < 30 or len(moving) < 30:
        raise ValueError("Insufficient experimentally resolved core coordinates")
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(
        "".join(seq1(r.resname) for r in fixed), "".join(seq1(r.resname) for r in moving)
    )[0]
    pairs = [
        (fixed[i], moving[j])
        for a, b in zip(*alignment.aligned)
        for i, j in zip(range(*a), range(*b))
    ]
    x = np.array([b["CA"].coord for a, b in pairs])
    y = np.array([a["CA"].coord for a, b in pairs])
    keep = np.ones(len(pairs), dtype=bool)
    initial = len(pairs)
    for _ in range(6):
        if keep.sum() < 30:
            raise ValueError("Too few homologous core pairs after outlier rejection")
        rotation, translation = kabsch(x[keep], y[keep])
        distances = np.linalg.norm(x @ rotation + translation - y, axis=1)
        revised = keep & (distances <= 2.0)
        if np.array_equal(keep, revised):
            break
        keep = revised
    if keep.sum() < 30 or keep.sum() < 0.5 * initial:
        raise ValueError("Homologous donor fit has insufficient retained core coverage")
    rotation, translation = kabsch(x[keep], y[keep])
    rmsd = float(
        np.sqrt(np.mean(np.sum((x[keep] @ rotation + translation - y[keep]) ** 2, axis=1)))
    )
    if rmsd > 2.0:
        raise ValueError("Homologous core fit exceeds 2 Angstrom")
    record = dict(
        method="BLOSUM62 global alignment; CA Kabsch; six rounds of 2 A outlier rejection",
        initial_pairs=initial,
        retained_pairs=int(keep.sum()),
        core_rmsd_angstrom=rmsd,
        fitted_residues=[
            dict(target=a.id[1], target_insertion=a.id[2], donor=b.id[1], donor_insertion=b.id[2])
            for (a, b), use in zip(pairs, keep)
            if use
        ],
    )
    return rotation, translation, record


def minimum_distance(chain):
    protein = np.array([a.coord for r in chain if r.id[0] == " " for a in r if a.element != "H"])
    ligand = np.array([a.coord for r in chain if r.resname == "ATP" for a in r if a.element != "H"])
    return float(np.linalg.norm(protein[:, None] - ligand, axis=2).min())


def relax_pocket(model, folder, chain_id):
    """Minimize pocket sidechains only; donor pose and protein backbone stay fixed.

    This is local model preparation in vacuum, not equilibration or docking.
    No dynamics, solvent, neutralizing ions or implicit-solvent assumptions enter
    this geometry step. The solvated physical model is built separately by prep.
    """
    from openmm import Context, LocalEnergyMinimizer, Platform, VerletIntegrator, unit
    from openmm.app import AmberInpcrdFile, AmberPrmtopFile, NoCutoff

    from colabmda.amber.leap import generate_tleap_script, run_tleap

    folder.mkdir(parents=True, exist_ok=True)
    source = folder / "unrelaxed.pdb"
    writer = PDBIO()
    writer.set_structure(model)
    writer.save(str(source))
    terminate_components(source)
    top, crd = folder / "pocket.prmtop", folder / "pocket.inpcrd"
    script = generate_tleap_script(SystemSpec.from_tag("atpmg"), source, top, crd, solvate=False)
    rc, log = run_tleap(script, folder, keep_script=True)
    (folder / "leap.log").write_text(log)
    if rc or not top.is_file() or not crd.is_file():
        raise RuntimeError("Pocket topology failed; see relaxation/leap.log")
    topology = AmberPrmtopFile(str(top))
    positions = np.asarray(AmberInpcrdFile(str(crd)).positions.value_in_unit(unit.nanometer))
    system = topology.createSystem(nonbondedMethod=NoCutoff, constraints=None, removeCMMotion=False)
    atoms = list(topology.topology.atoms())
    protein = [r for r in model[0][chain_id] if r.id[0] == " "]
    residues = list(topology.topology.residues())
    if len(residues) != len(protein) + 2:
        raise ValueError("Unexpected pocket topology residue count")
    ligand_ids = [a.index for a in atoms if a.residue.name == "ATP"]
    mg_ids = [a.index for a in atoms if a.residue.name == "MG"]
    if len(ligand_ids) != 43 or len(mg_ids) != 1:
        raise ValueError("Pocket topology must contain ATP43 and one Mg")
    backbone = {"N", "CA", "C", "O", "OXT", "H", "H1", "H2", "H3", "HA", "HA2", "HA3"}
    selected = set()
    for i, residue in enumerate(residues[: len(protein)]):
        side = [a for a in residue.atoms() if a.name not in backbone]
        heavy = [a.index for a in side if a.element.symbol != "H"]
        if (
            heavy
            and np.linalg.norm(positions[heavy, None] - positions[ligand_ids], axis=2).min() <= 0.5
        ):
            selected.update(a.index for a in side)
    if not selected:
        raise ValueError("No pocket sidechains selected")
    for atom in atoms:
        if atom.index not in selected:
            system.setParticleMass(atom.index, 0)
    integrator = VerletIntegrator(0.001)
    context = Context(system, integrator, Platform.getPlatformByName("CPU"), {"Threads": "2"})
    context.setPositions(positions)
    initial = (
        context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    )
    LocalEnergyMinimizer.minimize(context, 10.0, 1000)
    state = context.getState(getPositions=True, getEnergy=True)
    result = np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
    final = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    frozen = sorted(set(range(len(atoms))) - selected)
    if not np.isfinite(result).all() or not np.isfinite([initial, final]).all() or final > initial:
        raise ValueError("Pocket relaxation failed finite-energy/descent checks")
    if not np.allclose(result[frozen], positions[frozen], atol=1e-10, rtol=0):
        raise ValueError("Fixed backbone or cofactor coordinates changed")
    displacement = float(np.linalg.norm(result - positions, axis=1).max() * 10)
    if displacement > 5.0:
        raise ValueError("Pocket relaxation moved an atom more than 5 Angstrom")
    for original, amber in zip(protein, residues):
        lookup = {a.name: a.index for a in amber.atoms()}
        for atom in original:
            if atom.name not in lookup:
                raise ValueError(f"Missing protein atom in pocket topology: {atom.name}")
            atom.coord = result[lookup[atom.name]] * 10
    del context, integrator
    return dict(
        method="Amber ff19SB/ATP-B3/OPC Li-Merz12-6; vacuum local minimization; no dynamics",
        selected_atoms=len(selected),
        pocket_radius_angstrom=5.0,
        backbone_and_cofactors_fixed=True,
        max_iterations=1000,
        initial_energy_kj_mol=float(initial),
        final_energy_kj_mol=float(final),
        maximum_atom_displacement_angstrom=displacement,
    )


def assemble_public(args):
    """Build a new public-input model; never consume activation-study files."""
    from colabmda.core.workflow import check, fetch, sha, write_json

    if not re.fullmatch(r"[0-9A-Za-z]{4}", args.template):
        raise ValueError("Public donor must be a four-character archive ID")
    if not re.fullmatch(r"WT|[A-Z][1-9][0-9]*[A-Z]", args.variant):
        raise ValueError("Invalid variant")
    if args.state != "atpmg":
        raise ValueError("Public ATP/Mg transfer requires --state atpmg")
    if any(
        getattr(args, name, None)
        for name in ["receptor", "reference", "ligand", "ions", "mapping", "out"]
    ):
        raise ValueError(
            "Public mode uses state model paths; do not mix explicit reviewed-pose inputs"
        )
    apo = ProjectPaths(resolve_project_root(args.root), state="apo").structure_dir
    inputs = json.loads((apo / "mapping/input.json").read_text())
    target_tag = (
        f"{inputs['uniprot'].upper()}_{inputs['pdb'].upper()}"
        if "uniprot" in inputs and "pdb" in inputs
        else None
    )
    paths = ProjectPaths(resolve_project_root(args.root), state="atpmg", target=target_tag)
    parent = paths.structure_dir / "models" / args.variant
    out = parent / "complex.pdb"
    marker = parent / "ASSEMBLY_DONE"
    source = apo / "models" / args.variant / f"{args.variant}.pdb"
    check(SimpleNamespace(root=str(paths.root), state="apo", variant=args.variant, model=None))
    request = dict(
        model_sha256=sha(source),
        mapping_sha256=sha(apo / "mapping/input.json"),
        donor=args.template.upper(),
        chain=args.chain,
        core=args.core,
        max_rmsd=args.max_rmsd,
        relax_sidechains=args.relax_sidechains,
        implementation_sha256=sha(Path(__file__)),
        leap_sha256=sha(Path(__file__).parents[1] / "amber/leap.py"),
    )
    if out.exists():
        report = json.loads(out.with_suffix(".json").read_text())
        donor_rel = report["donor"].get("archive_path")
        if donor_rel and (paths.structure_dir / donor_rel).is_file():
            archive = paths.structure_dir / donor_rel
        elif (
            paths.structure_dir
            / "donors"
            / report["donor"]["pdb"]
            / report["donor"]["archive_filename"]
        ).is_file():
            archive = (
                paths.structure_dir
                / "donors"
                / report["donor"]["pdb"]
                / report["donor"]["archive_filename"]
            )
        else:
            archive = paths.structure_dir / "raw" / report["donor"]["archive_filename"]
        if (
            not marker.is_file()
            or report.get("request") != request
            or sha(out) != report["output_sha256"]
            or sha(archive) != report["donor"]["archive_sha256"]
        ):
            raise ValueError("Existing assembly inputs/outputs differ; preserved for inspection")
        print("Reused verified public ATP/Mg assembly")
        return 0
    marker.unlink(missing_ok=True)
    fetch(
        SimpleNamespace(
            root=str(paths.root),
            state="atpmg",
            pdb_id=inputs["pdb"],
            uniprot_id=inputs["uniprot"],
            chain=inputs["chain"],
            range=inputs["range"],
            mutations="" if args.variant == "WT" else args.variant,
        )
    )
    if args.chain != inputs["chain"]:
        raise ValueError("Public model chain differs from archive mapping")
    parent.mkdir(parents=True, exist_ok=True)
    receptor = parent / f"{args.variant}.pdb"
    if receptor.exists() and receptor.read_bytes() != source.read_bytes():
        raise ValueError("Different existing state protein preserved")
    receptor.write_bytes(source.read_bytes())
    check(SimpleNamespace(root=str(paths.root), state="atpmg", variant=args.variant, model=None))
    model = PDBParser(QUIET=True).get_structure("target", source)
    donor_dir = paths.structure_dir / "donors" / args.template.upper()
    donor_dir.mkdir(parents=True, exist_ok=True)
    donor, archive, url = fetch_archive(args.template, donor_dir)
    if args.chain not in donor[0]:
        raise ValueError("Requested donor chain absent")
    chain = model[0][args.chain]
    donor_chain = donor[0][args.chain]
    components = [copy.deepcopy(r) for r in donor_chain if r.resname in {"ATP", "MG"}]
    if [r.resname for r in components].count("ATP") != 1 or [r.resname for r in components].count(
        "MG"
    ) != 1:
        raise ValueError("Donor chain must have exactly one ATP and one Mg")
    atp = next(r for r in components if r.resname == "ATP")
    mg = next(r for r in components if r.resname == "MG")
    oxygens = [a.coord for a in atp if a.element == "O"]
    if len([a for a in atp if a.element != "H"]) != 31 or len(list(mg)) != 1:
        raise ValueError("Donor ATP/Mg composition mismatch")
    mg_distance = float(np.linalg.norm(np.array(oxygens) - next(iter(mg)).coord, axis=1).min())
    if not 1.5 <= mg_distance <= 3.0:
        raise ValueError("Donor Mg lacks an ATP oxygen contact at 1.5-3 A")
    rows = json.loads((apo / "mapping/residues.json").read_text())
    rotation, translation, fit = fit_donor(
        chain, donor_chain, {r["uniprot"] for r in rows if r["resolved"]}, args.core
    )
    if fit["core_rmsd_angstrom"] > args.max_rmsd:
        raise ValueError("Homologous donor fit exceeds --max-rmsd")
    for i, residue in enumerate(components, max(r.id[1] for r in chain) + 1):
        residue.detach_parent()
        residue.id = (f"H_{residue.resname}", i, " ")
        for atom in residue:
            atom.coord = atom.coord @ rotation + translation
        chain.add(residue)
    before = minimum_distance(chain)
    relaxation = (
        relax_pocket(model, parent / "relaxation", args.chain) if args.relax_sidechains else None
    )
    after = minimum_distance(chain)
    if after < 1.5:
        raise ValueError(f"Protein/ATP clash ({after:.3f} A); no accepted complex written")
    verify_sequence(
        receptor,
        args.chain,
        "".join(seq1(r.resname) for r in chain if r.id[0] == " "),
        inputs["range"][0],
    )
    temporary = out.with_suffix(".tmp.pdb")
    writer = PDBIO()
    writer.set_structure(model)
    writer.save(str(temporary))
    terminate_components(temporary)
    record = dict(
        request=request,
        coordinate_source="Public homologous donor transfer; modeled CDKL5 pose, not an experimental complex or legacy-study replica",
        donor=dict(
            pdb=args.template.upper(),
            chain=args.chain,
            archive_url=url,
            archive_filename=archive.name,
            archive_path=str(archive.relative_to(paths.structure_dir)),
            archive_sha256=sha(archive),
            magnesium_oxygen_distance_angstrom=mg_distance,
        ),
        variant=args.variant,
        fit=fit,
        relaxation=relaxation,
        minimum_protein_atp_distance_before_angstrom=before,
        minimum_protein_atp_distance_angstrom=after,
        input_sha256={
            str(source): sha(source),
            str(apo / "mapping/input.json"): sha(apo / "mapping/input.json"),
        },
        output_sha256=sha(temporary),
    )
    write_json(out.with_suffix(".json"), record)
    temporary.replace(out)
    marker.write_text(
        "Public donor fit, fixed backbone/cofactors, finite energy and ATP clash checks passed\n"
    )
    print(
        f"Assembled public {args.template.upper()} donor: {fit['retained_pairs']} core pairs, RMSD {fit['core_rmsd_angstrom']:.3f} A, minimum protein/ATP distance {after:.3f} A"
    )
    return 0
