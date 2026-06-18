import argparse
import json
import os

from openmm import MonteCarloBarostat, XmlSerializer, unit
from openmm.app import (
    PME,
    AmberInpcrdFile,
    AmberPrmtopFile,
    ForceField,
    HBonds,
    Modeller,
    PDBFile,
)
from pdbfixer import PDBFixer

from colabmda.openmm.modular.utils import configure_system, make_sim, save_stage_state

FORCE_FIELD_CHOICES = {
    ("ff19SB", "opc"): ("amber19-all.xml", "amber19/opc.xml", "AMBER ff19SB", "OPC", "tip4pew"),
    (
        "ff19SB",
        "tip3p",
    ): ("amber19-all.xml", "amber19/tip3p.xml", "AMBER ff19SB", "TIP3P", "tip3p"),
    (
        "ff19SB",
        "tip3pfb",
    ): ("amber19-all.xml", "amber19/tip3pfb.xml", "AMBER ff19SB", "TIP3P-FB", "tip3p"),
    (
        "ff19SB",
        "tip4pew",
    ): ("amber19-all.xml", "amber19/tip4pew.xml", "AMBER ff19SB", "TIP4P-Ew", "tip4pew"),
    ("ff14SB", "opc"): ("amber14-all.xml", "amber14/opc.xml", "AMBER ff14SB", "OPC", "tip4pew"),
    (
        "ff14SB",
        "tip3p",
    ): ("amber14-all.xml", "amber14/tip3p.xml", "AMBER ff14SB", "TIP3P", "tip3p"),
    (
        "ff14SB",
        "tip3pfb",
    ): ("amber14-all.xml", "amber14/tip3pfb.xml", "AMBER ff14SB", "TIP3P-FB", "tip3p"),
    (
        "ff14SB",
        "tip4pew",
    ): ("amber14-all.xml", "amber14/tip4pew.xml", "AMBER ff14SB", "TIP4P-Ew", "tip4pew"),
}


def _select_forcefield(protein_ff, water_model):
    key = (protein_ff, water_model.lower())
    if key not in FORCE_FIELD_CHOICES:
        valid = ", ".join(f"{p}/{w}" for p, w in sorted(FORCE_FIELD_CHOICES))
        raise ValueError(f"Unsupported force-field/water combination {key!r}. Valid: {valid}")
    protein_xml, water_xml, protein_label, water_label, solvent_geometry = FORCE_FIELD_CHOICES[key]
    return {
        "protein_xml": protein_xml,
        "water_xml": water_xml,
        "protein_label": protein_label,
        "water_label": water_label,
        "solvent_geometry": solvent_geometry,
    }


def _load_forcefield(*templates):
    try:
        return ForceField(*templates)
    except ValueError as exc:
        if "Could not locate file" in str(exc):
            raise ValueError(
                f"OpenMM could not load force-field XML files {templates}. "
                "Install a newer OpenMM for Amber19 support, or choose an older installed "
                "combination such as --protein-ff ff14SB --water-model opc."
            ) from exc
        raise


def _resolve_padding(positions, requested):
    if requested != "auto":
        padding_nm = float(requested)
    else:
        xyz = [p.value_in_unit(unit.nanometer) for p in positions]
        spans = [max(p[i] for p in xyz) - min(p[i] for p in xyz) for i in range(3)]
        longest = max(spans)
        # Padding is primarily constrained by the 1 nm cutoff. Extra margin is
        # conservative for unusually extended systems, which are more likely to move.
        padding_nm = 1.0 if longest <= 8.0 else 1.2 if longest <= 12.0 else 1.4
        estimated = [span + 2 * padding_nm for span in spans]
        print(
            "  • Auto box: solute spans "
            + " × ".join(f"{value:.2f}" for value in spans)
            + " nm; selected padding "
            + f"{padding_nm:.2f} nm; estimated box "
            + " × ".join(f"{value:.2f}" for value in estimated)
            + " nm"
        )
    if padding_nm < 1.0:
        raise ValueError(
            "Solvent padding must be at least the 1.0 nm nonbonded cutoff; "
            "use 1.0 nm for a compact box or 1.4 nm for conservative padding."
        )
    return padding_nm


def _as_abs_or_none(path):
    return os.path.abspath(path) if path else None


def _assert_compatible_cached_em(workdir, requested):
    manifest = os.path.join(workdir, "parameterization.json")
    if not os.path.exists(manifest):
        raise RuntimeError(
            "Existing EM state has no parameterization.json; use a fresh workdir for this protocol"
        )
    with open(manifest) as handle:
        recorded = json.load(handle)
    for key, value in requested.items():
        if recorded.get(key) != value:
            raise RuntimeError(
                f"Existing EM state was prepared with {key}={recorded.get(key)!r}, "
                f"but this run requested {value!r}. Use a fresh workdir."
            )


def run_em(
    workdir,
    pdbid,
    ligand=None,
    keep_mg=False,
    amber_prmtop=None,
    amber_inpcrd=None,
    padding_nm="auto",
    protein_ff="ff19SB",
    water_model="opc",
    small_molecule_ff="gaff-2.11",
):
    workdir = os.path.abspath(workdir)
    os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)

    cleaned_pdb = os.path.join(workdir, f"{pdbid}_cleaned.pdb")
    if not os.path.exists(cleaned_pdb):
        print(f"Error: cleaned PDB not found: {cleaned_pdb}")
        return False

    xml_system = os.path.join(workdir, "system.xml")
    pdb_saved = os.path.join(workdir, "solvated.pdb")
    chk_file = os.path.join(workdir, "em.chk")
    state_file = os.path.join(workdir, "em.state.xml")

    print(f"[DEBUG] Checking for EM files in: {workdir}")
    ff_choice = _select_forcefield(protein_ff, water_model)
    requested_em = {
        "requested_ligand": _as_abs_or_none(ligand),
        "requested_amber_prmtop": _as_abs_or_none(amber_prmtop),
        "requested_amber_inpcrd": _as_abs_or_none(amber_inpcrd),
        "requested_padding_nm": str(padding_nm),
        "protein": ff_choice["protein_label"],
        "water": ff_choice["water_label"],
        "protein_forcefield_xml": ff_choice["protein_xml"],
        "water_model_xml": ff_choice["water_xml"],
        "small_molecule_forcefield": small_molecule_ff if ligand else None,
    }
    if os.path.exists(state_file) and os.path.exists(xml_system):
        _assert_compatible_cached_em(workdir, requested_em)
        print(f"✔ EM already complete ({state_file} found). Skipping...")
        return True

    print("▶ Starting System Preparation & Minimization …")
    parameterization = {}
    if bool(amber_prmtop) != bool(amber_inpcrd):
        raise ValueError("--amber-prmtop and --amber-inpcrd must be supplied together")

    if amber_prmtop:
        print("  • Loading pre-parameterized AMBER topology/coordinates…")
        prmtop = AmberPrmtopFile(os.path.abspath(amber_prmtop))
        inpcrd = AmberInpcrdFile(os.path.abspath(amber_inpcrd))
        modeller = Modeller(prmtop.topology, inpcrd.positions)
        if inpcrd.boxVectors is not None:
            modeller.topology.setPeriodicBoxVectors(inpcrd.boxVectors)
        system = prmtop.createSystem(
            nonbondedMethod=PME,
            nonbondedCutoff=1.0 * unit.nanometer,
            constraints=HBonds,
            rigidWater=True,
        )
        parameterization = {
            "source": "amber-prmtop",
            "prmtop": os.path.abspath(amber_prmtop),
            "inpcrd": os.path.abspath(amber_inpcrd),
            "note": "Preserves VarMDyn/LEaP parameters, including specialized ligand parameters.",
        }
    elif ligand:
        # Dynamic imports to avoid dependency issues if openff/openmmforcefields is not installed
        from openff.toolkit.topology import Molecule
        from openmmforcefields.generators import SystemGenerator

        print(f"  • Loading ligand from: {ligand}")
        ligand_mol = Molecule.from_file(ligand)
        if isinstance(ligand_mol, list):
            if len(ligand_mol) != 1:
                raise ValueError("Ligand file must contain exactly one molecule")
            ligand_mol = ligand_mol[0]
        if len(ligand_mol.conformers) == 0:
            raise ValueError("Ligand file must contain a 3D conformer aligned to the receptor")
        charge_source = "input"
        if ligand_mol.partial_charges is None or all(
            abs(float(value)) < 1e-8 for value in ligand_mol.partial_charges.magnitude
        ):
            print("  • Assigning ligand AM1-BCC partial charges…")
            ligand_mol.assign_partial_charges("am1bcc")
            charge_source = "OpenFF AM1-BCC"
        print(
            f"  • Ligand formal charge: {ligand_mol.total_charge}; partial charges: {charge_source}"
        )

        # Filter cleaned_pdb to keep only protein residues, water, and optionally MG
        # Standard residues, protonation states, capping groups, and ions
        protein_residues = {
            "ALA",
            "ARG",
            "ASN",
            "ASP",
            "CYS",
            "GLN",
            "GLU",
            "GLY",
            "HIS",
            "ILE",
            "LEU",
            "LYS",
            "MET",
            "PHE",
            "PRO",
            "SER",
            "THR",
            "TRP",
            "TYR",
            "VAL",
            # Protonation states
            "HID",
            "HIE",
            "HIP",
            "CYX",
            "ASH",
            "GLH",
            "LYN",
            # Capping groups
            "ACE",
            "NME",
            "NH2",
            # Selenomethionine
            "MSE",
        }
        allowed_resnames = set(protein_residues)
        allowed_resnames.update(["HOH", "WAT", "NA", "CL", "K", "ZN", "CA"])
        if keep_mg:
            allowed_resnames.add("MG")
            print("  • Preserving Magnesium (MG) ions from starting structure")

        filtered_lines = []
        with open(cleaned_pdb) as f_in:
            for line in f_in:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    resname = line[17:20].strip()
                    if resname in allowed_resnames:
                        filtered_lines.append(line)
                elif line.startswith("TER") or line.startswith("CONECT") or line.startswith("END"):
                    filtered_lines.append(line)

        # If keep_mg is true, extract MG atoms from raw PDB if they are not in cleaned_pdb
        if keep_mg:
            has_mg = any(
                line[17:20].strip() == "MG"
                for line in filtered_lines
                if line.startswith("ATOM") or line.startswith("HETATM")
            )
            if not has_mg:
                # Find raw PDB in workdir, or fallback to structures directory
                raw_pdb = os.path.join(workdir, f"{pdbid}.pdb")
                candidate_paths = [raw_pdb]

                # Extract base PDB ID (e.g. 4ldj from 4ldj_gdp)
                base_pdb = pdbid.split("_")[0]

                # Locate the project root from either data/sim (current), sim, or
                # simulations (legacy). A fixed parent offset is incorrect for
                # data/sim because that layout has one additional path segment.
                from pathlib import Path as pathlib_Path

                workdir_path = pathlib_Path(workdir)
                root_dir = None
                for ancestor in (workdir_path, *workdir_path.parents):
                    if ancestor.name not in {"sim", "simulations"}:
                        continue
                    if ancestor.name == "sim" and ancestor.parent.name == "data":
                        root_dir = ancestor.parent.parent
                    else:
                        root_dir = ancestor.parent
                    break

                if root_dir is not None:
                    # Priority: data/str > str > structures (legacy)
                    candidate_paths.append(
                        str(root_dir / "data" / "str" / base_pdb / f"{base_pdb}_orig.pdb")
                    )
                    candidate_paths.append(
                        str(root_dir / "data" / "str" / base_pdb / f"{base_pdb}.pdb")
                    )
                    candidate_paths.append(
                        str(root_dir / "str" / base_pdb / f"{base_pdb}_orig.pdb")
                    )
                    candidate_paths.append(
                        str(root_dir / "structures" / base_pdb / f"{base_pdb}_orig.pdb")
                    )
                    candidate_paths.append(str(root_dir / "str" / base_pdb / f"{base_pdb}.pdb"))
                    candidate_paths.append(
                        str(root_dir / "structures" / base_pdb / f"{base_pdb}.pdb")
                    )
                    candidate_paths.append(str(root_dir / "data" / "str" / f"{base_pdb}_orig.pdb"))
                    candidate_paths.append(str(root_dir / "str" / f"{base_pdb}_orig.pdb"))
                    candidate_paths.append(str(root_dir / "structures" / f"{base_pdb}_orig.pdb"))

                # Also check local CWD-relative str and structures folders
                candidate_paths.append(str(pathlib_Path("str") / base_pdb / f"{base_pdb}_orig.pdb"))
                candidate_paths.append(
                    str(pathlib_Path("structures") / base_pdb / f"{base_pdb}_orig.pdb")
                )

                found_raw = None
                for path in candidate_paths:
                    if os.path.exists(path):
                        found_raw = path
                        break

                if found_raw:
                    mg_lines = []
                    with open(found_raw) as f_raw:
                        for line in f_raw:
                            if line.startswith("ATOM") or line.startswith("HETATM"):
                                resname = line[17:20].strip()
                                if resname == "MG":
                                    mg_lines.append(line)
                    if mg_lines:
                        print(
                            f"  • Extracted {len(mg_lines)} Magnesium (MG) ions from raw PDB: {found_raw}"
                        )
                        # Insert before END or CONECT lines
                        insert_idx = len(filtered_lines)
                        for idx, line in enumerate(filtered_lines):
                            if line.startswith("END") or line.startswith("CONECT"):
                                insert_idx = idx
                                break
                        filtered_lines[insert_idx:insert_idx] = mg_lines

        temp_cleaned_pdb = cleaned_pdb + ".tmp_filter"
        with open(temp_cleaned_pdb, "w") as f_out:
            f_out.writelines(filtered_lines)

        # Run PDBFixer on the filtered structure without removing kept heterogens
        fixer = PDBFixer(filename=temp_cleaned_pdb)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(pH=7.0)

        # Clean up temp file
        if os.path.exists(temp_cleaned_pdb):
            os.remove(temp_cleaned_pdb)

        # Create combined Modeller object
        modeller = Modeller(fixer.topology, fixer.positions)
        lig_openmm_top = ligand_mol.to_topology().to_openmm()
        lig_openmm_pos = ligand_mol.conformers[0].to_openmm()
        modeller.add(lig_openmm_top, lig_openmm_pos)
        resolved_padding_nm = _resolve_padding(modeller.positions, padding_nm)

        # Initialize SystemGenerator with the selected small-molecule force field
        templates = [ff_choice["protein_xml"], ff_choice["water_xml"]]
        system_generator = SystemGenerator(
            forcefields=templates,
            small_molecule_forcefield=small_molecule_ff,
            molecules=[ligand_mol],
            periodic_forcefield_kwargs={
                "constraints": HBonds,
                "rigidWater": True,
                "nonbondedMethod": PME,
                "nonbondedCutoff": 1.0 * unit.nanometer,
            },
        )

        print(
            f"  • Solvating in {ff_choice['water_label']} with {resolved_padding_nm:.2f} nm padding and neutralizing ions…"
        )
        modeller.addSolvent(
            system_generator.forcefield,
            # OpenMM's Modeller provides a generic pre-equilibrated four-site
            # geometry via tip4pew; the ForceField assigns OPC parameters and
            # minimization relaxes the generated box.
            model=ff_choice["solvent_geometry"],
            padding=resolved_padding_nm * unit.nanometer,
            neutralize=True,
        )

        system = system_generator.create_system(modeller.topology)
        parameterization = {
            "source": "openmmforcefields",
            "protein": ff_choice["protein_label"],
            "water": ff_choice["water_label"],
            "ligand": small_molecule_ff,
            "ligand_partial_charges": charge_source,
            "formal_charge": str(ligand_mol.total_charge),
            "warning": "Generic small-molecule parameters are not equivalent to curated ATP-Mg parameters.",
        }
    else:
        # Legacy protein-water setup pathway
        fixer = PDBFixer(filename=cleaned_pdb)
        fixer.removeHeterogens(keepWater=True)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(pH=7.0)

        modeller = Modeller(fixer.topology, fixer.positions)
        resolved_padding_nm = _resolve_padding(modeller.positions, padding_nm)
        ff = _load_forcefield(ff_choice["protein_xml"], ff_choice["water_xml"])

        print(
            f"  • Solvating in {ff_choice['water_label']} with {resolved_padding_nm:.2f} nm padding and neutralizing ions…"
        )
        modeller.addSolvent(
            ff,
            model=ff_choice["solvent_geometry"],
            padding=resolved_padding_nm * unit.nanometer,
            neutralize=True,
        )

        system = ff.createSystem(
            modeller.topology,
            nonbondedMethod=PME,
            nonbondedCutoff=1.0 * unit.nanometer,
            constraints=HBonds,
            rigidWater=True,
        )
        parameterization = {
            "source": "OpenMM XML",
            "protein": ff_choice["protein_label"],
            "water": ff_choice["water_label"],
        }

    system.addForce(MonteCarloBarostat(1 * unit.bar, 300 * unit.kelvin, 25))
    configure_system(system)

    dt = 2.0 * unit.femtoseconds
    sim = make_sim(modeller.topology, system, dt)
    sim.context.setPositions(modeller.positions)

    print("  • Minimizing…")
    sim.minimizeEnergy()

    # Save files
    with open(xml_system, "w") as f:
        f.write(XmlSerializer.serializeSystem(sim.context.getSystem()))
    with open(pdb_saved, "w") as f:
        PDBFile.writeFile(sim.topology, sim.context.getState(getPositions=True).getPositions(), f)
    sim.saveCheckpoint(chk_file)
    save_stage_state(sim, state_file)
    parameterization.update(
        {
            **requested_em,
            "temperature_kelvin": 300,
            "pressure_bar": 1,
            "cutoff_nm": 1.0,
            "switching_distance_nm": 0.8,
            "solvent_padding_nm": None if amber_prmtop else resolved_padding_nm,
            "hydrogen_bond_constraints": True,
        }
    )
    with open("parameterization.json", "w") as handle:
        json.dump(parameterization, handle, indent=2)
        handle.write("\n")
    print("✔ Minimization complete (em.state.xml saved; em.chk retained for same-context recovery)")
    return True


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("workdir")
    p.add_argument("--pdbid", required=True)
    p.add_argument("--ligand", default=None)
    p.add_argument("--keep-mg", action="store_true")
    p.add_argument("--amber-prmtop")
    p.add_argument("--amber-inpcrd")
    p.add_argument("--padding-nm", default="auto")
    p.add_argument("--protein-ff", choices=["ff19SB", "ff14SB"], default="ff19SB")
    p.add_argument("--water-model", choices=["opc", "tip3p", "tip3pfb", "tip4pew"], default="opc")
    p.add_argument("--small-molecule-ff", default="gaff-2.11")
    args = p.parse_args()
    run_em(
        args.workdir,
        args.pdbid,
        ligand=args.ligand,
        keep_mg=args.keep_mg,
        amber_prmtop=args.amber_prmtop,
        amber_inpcrd=args.amber_inpcrd,
        padding_nm=args.padding_nm,
        protein_ff=args.protein_ff,
        water_model=args.water_model,
        small_molecule_ff=args.small_molecule_ff,
    )
