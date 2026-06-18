import argparse
import os

from openmm import MonteCarloBarostat, XmlSerializer, unit
from openmm.app import PME, ForceField, HBonds, Modeller, PDBFile
from pdbfixer import PDBFixer

from colabmda.openmm.modular.utils import make_sim


def run_em(workdir, pdbid, ligand=None, keep_mg=False):
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

    print(f"[DEBUG] Checking for EM files in: {workdir}")
    if os.path.exists(chk_file) and os.path.exists(xml_system):
        print(f"✔ EM already complete ({chk_file} found). Skipping...")
        return True

    print("▶ Starting System Preparation & Minimization …")

    if ligand:
        # Dynamic imports to avoid dependency issues if openff/openmmforcefields is not installed
        from openff.toolkit.topology import Molecule
        from openmmforcefields.generators import SystemGenerator

        print(f"  • Loading ligand from: {ligand}")
        ligand_mol = Molecule.from_file(ligand)

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

                # Check root/str/base_pdb/base_pdb_orig.pdb or root/structures/...
                # Since workdir is typically <root>/sim/<name>/<replica> or <root>/simulations/... parents[2] is <root>
                from pathlib import Path as pathlib_Path

                workdir_path = pathlib_Path(workdir)
                if len(workdir_path.parents) >= 3:
                    root_dir = workdir_path.parents[2]
                    candidate_paths.append(
                        str(root_dir / "str" / base_pdb / f"{base_pdb}_orig.pdb")
                    )
                    candidate_paths.append(
                        str(root_dir / "structures" / base_pdb / f"{base_pdb}_orig.pdb")
                    )
                    candidate_paths.append(
                        str(root_dir / "str" / base_pdb / f"{base_pdb}.pdb")
                    )
                    candidate_paths.append(
                        str(root_dir / "structures" / base_pdb / f"{base_pdb}.pdb")
                    )
                    candidate_paths.append(str(root_dir / "str" / f"{base_pdb}_orig.pdb"))
                    candidate_paths.append(str(root_dir / "structures" / f"{base_pdb}_orig.pdb"))

                # Also check local CWD-relative str and structures folders
                candidate_paths.append(
                    str(pathlib_Path("str") / base_pdb / f"{base_pdb}_orig.pdb")
                )
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

        # Initialize SystemGenerator with GAFF2 + AM1-BCC
        templates = ["amber14-all.xml", "amber14/tip3p.xml"]
        system_generator = SystemGenerator(
            forcefields=templates,
            small_molecule_forcefield="gaff-2.11",
            molecules=[ligand_mol],
            periodic_forcefield_kwargs={
                "constraints": HBonds,
                "rigidWater": True,
                "nonbondedMethod": PME,
                "nonbondedCutoff": 1.0 * unit.nanometer,
            },
        )

        print("  • Solvating system and adding neutralizing ions…")
        modeller.addSolvent(
            system_generator.forcefield,
            model="tip3p",
            padding=1.0 * unit.nanometer,
            ionicStrength=0.15 * unit.molar,
        )

        system = system_generator.create_system(modeller.topology)
    else:
        # Legacy protein-water setup pathway
        fixer = PDBFixer(filename=cleaned_pdb)
        fixer.removeHeterogens(keepWater=True)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(pH=7.0)

        modeller = Modeller(fixer.topology, fixer.positions)
        ff = ForceField("amber14-all.xml", "amber14/tip3p.xml")

        print("  • Solvating system and adding neutralizing ions…")
        modeller.addSolvent(
            ff, model="tip3p", padding=1.0 * unit.nanometer, ionicStrength=0.15 * unit.molar
        )

        system = ff.createSystem(
            modeller.topology,
            nonbondedMethod=PME,
            nonbondedCutoff=1.0 * unit.nanometer,
            constraints=HBonds,
        )

    system.addForce(MonteCarloBarostat(1 * unit.atmosphere, 300 * unit.kelvin, 25))

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
    print("✔ Minimization complete (em.chk saved)")
    return True


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("workdir")
    p.add_argument("--pdbid", required=True)
    p.add_argument("--ligand", default=None)
    p.add_argument("--keep-mg", action="store_true")
    args = p.parse_args()
    run_em(args.workdir, args.pdbid, ligand=args.ligand, keep_mg=args.keep_mg)
