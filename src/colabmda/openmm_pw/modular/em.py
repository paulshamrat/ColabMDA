import argparse
import os

from openmm import MonteCarloBarostat, XmlSerializer, unit
from openmm.app import PME, ForceField, HBonds, Modeller, PDBFile
from pdbfixer import PDBFixer

from colabmda.openmm_pw.modular.utils import make_sim


def run_em(workdir, pdbid):
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
    else:
        print(
            f"[DEBUG] Files not found. em.chk: {os.path.exists(chk_file)}, system.xml: {os.path.exists(xml_system)}"
        )
        if os.path.exists(workdir):
            print(f"[DEBUG] Directory contents: {os.listdir(workdir)}")

    print("▶ Starting System Preparation & Minimization …")
    fixer = PDBFixer(filename=cleaned_pdb)
    fixer.removeHeterogens(keepWater=True)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.0)

    modeller = Modeller(fixer.topology, fixer.positions)
    ff = ForceField("amber14-all.xml", "amber14/tip3p.xml")
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
    args = p.parse_args()
    run_em(args.workdir, args.pdbid)
