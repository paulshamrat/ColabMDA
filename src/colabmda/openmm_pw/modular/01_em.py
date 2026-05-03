import os, sys, argparse
from pdbfixer import PDBFixer
from openmm.app import Modeller, ForceField, PME, HBonds, PDBFile
from openmm import XmlSerializer, unit, MonteCarloBarostat
from colabmda.openmm_pw.modular.utils import make_sim

def run_em(workdir, pdbid):
    os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)
    
    cleaned_pdb = os.path.join(workdir, f"{pdbid}_cleaned.pdb")
    if not os.path.exists(cleaned_pdb):
        print(f"Error: cleaned PDB not found: {cleaned_pdb}")
        return False

    xml_system = 'system.xml'
    pdb_saved  = 'solvated.pdb'
    chk_file   = 'em.chk'

    print("▶ Starting System Preparation & Minimization …")
    fixer = PDBFixer(filename=cleaned_pdb)
    fixer.removeHeterogens(keepWater=True)
    fixer.findMissingResidues(); fixer.findMissingAtoms()
    fixer.addMissingAtoms(); fixer.addMissingHydrogens(pH=7.0)

    modeller = Modeller(fixer.topology, fixer.positions)
    ff = ForceField('amber14-all.xml','amber14/tip3p.xml')
    modeller.addSolvent(ff, model='tip3p', padding=1.0*unit.nanometer, ionicStrength=0.15*unit.molar)

    system = ff.createSystem(modeller.topology,
                             nonbondedMethod=PME,
                             nonbondedCutoff=1.0*unit.nanometer,
                             constraints=HBonds)
    system.addForce(MonteCarloBarostat(1*unit.atmosphere, 300*unit.kelvin, 25))

    dt = 2.0*unit.femtoseconds
    sim = make_sim(modeller.topology, system, dt)
    sim.context.setPositions(modeller.positions)

    print("  • Minimizing…")
    sim.minimizeEnergy()
    
    # Save files
    with open(xml_system,'w') as f:
        f.write(XmlSerializer.serializeSystem(sim.context.getSystem()))
    with open(pdb_saved,'w') as f:
        PDBFile.writeFile(sim.topology,
                          sim.context.getState(getPositions=True).getPositions(),
                          f)
    sim.saveCheckpoint(chk_file)
    print("✔ Minimization complete (em.chk saved)")
    return True

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("workdir")
    p.add_argument("--pdbid", required=True)
    args = p.parse_args()
    run_em(args.workdir, args.pdbid)
