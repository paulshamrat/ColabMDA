import os, sys, argparse
from openmm.app import PDBFile, StateDataReporter
from openmm import XmlSerializer, unit
from colabmda.openmm_pw.modular.utils import make_sim

def run_npt(workdir, pdbid, equil_time_ps=100.0):
    os.chdir(workdir)
    
    xml_system = 'system.xml'
    pdb_saved  = 'solvated.pdb'
    chk_in     = 'nvt.chk'
    chk_out    = 'npt.chk'
    log_file   = 'npt.log'

    if not all(os.path.exists(f) for f in [xml_system, pdb_saved, chk_in]):
        print("Error: Missing nvt.chk or system files. Run NVT first.")
        return False

    print(f"▶ Starting NPT Equilibration ({equil_time_ps} ps) …")
    with open(xml_system) as f:
        system = XmlSerializer.deserializeSystem(f.read())
    pdb = PDBFile(pdb_saved)
    
    dt = 2.0*unit.femtoseconds
    sim = make_sim(pdb.topology, system, dt)
    
    with open(chk_in, 'rb') as f:
        sim.loadCheckpoint(f)
    
    eq_steps = max(1, int((equil_time_ps*unit.picoseconds)/dt))
    sim.reporters.append(StateDataReporter(log_file, max(1, eq_steps//50), 
                                           step=True, time=True, temperature=True, 
                                           potentialEnergy=True, volume=True, density=True))
    
    sim.step(eq_steps)
    
    # Save final equilibrated PDB
    with open('equilibrated.pdb','w') as f:
        PDBFile.writeFile(sim.topology,
                          sim.context.getState(getPositions=True).getPositions(),
                          f)
                          
    sim.saveCheckpoint(chk_out)
    print("✔ NPT complete (npt.chk saved)")
    return True

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("workdir")
    p.add_argument("--pdbid", required=True)
    p.add_argument("--equil-time", type=float, default=100.0)
    args = p.parse_args()
    run_npt(args.workdir, args.pdbid, args.equil_time)
