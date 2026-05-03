import os, sys, shutil
from openmm.app import Simulation
from openmm import Platform, unit, LangevinMiddleIntegrator

def pick_platform():
    for name in ["CUDA", "OpenCL", "CPU"]:
        try:
            plat = Platform.getPlatformByName(name)
            # print(f"[OpenMM] Using platform: {name}")
            return plat
        except Exception:
            pass
    raise RuntimeError("No OpenMM platform available")

def make_sim(top, sys_, dt, temp=300):
    integ = LangevinMiddleIntegrator(temp*unit.kelvin, 1/unit.picosecond, dt)
    plat  = pick_platform()
    return Simulation(top, sys_, integ, plat)

def atomic_rename(tmp_path, final_path):
    if os.path.exists(tmp_path):
        os.replace(tmp_path, final_path)

def safe_remove_if_empty(path):
    if os.path.isfile(path) and os.path.getsize(path) == 0:
        try:
            os.remove(path)
        except Exception:
            pass

def sync_outputs(workdir, sync_dir, extra_files=None):
    if not sync_dir:
        return
    os.makedirs(sync_dir, exist_ok=True)
    essentials = ["system.xml", "solvated.pdb", "em.chk", "nvt.chk", "npt.chk", "prod.chk", "nvt.log", "npt.log", "prod_full.log"]
    if extra_files:
        essentials += extra_files
    for fn in essentials:
        src = os.path.join(workdir, fn)
        if os.path.exists(src):
            try:
                shutil.copy2(src, os.path.join(sync_dir, fn))
            except Exception:
                pass
