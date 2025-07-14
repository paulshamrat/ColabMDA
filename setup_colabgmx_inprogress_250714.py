#!/usr/bin/env python3
"""
setup_gromacs_colab.py

Automates installation of GROMACS 2023.x on Google Drive via Google Colab.
"""

import os
import subprocess
import sys

from google.colab import drive

def run(cmd, check=True):
    """Helper to run a shell command."""
    print(f"\n$ {cmd}")
    result = subprocess.run(cmd, shell=True, check=check, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(result.stdout)
    return result

def main():
    # 1. Mount Google Drive
    print("Mounting Google Drive…")
    drive.mount('/content/drive', force_remount=True)
    
    # 2. Create installation directory and cd into it
    gmx_dir = "/content/drive/MyDrive/gmx"
    print(f"Creating directory {gmx_dir} and changing into it…")
    os.makedirs(gmx_dir, exist_ok=True)
    os.chdir(gmx_dir)
    print("Current directory:", os.getcwd())
    
    # 3. Check GPU
    run("nvidia-smi")
    
    # 4. Ensure a compatible CMake version (downgrading if needed)
    print("Installing / downgrading CMake to v3.26.4…")
    run("pip install --upgrade cmake")
    run("wget -q https://github.com/Kitware/CMake/releases/download/v3.26.4/cmake-3.26.4-linux-x86_64.sh")
    run("chmod +x cmake-3.26.4-linux-x86_64.sh")
    run("./cmake-3.26.4-linux-x86_64.sh --skip-license --prefix=/usr/local")
    run("cmake --version")
    
    # 5. Install system dependencies
    print("Updating apt and installing dependencies…")
    run("apt update -y && apt upgrade -y")
    run("apt install -y gcc build-essential libfftw3-dev nvidia-cuda-toolkit wget tar")
    
    # 6. Build and install hwloc (for topology queries)
    print("Downloading and installing hwloc 1.11.13…")
    run("wget -q https://download.open-mpi.org/release/hwloc/v1.11/hwloc-1.11.13.tar.gz")
    run("tar xzf hwloc-1.11.13.tar.gz")
    os.chdir("hwloc-1.11.13")
    run("./configure")
    run("make -j4")
    run("make install")
    os.chdir(gmx_dir)
    
    # 7. Download GROMACS source
    gmx_version = "2023.1"
    print(f"Downloading GROMACS {gmx_version}…")
    run(f"wget -q ftp://ftp.gromacs.org/gromacs/gromacs-{gmx_version}.tar.gz")
    run(f"tar xfz gromacs-{gmx_version}.tar.gz")
    
    # 8. Configure GROMACS build
    build_dir = f"{gmx_dir}/gromacs-{gmx_version}/build"
    print(f"Configuring GROMACS build in {build_dir}…")
    os.makedirs(build_dir, exist_ok=True)
    cmake_cmd = (
        f"cmake -B {build_dir} "
        f"-S {gmx_dir}/gromacs-{gmx_version} "
        "-DGMX_BUILD_OWN_FFTW=ON "
        "-DREGRESSIONTEST_DOWNLOAD=ON "
        "-DGMX_GPU=CUDA "
        "-DGMX_MPI=OFF "
        "-DGMX_OPENMP=OFF "
        "-DGMX_DOUBLE=OFF "
        f"-DCMAKE_INSTALL_PREFIX={gmx_dir}/gromacs "
        "-DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda "
        "-DCMAKE_C_COMPILER=gcc "
        "-DCMAKE_CXX_COMPILER=g++"
    )
    run(cmake_cmd)
    
    # 9. Build and install GROMACS
    print("Building and installing GROMACS…")
    run(f"cmake --build {build_dir} -j4")
    run(f"cmake --build {build_dir} --target install")
    
    # 10. Verify installation
    print("Sourcing GMXRC and checking gmx command…")
    verify_script = f"source {gmx_dir}/gromacs/bin/GMXRC && gmx --version"
    subprocess.call(verify_script, shell=True, executable="/bin/bash")

if __name__ == "__main__":
    main()
