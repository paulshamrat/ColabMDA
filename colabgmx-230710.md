
Pre-compiled GROMACS installation

# Download
!wget https://github.com/paulshamrat/ColabMDA/raw/main/gromacs2021.5.tar

# successful untar
import subprocess

# Specify the tar file path
tar_file_path = '/content/gromacs2021.5.tar'

# Specify the target directory for untarring
target_directory = '/content'

# Untar the file
subprocess.run(['tar', '-xf', tar_file_path, '-C', target_directory])

print("Tar file successfully untarred.")


# Checking that our GROMACS works
%%bash
source /content/gromacs/bin/GMXRC
gmx -h

# We will be needed biopython for pdb download
! pip install biopython

Mount Google Drive

# Mount google drive
from google.colab import drive
drive.mount('/content/gdrive')


# Change Working Directory
import os
dir = "/content/gdrive/MyDrive/works/1aki/"
if not os.path.exists(dir):
    os.makedirs(dir)
%cd /content/gdrive/MyDrive/works/1aki/1aki_1/
!pwd


# Change Working Directory
import os
dir = "/content/gdrive/MyDrive/works/1aki/"
if not os.path.exists(dir):
    os.makedirs(dir)
%cd /content/gdrive/MyDrive/works/1aki/
!pwd

Run MD Simulation

# Download the mdp file first and make necessary changes in mdp parameters
%%bash
source /content/gromacs/bin/GMXRC
wget -O ions.mdp https://raw.githubusercontent.com/paulshamrat/ColabMDA/main/mdps/ions.mdp
wget -O md.mdp https://raw.githubusercontent.com/paulshamrat/ColabMDA/main/mdps/md.mdp
wget -O minim.mdp https://raw.githubusercontent.com/paulshamrat/ColabMDA/main/mdps/minim.mdp
wget -O npt.mdp https://raw.githubusercontent.com/paulshamrat/ColabMDA/main/mdps/npt.mdp
wget -O nvt.mdp https://raw.githubusercontent.com/paulshamrat/ColabMDA/main/mdps/nvt.mdp


#Importing your PDB file using biopython
import os
from Bio.PDB import *
pdbid = ['1aki']
pdbl = PDBList()
for s in pdbid:
  pdbl.retrieve_pdb_file(s, pdir='.', file_format ="pdb", overwrite=True)
  os.rename("pdb"+s+".ent", s+".pdb")



# We will constantly need to source GMXRC for GROMACS to work
%%bash
source /content/gromacs/bin/GMXRC

# Strip out Water
grep -v HOH 1aki.pdb > 1AKI_clean.pdb

# cleaned file
gmx pdb2gmx -f 1AKI_clean.pdb -o 1AKI_processed.gro -water spce -ignh -ff amber99sb-ildn

# add box
gmx editconf -f 1AKI_processed.gro -o 1AKI_newbox.gro -c -d 1.0 -bt cubic

# solvate
gmx solvate -cp 1AKI_newbox.gro -cs spc216.gro -o 1AKI_solv.gro -p topol.top

# final sove
gmx grompp -f ions.mdp -c 1AKI_solv.gro -p topol.top -o ions.tpr

#This is a trick to provide interactive options to gmx
echo "SOL" > options
echo " " >> options
# Now we have an atomic-level description of our system in the binary file ions.tpr. We will pass this file to genion:
# When prompted, choose group 13 "SOL" for embedding ions.
gmx genion -s ions.tpr -o 1AKI_solv_ions.gro -p topol.top -pname NA -nname CL -neutral < options

# Assemble the binary input using grompp using [](this) input parameter file:
gmx grompp -f minim.mdp -c 1AKI_solv_ions.gro -p topol.top -o em.tpr

#Using grompp to prepare our minimization MD
gmx grompp -f minim.mdp -c 1AKI_solv_ions.gro -p topol.top -o em.tpr

#Run our minimization
gmx mdrun -v -deffnm em -nb gpu -ntmpi 1

#This is a trick to provide interactive options to gmx
echo "Potential" > options
echo " " >> options

# Let's do a bit of analysis. The em.edr file contains all of the energy terms that GROMACS collects during EM.
# You can analyze any .edr file using the GROMACS energy module:
gmx energy -f em.edr -o potential.xvg < options
# At the prompt, type "10 0" to select Potential (10); zero (0) terminates input.

# We will call grompp and mdrun just as we did at the EM step:
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt -nb gpu -ntmpi 1

#This is a trick to provide interactive options to gmx
echo "Temperature" > options
echo " " >> options
# Let's analyze the temperature progression, again using energy:
gmx energy -f nvt.edr -o temperature.xvg < options
# Type "16 0" at the prompt to select the temperature of the system and exit.

# We will call grompp and mdrun just as we did for NVT equilibration.
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt -nb gpu -ntmpi 1


#This is a trick to provide interactive options to gmx
echo "Pressure" > options
echo " " >> options
# Let's analyze the pressure progression, again using energy:
gmx energy -f npt.edr -o pressure.xvg < options
# Type "18 0" at the prompt to select the pressure of the system and exit.


#This is a trick to provide interactive options to gmx
echo "Density" > options
echo " " >> options

# Let's take a look at density as well
gmx energy -f npt.edr -o density.xvg < options
# this time using energy and entering "24 0" at the prompt.




Production MD
# We will run a 1-ns MD simulation,
%%bash
source /content/gromacs/bin/GMXRC
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr

# Assuming you have one GPU available, the mdrun command to make use of it is as simple as:
%%bash
source /content/gromacs/bin/GMXRC
gmx mdrun -deffnm md_0_1 -nb gpu -ntmpi 1

Appending Simulation
%%bash
#Append
source /content/gromacs/bin/GMXRC
gmx mdrun -v -cpi md_1.cpt -noappend -deffnm md_1 -gpu_id 0


Concatenation of XTC files
%%bash
#Concatenation of xtc files
source /content/gromacs/bin/GMXRC
gmx trjcat -f md_1.xtc md_1.part0002.xtc md_1.part0003.xtc md_1.part0004.xtc md_1.part0005.xtc md_1.part0006.xtc md_1.part0007.xtc md_1.part0008.xtc md_1.part0009.xtc md_1.part0010.xtc md_1.part0011.xtc md_1.part0012.xtc md_1.part0013.xtc -o md_1_all.xtc


