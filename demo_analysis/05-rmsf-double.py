# RMSF -CA
import MDAnalysis as mda
from MDAnalysis.analysis import rms
from MDAnalysis.analysis.rms import RMSF
import matplotlib.pyplot as plt

# Load the two trajectory and topology files
u1 = mda.Universe('/content/gdrive/MyDrive/1aki/1AKI_solv_ions.gro', '/content/gdrive/MyDrive/1aki/md_0_1.xtc')
u2 = mda.Universe('/content/gdrive/MyDrive/1aki/1AKI_solv_ions.gro', '/content/gdrive/MyDrive/1aki/md_0_1.xtc')

# Select the C-alpha atoms
calpha1 = u1.select_atoms('protein and name CA')
calpha2 = u2.select_atoms('protein and name CA')

# Align the protein to the reference structure
ref1 = mda.Universe('/content/gdrive/MyDrive/1aki/1AKI_solv_ions.gro')
ref2 = mda.Universe('/content/gdrive/MyDrive/1aki/1AKI_solv_ions.gro')
R1 = rms.RMSD(u1, ref1, select='protein and name CA', center=True, superposition=True)
R1.run()
R2 = rms.RMSD(u2, ref2, select='protein and name CA', center=True, superposition=True)
R2.run()

# Calculate RMSF for each trajectory
RMSF1 = RMSF(calpha1).run()
RMSF2 = RMSF(calpha2).run()

# Plot RMSF values of both trajectories on the same plot
fig, ax = plt.subplots()
ax.plot(RMSF1.rmsf, label='data_1')
ax.plot(RMSF2.rmsf, label='data_2')
ax.set_xlabel('Residue')
ax.set_ylabel('RMSF (nm)')
ax.legend()

# Set the y-axis limits based on the range of your RMSF values
ymin = 0
ymax = max(RMSF1.rmsf.max(), RMSF2.rmsf.max()) + 0.1
ax.set_ylim(ymin, ymax)

# Save the plot as a high-resolution PNG image
fig.savefig('rmsf_ca_230425.png', dpi=300)

plt.show()
