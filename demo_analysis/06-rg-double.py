# RG
import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

# Load the trajectory and topology files for both systems
u1 = mda.Universe('/content/gdrive/MyDrive/1aki/1AKI_solv_ions.gro', '/content/gdrive/MyDrive/1aki/md_0_1.xtc')
u2 = mda.Universe('/content/gdrive/MyDrive/1aki/1AKI_solv_ions.gro', '/content/gdrive/MyDrive/1aki/md_0_1.xtc')

# Select only protein atoms
protein_sel1 = u1.select_atoms('protein')
protein_sel2 = u2.select_atoms('protein')

# Initialize arrays to store Rg values and time
Rg1 = np.zeros(len(u1.trajectory))
Rg2 = np.zeros(len(u2.trajectory))
time1 = np.zeros(len(u1.trajectory))
time2 = np.zeros(len(u2.trajectory))

# Loop over all frames in trajectory and calculate Rg
for ts in u1.trajectory:
    Rg1[ts.frame] = protein_sel1.radius_of_gyration()
    time1[ts.frame] = u1.trajectory.time

for ts in u2.trajectory:
    Rg2[ts.frame] = protein_sel2.radius_of_gyration()
    time2[ts.frame] = u2.trajectory.time

# Plot Rg values of both systems on the same plot
fig, ax = plt.subplots()
ax.plot(time1, Rg1, label='data_1')
ax.plot(time2, Rg2, label='data_2')
ax.set_xlabel('Time (ps)')
ax.set_ylabel('Radius of gyration (nm)')
ax.legend()

# Save the plot as a high-resolution PNG image
fig.savefig('rg_230425.png', dpi=300)

plt.show()
