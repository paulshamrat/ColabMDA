import MDAnalysis as mda
from MDAnalysis.analysis import rms
import matplotlib.pyplot as plt

# load the trajectories and topologies
u1 = mda.Universe('/content/gdrive/MyDrive/1aki/1AKI_solv_ions.gro', '/content/gdrive/MyDrive/1aki/md_0_1.xtc')
u2 = mda.Universe('/content/gdrive/MyDrive/1aki/1AKI_solv_ions.gro', '/content/gdrive/MyDrive/1aki/md_0_1.xtc')

# select the protein atoms for RMSD calculation
ref1 = u1.select_atoms('protein')
ref2 = u2.select_atoms('protein')

# calculate the RMSD
R1 = rms.RMSD(u1, ref1, select='backbone', groupselections=['protein'])
R1.run()
R2 = rms.RMSD(u2, ref2, select='backbone', groupselections=['protein'])
R2.run()

# plot the RMSD
fig, ax = plt.subplots()
ax.plot(R1.rmsd[:,1], R1.rmsd[:,2], label='data_1')
ax.plot(R2.rmsd[:,1], R2.rmsd[:,2], label='data_2')
ax.legend()
ax.set_xlabel('Time (ps)')
ax.set_ylabel('RMSD (Å)')

# save the plot as a PNG file
plt.savefig('rmsd_plot_230425.png', dpi=300)

# show the plot
plt.show()
