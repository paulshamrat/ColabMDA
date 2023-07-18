#Rg
# load the trajectory and topology files
u = mda.Universe('/content/gdrive/MyDrive/1aki/1AKI_solv_ions.gro', '/content/gdrive/MyDrive/1aki/md_0_1.xtc')

# select protein atoms for analysis
protein = u.select_atoms('protein')

# calculate the radius of gyration
com = np.array([protein.center_of_mass()])
Rg_list = []
time_list = []
for ts in u.trajectory:
    Rg = np.sqrt(np.sum((protein.positions - com)**2)/len(protein))
    Rg_list.append(Rg)
    time_list.append(ts.time)

# plot the radius of gyration
plt.plot(time_list, Rg_list)
plt.xlabel('Time (ps)')
plt.ylabel('Radius of gyration (Å)')

# save the plot as a PNG file
plt.savefig('rg_plot.png', dpi=300)

# show the plot
plt.show()
