#====================================
# 01 Reading Data
#====================================

# Access Data
import MDAnalysis as mda
u = mda.Universe('/content/gdrive/MyDrive/1aki/6wzu_solv_ions.pdb', '/content/gdrive/MyDrive/1aki/md_1_all.xtc')

# Write Data
import MDAnalysis as mda
u = mda.Universe('/content/gdrive/MyDrive/1aki/6wzu_solv_ions.pdb', '/content/gdrive/MyDrive/1aki/md_1_all.xtc')
ag = u.select_atoms("name CA")
ag.write("c-alpha.pdb")

# Pass in the frames keyword to write out trajectories.
ag.write('c-alpha_all.xtc', frames='all')

# Slice or index the trajectory to choose which frames to write:
ag.write('c-alpha_skip2.trr', frames=u.trajectory[::2])
ag.write('c-alpha_some.dcd', frames=u.trajectory[[0,2,3]])

# make xtc to dcd with all frames
ag.write('c-alpha_all.trr', frames='all')
ag.write('c-alpha_all.dcd', frames='all')

# Alternatively, iterate over the trajectory frame-by-frame with Writer(). This requires you to pass in the number of atoms to write.
with mda.Writer('c-alpha.xyz', ag.n_atoms) as w:
    for ts in u.trajectory:
        w.write(ag)


import mdtraj as md
t = md.load('c-alpha_all.xtc', top='c-alpha.pdb')
print(t)

# Save as h5 files
atoms_to_keep = [a.index for a in t.topology.atoms if a.name == 'CA']
t.restrict_atoms(atoms_to_keep)  # this acts inplace on the trajectory
t.save('CA-only.h5')





#====================================
# 02 Atom Selection
#====================================

# Access h5 files
from __future__ import print_function
import mdtraj as md

traj = md.load('CA-only.h5')
print(traj)

# Atom count
print('How many atoms?    %s' % traj.n_atoms)
print('How many residues? %s' % traj.n_residues)

# idx
frame_idx = 4 # zero indexed frame number
atom_idx = 9 # zero indexed atom index
print('Where is the fifth atom at the tenth frame?')
print('x: %s\ty: %s\tz: %s' % tuple(traj.xyz[frame_idx, atom_idx,:]))

# topology
topology = traj.topology
print(topology)

# atom count
print('Fifth atom: %s' % topology.atom(4))
print('All atoms: %s' % [atom for atom in topology.atoms])

# print residue
print('Second residue: %s' % traj.topology.residue(1))
print('All residues: %s' % [residue for residue in traj.topology.residues])

# topology
atom = topology.atom(10)
print('''Hi! I am the %sth atom, and my name is %s.
I am a %s atom with %s bonds.
I am part of an %s residue.''' % ( atom.index, atom.name, atom.element.name, atom.n_bonds, atom.residue.name))






#====================================
# 03 Put Everything Together 
#====================================

print([atom.index for atom in topology.atoms if atom.element.symbol is 'C' and atom.is_sidechain])
print([residue for residue in topology.chain(0).residues if residue.index % 2 == 0])





#====================================
# 04 Atom Selection Language
#====================================

# topology
print(topology.select('resid 1 to 2'))
print(topology.select('name N and backbone'))

selection = topology.select_expression('name CA and resid 1 to 2')
print(selection)





#====================================
# 05 Clustering
#====================================

# import clustering modules
from __future__ import print_function
%matplotlib inline
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy
from scipy.spatial.distance import squareform

# load h5 file
traj = md.load('CA-only.h5')

# compute pairwise rmsd between conformations
distances = np.empty((traj.n_frames, traj.n_frames))
for i in range(traj.n_frames):
    distances[i] = md.rmsd(traj, traj, i)
print('Max pairwise rmsd: %f nm' % np.max(distances))


# Clustering only accepts reduced form. Squareform's checks are too stringent
# when calculating a massinve numer of frames initially it was 1e-6; changed to 1e6 and the plot generated.
assert np.all(distances - distances.T < 1e6)
reduced_distances = squareform(distances, checks=False)

# linkage
linkage = scipy.cluster.hierarchy.linkage(reduced_distances, method='average')

# Plotting
plt.title('RMSD Average linkage hierarchical clustering')
_ = scipy.cluster.hierarchy.dendrogram(linkage, no_labels=True, count_sort='descendent')

# save the plot as a PNG file
plt.savefig('RMSD Average linkage hierarchical clustering.png', dpi=300)

# show the plot
plt.show()





#====================================
# 05 PCA: Principal Component Analysis
#====================================

# importing things
%matplotlib inline
from __future__ import print_function
import mdtraj as md
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

# load h5 files
traj = md.load('CA-only.h5')
traj

# components
pca1 = PCA(n_components=2)
traj.superpose(traj, 0)

# cartesian
reduced_cartesian = pca1.fit_transform(traj.xyz.reshape(traj.n_frames, traj.n_atoms * 3))
print(reduced_cartesian.shape)





#====================================
# 6.1 Cartesian Coordinate PCA
#====================================

# Cartesian Coordinate PCA
plt.figure()
plt.scatter(reduced_cartesian[:, 0], reduced_cartesian[:,1], marker='x', c=traj.time)
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title('Cartesian coordinate PCA')
cbar = plt.colorbar()
cbar.set_label('Time [ps]')

# save the plot as a PNG file
plt.savefig('Cartesian coordinate PCA', dpi=300)

# show the plot
plt.show()






#====================================
# 6.2 Pairwise Distance PCA
#====================================

# Pairwise Distance PCA
pca2 = PCA(n_components=2)

from itertools import combinations
# this python function gives you all unique pairs of elements from a list

atom_pairs = list(combinations(range(traj.n_atoms), 2))
pairwise_distances = md.geometry.compute_distances(traj, atom_pairs)
print(pairwise_distances.shape)
reduced_distances = pca2.fit_transform(pairwise_distances)

# Save plot
plt.figure()
plt.scatter(reduced_distances[:, 0], reduced_distances[:,1], marker='x', c=traj.time)
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title('Pairwise distance PCA')
cbar = plt.colorbar()
cbar.set_label('Time [ps]')

# save the plot as a PNG file
plt.savefig('Pairwise distance PCA', dpi=300)

# show the plot
plt.show()





#====================================
# 7. SASA
#====================================


# Import Libraries
from __future__ import print_function
%matplotlib inline
import numpy as np
import mdtraj as md

# shrake
help(md.shrake_rupley)

# sasa data shape
trajectory = md.load('CA-only.h5')
sasa = md.shrake_rupley(trajectory)

print(trajectory)
print('sasa data shape', sasa.shape)


# total sasa
total_sasa = sasa.sum(axis=1)
print(total_sasa.shape)

# total sasa
from matplotlib.pylab import *

plot(trajectory.time, total_sasa)
xlabel('Time [ps]', size=16)
ylabel('Total SASA (nm)^2', size=16)


# save the plot as a PNG file
plt.savefig('Total_SASA.png', dpi=300)

# show the plot
plt.show()

# sasa autocorrelation
def autocorr(x):
    "Compute an autocorrelation with numpy"
    x = x - np.mean(x)
    result = np.correlate(x, x, mode='full')
    result = result[result.size//2:]
    return result / result[0]

semilogx(trajectory.time, autocorr(total_sasa))
xlabel('Time [ps]', size=16)
ylabel('SASA autocorrelation', size=16)
show()