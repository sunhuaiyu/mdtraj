# working folder: 
# hsun@br1:/data1/hsun/SUMO3_20170928
# This is a 400ns trajectory of SUMO3, with 4000 frames
# All frames are aligned using VMD/RMSD visualization tools and
# saved to pdb file SUMO3_2mp2_md.pdb

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy
from scipy.spatial.distance import squareform

# This takes a few min to load; file size 198563303

#tr = md.load('SUMO1_2asq_md400ns_chainA.pdb')
#tr.save_hdf5('SUMO1_2asq_md400ns_chainA.h5')
tr = md.load('SUMO1_2asq_md400ns_chainA.h5')

# calculate pairwise distance of atoms 89-298 among all trajectories
# Atoms 89-298 are all *heavy* atoms of Ser32-Gly56 on SUMO1, SIM-binding site
# Hydrogen atoms are removed in the traj file (vmd: chain A and noh)
atoms_pick = range(89, 299)
tr.superpose(tr, frame=0, atom_indices=atoms_pick)  # do this everytime


# This takes a few min to calculate; will need multi CPU
distances = np.empty((tr.n_frames, tr.n_frames))
for i in range(tr.n_frames):
    distances[i] = md.rmsd(tr, tr, i, atom_indices=atoms_pick)

print('Max pairwise rmsd: %f nm' % np.max(distances))
'''Max pairwise rmsd: 0.366354 nm'''
# This number is the same regardless of superposing
'''
distances.tofile('RMSD_dist_SUMO1_S32-G56.txt', sep='\t')

distances = np.fromfile('RMSD_dist_SUMO1_S32-G56.txt', 
                        sep='\t').reshape((tr.n_frames, tr.n_frames))
'''


#RMSD over time
plot(distances[0])
xlabel('Frame'); ylabel('RMSD w.r.t. Frame[0] (nm)')
savefig('RMSD_SUMO1_S32-G56_trajectory.jpg', dpi=300)



# Clustering only accepts reduced form. Squareform's checks are too stringent
# assert np.all(distances - distances.T < 1e-6)
reduced_distances = squareform(distances, checks=False)
linkage_matrix = scipy.cluster.hierarchy.ward(reduced_distances)
# linkage_matrix encode clustering result
# used Ward's metrics

# dendrogram
plt.title('RMSD average linkage hierarchical clustering')
scipy.cluster.hierarchy.set_link_color_palette(['grey', 'b', 'c', 'g', 'k', 
                                                'm', 'r', 'y', 'indigo'])
                                                
hierarchy = scipy.cluster.hierarchy.dendrogram(
        linkage_matrix, 
        no_labels=True, 
        color_threshold=1.5, 
        count_sort='descendent',
        above_threshold_color='indigo')
yticks([])
savefig('linkage_RMSD_clustering.jpg', dpi=300)
close()



# Which cluster does each frame belong to? 
labels = scipy.cluster.hierarchy.fcluster(linkage_matrix, t=5, 
                                          criterion='maxclust')

# assign color to each frame on the RMSD plot
colors = ['grey', 'b', 'c', 'g', 'k', 'm', 'r', 'y']
for i in range(tr.n_frames):
    plot(i, distances[0][i], '.', markersize=2, color=colors[labels[i]-1])
xlabel('Frame No.'); ylabel('RMSD w.r.t. Frame[0] (nm)')
savefig('RMSD_traj_colored.jpg', dpi=300)
close()

# How many frames are there in each cluster?
n_clusters, first_idx, counts = np.unique(labels, return_counts=True, 
                                          return_index=True)

'''
(array([1, 2, 3, 4, 5], dtype=int32),
 array([ 118, 1459,    0, 1350, 3548]),
 array([1334, 1531,  118,  566,  453]))
'''

# compute the centroid frame of each cluster
clusters_frame = {i:[] for i in n_clusters}
# sort all frames into clusters
for i in range(tr.n_frames):
    clusters_frame[labels[i]].append(i)

all_centroids = [clusters_frame[c][np.exp(- distances[clusters_frame[c]] / 
                    distances[clusters_frame[c]].std()).sum(axis=1).argmax()] 
                for c in n_clusters ]
'''[543, 3300, 98, 1809, 3727]'''


# take the atomic slice of each centroid frame
SIM_site = tr[all_centroids].atom_slice(atoms_pick)
SIM_site.superpose(SIM_site, frame=0)
SIM_site.save_pdb('SUMO1_2asq_md400ns_SIM_site_rep.pdb')

for i in all_centroids:
    tr[i].save_pdb('SUMO1_2asq_md400ns_frame{}.pdb'.format(i))


# 20181217: 
# determine how the frames in original 2asq (NMR, 10 models) 
# fit in the clustered trajectory frames
# 2asq processed in VMD, saving a pdb with "chain A and noh"
# paste the same unitcell from MD trajectory on top of the pdb file
# named SUMO1_2asq_chainA.pdb

tr2 = md.load('SUMO1_2asq_chainA.pdb')

# join to the original traj
tr = md.load('SUMO1_2asq_md400ns_chainA.h5')
tr += tr2 

distances = np.empty((tr.n_frames, tr.n_frames))
for i in range(tr.n_frames):
    distances[i] = md.rmsd(tr, tr, i, atom_indices=atoms_pick)

reduced_distances = squareform(distances, checks=False)
linkage_matrix = scipy.cluster.hierarchy.ward(reduced_distances)
plt.title('RMSD average linkage hierarchical clustering')
scipy.cluster.hierarchy.set_link_color_palette(['grey', 'b', 'c', 'g', 'k', 
                                                'm', 'r', 'y', 'indigo'])
                                                
hierarchy = scipy.cluster.hierarchy.dendrogram(
        linkage_matrix, 
        no_labels=True, 
        color_threshold=1.5, 
        count_sort='descendent',
        above_threshold_color='indigo')
yticks([])

labels = scipy.cluster.hierarchy.fcluster(linkage_matrix, t=5, 
                                          criterion='maxclust')

# assign color to each frame on the RMSD plot
colors = ['grey', 'b', 'c', 'g', 'k', 'm', 'r', 'y']
for i in range(tr.n_frames):
    plot(i, distances[0][i], '.', markersize=2, color=colors[labels[i]-1])
xlabel('Frame No.'); ylabel('RMSD w.r.t. Frame[0] (nm)')

n_clusters, first_idx, counts = np.unique(labels, return_counts=True, 
                                          return_index=True)


'''
(array([1, 2, 3, 4, 5], dtype=int32),
 array([ 118, 1459,    0, 1350, 3548]),
 array([1334, 1531,  118,  576,  453])) 
'''
# the 10 models are assigned to only 1 cluster (No. 4) 
# cluster #4 used to contain 566 frames




