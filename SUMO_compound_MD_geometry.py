# geometry at the SIM binding site of SUMO1

import mdtraj as md  #mdtraj 1.9
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from scipy.spatial import ConvexHull
from glob import glob
from multiprocessing import Pool

def convexhull_volume(traj, SUMO):
    '''
    Compute convexhull volume of given residues along a md trajectory
    traj: mdtraj trajectory
    SUMO: integer 1 or 2
    return a numpy array of volumes calculated from each frame
    '''
    
    # SUMO1: I34, F36, K37, K39, K46, Y51, R54
    # SUMO2:
    refseq_num = {1: ['ILE34', 'PHE36', 'LYS37', 'LYS39', 'LYS46', 'TYR51', 'ARG54'],
                  2: ['VAL30', 'PHE32', 'LYS33', 'LYS35', 'LYS42', 'TYR47', 'ARG50'] }

    top = traj.topology
    xyz = traj.xyz #coordinates of all atoms
    selected_residues = refseq_num[SUMO]
    
    # mapping RefSeq residue names to mdtraj residue indices
    residues = [str(i) for i in top.residues]
    residue_indices = [residues.index(i) for i in selected_residues] 

    # obtain atom_indices
    atom_indices = [top.select('resid == {}'.format(i)) for i in residue_indices]
    atom_indices = np.concatenate(atom_indices)

    return np.array([ConvexHull(i).volume for i in xyz[:, atom_indices, :]])


def SUMO_ligand_dist(traj):
    '''
    given a SUMO1-compound trajectory, 
    compute distance between center of compound and F36CG
    '''
    #coordinates for the Cgamma of SUMO1_F36, SUMO2/3_F31
    atom_ix = traj.topology.select('residue==36 and name==CG')[0]
    a = traj.xyz[:, atom_ix]

    # ligand all atom coordinatess:
    lig = traj.atom_slice(traj.topology.select('chainid==1'))
    # ligand center of mass: 
    b = md.compute_center_of_mass(lig)

    return (((a - b) ** 2).sum(1)) ** 0.5


def atom_pair_ix(traj, pair='F36CG_R54CZ'):
    top = traj.topology
    s = pair.split('_')
    pair_ix = top.select_pairs('residue=={0} and name=={1}'.format(s[0][1:3], s[0][3:]),
                               'residue=={0} and name=={1}'.format(s[1][1:3], s[1][3:]))
    return pair_ix

def atom_pair_dist2(traj_list, pair='F36CG_R54CZ'):
    pair_ix = atom_pair_ix(traj_list[0], pair=pair)    
    dist = array([md.compute_distances(i, atom_pairs=pair_ix, periodic=False).ravel() 
                    for i in traj_list])
    return dist
    
    

T = 0.7
# Histograms of atomic distance involving trajectories bound with SIM, A12 and B27
# dist_dict computed in SUMO_compound_mdtraj_analysis.py

A12, B27 = 'PHG00686', 'SEW05414'

pair = 'F36CG_Y51CG'
#pair = 'K37NZ_R54CZ'
#pair = 'F36CG_R54CZ'


tr2uyz = [md.load('SUMO1_2uyz_{}_400ns.h5'.format(i+1))[::10] for i in range(12)]
dist2uyz = atom_pair_dist2(tr2uyz, pair).ravel()
vol2uyz = np.array([convexhull_volume(i, 1) for i in tr2uyz]).ravel()

F36_A12 = [dist_dict['SUMO1_2uyz_{0}_F{1}_5000ns.h5'.format(A12, i+1)] for i in range(12)]
trA12 = md.join([traj_dict['SUMO1_2uyz_{0}_F{1}_5000ns.h5'.format(A12, i+1)][F36_A12[i] < T] 
                   for i in range(12)])
distA12 = atom_pair_dist2([trA12], pair).ravel()
volA12 = convexhull_volume(trA12, 1)

F36_B27 = [dist_dict['SUMO1_2uyz_{0}_F{1}_5000ns.h5'.format(B27, i+1)] for i in range(12)]
trB27 = md.join([traj_dict['SUMO1_2uyz_{0}_F{1}_5000ns.h5'.format(B27, i+1)][F36_B27[i] < T] 
                   for i in range(12)])
distB27 = atom_pair_dist2([trB27], pair).ravel()
volB27 = convexhull_volume(trB27, 1)

trSIM = [md.load('SUMO1_2asq_wSIM_{}_400ns.h5'.format(i+1)) for i in range(9)]
distSIM = atom_pair_dist2(trSIM, pair).ravel()
volSIM = np.array([convexhull_volume(i, 1) for i in  trSIM]).ravel()


figure(figsize=(3.6, 4.8))
hist(distSIM, color='tab:blue', alpha=0.6, bins=100, linewidth=1, orientation='horizontal')
hist(distA12, color= 'gray', alpha=0.6, bins=100, linewidth=1, orientation='horizontal')
legend(['SIM', 'A12'], fontsize=15, frameon=False)
name = 'SUMO1_SIM_A12_' + pair
#title(name);  
#ylim(0.4, 1.4) #pair = 'F36CG_Y51CG'
#xlim(0, 1500); ylim(0.5, 3.0) #pair = 'K37NZ_R54CZ'
tight_layout()
savefig(name + '.jpg', dpi=600) 


figure(figsize=(3.6, 4.8))
hist(distSIM, color='tab:blue', alpha=0.6, bins=100, linewidth=1, orientation='horizontal')
hist(distB27, color= 'gray', alpha=0.6, bins=100, linewidth=1, orientation='horizontal')
legend(['SIM', 'B27'], fontsize=15, frameon=False)
name = 'SUMO1_SIM_B27_' + pair
#title(name); 
#ylim(0.4, 1.4) #pair = 'F36CG_Y51CG'
#xlim(0, 1500); ylim(0.5, 3.0) #pair = 'K37NZ_R54CZ'
tight_layout()
savefig(name + '.jpg', dpi=600) 



figure(figsize=(3.6, 4.8))
hist(volSIM, color='tab:blue', alpha=0.6, bins=100, linewidth=1, orientation='horizontal')
hist(volA12, color= 'gray', alpha=0.6, bins=100, linewidth=1, orientation='horizontal')
legend(['SIM', 'A12'], fontsize=15, frameon=False)
name = 'SUMO1_SIM_A12_convexhull' 
#title(name);  
ylim(1.2, 2.4)
tight_layout()
savefig(name + '.jpg', dpi=600) 


figure(figsize=(3.6, 4.8))
hist(volSIM, color='tab:blue', alpha=0.6, bins=100, linewidth=1, orientation='horizontal')
hist(volB27, color= 'gray', alpha=0.6, bins=100, linewidth=1, orientation='horizontal')
legend(['SIM', 'B27'], fontsize=15, frameon=False)
name = 'SUMO1_SIM_B27_convexhull'
#title(name); 
ylim(1.2, 2.4)
tight_layout()
savefig(name + '.jpg', dpi=600) 

