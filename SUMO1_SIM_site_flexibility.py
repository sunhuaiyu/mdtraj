# to estimate the volume of SIM-docking site on SUMO1 throughout an MD trajectory

import mdtraj as md  #mdtraj 1.9
import numpy as np
import matplotlib.pyplot as plt
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


def atom_pair_dist(traj, SUMO, pair='36-CG, 54-NE'):
    '''
    given an MD trajectory of SUMO, 
    return a numpy array of distances -- e.g.:
    between SUMO1's [F36CG, Y51CG] and [K37NZ and R54NE], or
    between SUMO2's [F32CG, Y47CG] and [K33NZ and R50NE]
    '''
    pairs = {1:[[[36, 'CG'], [51, 'CG']], [[34, 'CB'], [54, 'NE']], 
                [[37, 'NZ'], [46, 'NZ']], [[39, 'NZ'], [46, 'NZ']],
                [[37, 'NZ'], [54, 'NE']], [[36, 'CG'], [54, 'NE']]],
                
             2:[[[32, 'CG'], [47, 'CG']], [[33, 'NZ'], [50, 'NE']]]}
     
    top = traj.topology
    atom_pairs = [[top.select('residue=={0} and name=={1}'.format(*i[0]))[0],
                   top.select('residue=={0} and name=={1}'.format(*i[1]))[0]]
                   for i in pairs[SUMO]] 
                                   
    return md.compute_distances(traj, atom_pairs=atom_pairs, periodic=False)


def rmsf(traj):
    r = traj.superpose(traj[0]).xyz
    r_mean = r.mean(axis=0)
    rmsf_all = (((r - r_mean)**2).sum(axis=2).mean(axis=0))**0.5
    return rmsf_all[traj.topology.select('name==CA')]

def rmsd(traj):
    r = traj.superpose(traj[0]).xyz
    return (((r - r[0])**2).sum(axis=2).mean(axis=1))**0.5


'''
trajectory files:
SUMO1_2asq_noSIM_[123456]_400ns.h5
SUMO1_2asq_wSIM_[123456]_400ns.h5
'''

nstr = [md.load('SUMO1_2asq_noSIM_{}_400ns.h5'.format(i)) for i in arange(1, 7)]
wstr = [md.load('SUMO1_2asq_wSIM_{}_400ns.h5'.format(i)) for i in arange(1, 7)]

#convexhull volume
volns = array([convexhull_volume(i, 1) for i in nstr]).T 
volws = array([convexhull_volume(i, 1) for i in wstr]).T 

figure(1)
plot(volns, 'C1.', markersize=1)
plot(volws, 'C0.', markersize=1, alpha=0.5) 
ylim(1, 2.2)
savefig('SUMO1_2asq_volume.jpg', dpi=600)

figure(2)
hist(volns.ravel(), color='C1', bins=40, linewidth=1, 
     orientation='horizontal')
hist(volws.ravel(), color='C0', alpha=0.6, bins=40, linewidth=1, 
     orientation='horizontal')
legend(['no SIM', 'with SIM'], frameon=0) 
ylim(1, 2.2)
savefig('SUMO1_2asq_volume_hist.jpg', dpi=600)


#F36-Y51
distws = array([md.compute_distances(i, atom_pairs=[[126, 252]], periodic=False).ravel() 
                for i in wstr]).T 
distns = array([md.compute_distances(i, atom_pairs=[[126, 252]], periodic=False).ravel() 
                for i in nstr]).T

figure(1)
plot(distns, 'C1.', markersize=1)
plot(distws, 'C0.', markersize=1, alpha=0.5) 
ylim(0.5, 1.25)
savefig('SUMO1_2asq_dist_F35CG_Y51CG.jpg', dpi=600)

figure(2)
hist(distns.ravel(), color='C1', bins=40, linewidth=1, 
     orientation='horizontal')
hist(distws.ravel(), color='C0', alpha=0.6, bins=40, linewidth=1, 
     orientation='horizontal')
legend(['no SIM', 'with SIM'], frameon=0) 
ylim(0.5, 1.25)
savefig('SUMO1_2asq_dist_F35CG_Y51CG_hist.jpg', dpi=600)




