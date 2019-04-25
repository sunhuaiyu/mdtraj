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



