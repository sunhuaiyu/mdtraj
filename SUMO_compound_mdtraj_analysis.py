import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md
from glob import glob
from multiprocessing import Pool

def SUMO_ligand_dist(tr):
    #coordinates for the Cgamma of SUMO1_F36, SUMO2_F31, or SUMO3_F31:
    select_str = '(resname==PHE and (resid==15 or resid==30 or resid==17)) and (name==CG)'
    atom_ix = tr.topology.select(select_str)[0]
    a = tr.xyz[:, atom_ix]

    # ligand all atom coordinatess:
    lig = tr.atom_slice(tr.topology.select('chainid==1'))
    # ligand center of mass: 
    b = md.compute_center_of_mass(lig)

    # distance between K37/K32_CA and ligand center of mass:
    return (((a - b) ** 2).sum(1)) ** 0.5


# read trajectory file in HDF5 format (*.h5), compute SUMO_ligand_dist
def name2traj(file_name):
    tr = md.load(file_name)
    if tr.n_frames > 10000:
            tr = tr[::10]
    return tr


# given trajectory file name in HDF5 format, plot SUMO_ligand_dist
def plot_dist(traj_name):
    plt.plot(SUMO_ligand_dist(name2traj(traj_name)), linewidth=1)
    plt.ylim(0, 4.5)
    title = traj_name.split('.')[0]
    plt.title(title)
    plt.savefig(title + '.jpg', dpi=600)
    plt.close()


# calculate fraction of frames where the distance is less than a cut-off 
compound = ['PHG00686', 'SEW05414', 'HTS12268', 'BTB13496']
compound2traj_name = {i: glob('SUMO1_2uyz_{}_F*_5000ns.h5'.format(i)) for i in compound}
traj_files = sum(list(compound2traj_name.values()))

# traj_dict contains all loaded trajectories 
# dist_dict contains all calculated distances; 
# accelerate calculation with multiprocessing

def D(file_name):
    tr = name2traj(file_name)
    d  = SUMO_ligand_dist(tr)
    return [tr, d]

DD = Pool(48).map(D, traj_files) 
traj_dict = {i[0]:i[1][0] for i in zip(traj_files, DD)}
dist_dict = {i[0]:i[1][1] for i in zip(traj_files, DD)}


# distance (nm) threshold
T = 0.7

# calculate the fraction of trajectories with compound at SIM-binding site
for cp in compound:
    all_dist = np.array([dist_dict[i] for i in compound2traj_name[cp]]).ravel()
    bound_frames, total_frames = sum(all_dist < T), len(all_dist)
    fraction = bound_frames/total_frames
    print(cp, round(fraction, 3), total_frames//1000)
    

# plotting: stack all distance plot together for each compound
for cp in compound:
    n = len(compound2traj_name[cp])
    fig, axs = plt.subplots(nrows=n, ncols=1, sharex=True)    
    fig.set_figheight(n)
    fig.set_figwidth(4) 
    axs[0].set_title(cp)
    
    for i in np.arange(n):
        dc = dist_dict['SUMO1_2uyz_{0}_F{1}_5000ns.h5'.format(cp, i+1)]
        bound = dc < T
        unbound = np.invert(bound)
        length = dc.shape[0]

        axs[i].plot(np.arange(length)[unbound], dc[unbound], 
                    'C1.', markersize=0.5, alpha=0.6)     
        axs[i].plot(np.arange(length)[bound], dc[bound], 
                    'C0.', markersize=0.5, alpha=0.6)

        axs[i].set_ylim(0, 4.5)
    
    fig.subplots_adjust(hspace=0)
    fig.savefig('SUMO1_2uyz_{}_dist_all_traj.jpg'.format(cp),
                dpi=600, bbox_inches='tight')


# extract a centroid frame from each traj ending with significant binding;
# for each compound, superpose all centroids along the SIM-binding pocket
# and save as one pdb file
centroids = {cp:[] for cp in compound}
for cp in compound:
    n = len(compound2traj_name[cp])
    for i in np.arange(n):
        file_name = 'SUMO1_2uyz_{0}_F{1}_5000ns.h5'.format(cp, i+1)
        dc = dist_dict[file_name]
        bound = dc < T
        
        if sum(bound) > 1000:
            tr = traj_dict[file_name][bound]
            protein_atoms = tr.topology.select('residue 32 to 56')
            compound_atoms = tr.topology.select('chainid==1')
            atoms_ix = np.concatenate((protein_atoms, compound_atoms))
            tr.superpose(tr, frame=0, atom_indices=atoms_ix)

            m = np.empty((tr.n_frames, tr.n_frames)) # rmsd matrix
            for i in range(tr.n_frames): 
                m[i] = md.rmsd(tr, tr, i, atom_indices=atoms_ix)

            #compute the centroid frame: the one closest to mean frame
            centroid_ix = np.exp(-m/m.std()).sum(1).argmax()
            centroids[cp].append(tr[centroid_ix])
            
            print(file_name)
            
    centroids_tr = md.join(centroids[cp])
    centroids_tr.superpose(centroids_tr, frame=0, atom_indices=protein_atoms)
    centroids_tr.save_pdb('SUMO1_2uyz_{}_bound_centroids.pdb'.format(cp))

# compute  RMSD among bound_centroids
from scipy.spatial.distance import squareform 
for cp in compound:
    tr = md.load('SUMO1_2uyz_{}_bound_centroids.pdb'.format(cp))
    m = array([md.rmsd(tr, tr, i, atom_indices=protein_atoms) for i in range(len(tr))])
    m = squareform(m, checks=False)
    print(cp, min(m), max(m))
    

# compute atomic distances 
T = 0.7
tr2uyz = md.join([md.load('SUMO1_2uyz_{}_400ns.h5'.format(i+1)) for i in range(12)])

cp = 'PHG00686'
d = [dist_dict['SUMO1_2uyz_{0}_F{1}_5000ns.h5'.format(cp, i+1)] for i in range(12)]

tr1cp = md.join([traj_dict['SUMO1_2uyz_{0}_F{1}_5000ns.h5'.format(cp, i+1)][d[i] < T] for i in range(12)])

def atom_pair_dist3(cp, pair='F36CG_R54CZ'):
    top = tr2uyz[0].topology
    s = pair.split('_')
    pair_ix = top.select_pairs('residue=={0} and name=={1}'.format(s[0][1:3], s[0][3:]),
                               'residue=={0} and name=={1}'.format(s[1][1:3], s[1][3:]))
    
    dist2uyz = md.compute_distances(tr2uyz, atom_pairs=pair_ix, periodic=False) 
                    
    dist1cp = md.compute_distances(tr1cp, atom_pairs=pair_ix, periodic=False) 
                

    fig = plt.figure(figsize=(10, 4.8))

    gs = GridSpec(1, 2, width_ratios=[2, 1]) 
    ax0, ax1 = plt.subplot(gs[0]), plt.subplot(gs[1])
    
    ax0.plot(dist2uyz, 'C1.', markersize=1)
    ax0.plot(dist1cp, 'C0.', markersize=1, alpha=0.5)
    ax0.tick_params(labelsize=15)
    
    ax1.hist(dist2uyz, color='C1', bins=100, linewidth=1, 
         orientation='horizontal')
    ax1.hist(dist1cp, color='C0', alpha=0.6, bins=100, linewidth=1, 
         orientation='horizontal')
    ax1.tick_params(labelsize=15)
    ax1.legend(['no compound', 'with {}'.format(cp)], fontsize=15, frameon=0)     
    
    fig.tight_layout()     
    fig.savefig('SUMO1_2uyz_{0}_dist_{1}.jpg'.format(cp, pair), dpi=600)

