import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
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
    return (((a - b) ** 2).sum(1))**0.5


# read trajectory file in HDF5 format (*.h5), compute SUMO_ligand_dist
def name2dist(traj_name):
    tr = md.load(traj_name)
    if tr.n_frames > 10000:
            tr = tr[::10]
    return SUMO_ligand_dist(tr)   


# given trajectory file name in HDF5 format, plot SUMO_ligand_dist
def plot_dist(traj_name):
    plt.plot(name2dist(traj_name), linewidth=1)
    plt.ylim(0, 4.5)
    title = traj_name.split('.')[0]
    plt.title(title)
#    plt.savefig(title + 'jpg', dpi=600)
#    plt.close()


# calculate fraction of frames where the distance is less than a cut-off 
compound = ['PHG00686', 'SEW05414', 'HTS12268', 'BTB13496']
traj_dict = {i: glob('SUMO1_2uyz_{}_F*_5000ns.h5'.format(i)) for i in compound}
traj_files = sum(list(traj_dict.values()))

# d contains all distances of all trajectories
d = dict(zip(traj_files, Pool(48).map(name2dist, traj_files)))

# distance (nm) threshold
T = 0.7

# calculate the fraction of trajectories with compound at SIM-binding site
for cp in compound:
    all_dist = np.array([d[i] for i in traj_dict[cp]]).ravel()
    bound_frames, total_frames = sum(all_dist < T), len(all_dist)
    fraction = bound_frames/total_frames
    print(cp, round(fraction, 3), total_frames//1000)
    

# stack all distance plot together for each compound
for cp in compound:
    n = len(traj_dict[cp])
    fig, axs = plt.subplots(nrows=n, ncols=1, sharex=True)    
    fig.set_figheight(n)
    fig.set_figwidth(4) 
    axs[0].set_title(cp)
    
    for i in arange(n):
        dc = d['SUMO1_2uyz_{0}_F{1}_5000ns.h5'.format(cp, i+1)]
        bound = dc < T
        unbound = np.invert(bound)
        length = dc.shape[0]
        axs[i].plot(arange(length)[bound], dc[bound], 'C0.', markersize=0.5)
        axs[i].plot(arange(length)[unbound], dc[unbound], 'C1.', markersize=0.5)
        axs[i].set_ylim(0, 4.5)
    
    fig.subplots_adjust(hspace=0)
    fig.savefig('SUMO1_2uyz_{}_dist_all_traj.jpg'.format(cp),
                dpi=600, bbox_inches='tight')


  
