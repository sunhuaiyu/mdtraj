import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md
from glob import glob

def SUMO_ligand_dist(tr):
    #coordinates for the Cgamma of SUMO1_F36, SUMO2_F31, or SUMO3_F31:
    select_str = '(resname==PHE and (resid==15 or resid==30 or resid==17)) and (name==CG)'
    sel = tr.topology.select(select_str)
    a = md.compute_center_of_mass(tr.atom_slice(sel))

    # ligand all atom coordinatess:
    lig = tr.atom_slice(tr.topology.select('chainid==1'))
    # ligand center of mass: 
    b = md.compute_center_of_mass(lig)

    # distance between K37/K32_CA and ligand center of mass:
    return np.sqrt(((a - b) ** 2).sum(1))

# read trajectory in HDF5 format
for traj in glob('*.h5'):
    tr = md.load(traj)
    if tr.n_frames > 10000:
        tr = tr[::10]
    plt.plot(SUMO_ligand_dist(tr), linewidth=1)
    plt.ylim(0, 4.5)
    plt.savefig(traj[:-2] + 'jpg', dpi=600)
    plt.close()



# calculate fraction of frames where the distance is less than a cut-off 
from multiprocessing import Pool

def map_fn(traj_name):
    tr = md.load(traj_name)
    if tr.n_frames > 10000:
            tr = tr[::10]
    dist = SUMO_ligand_dist(tr)   
    return [sum(dist < 0.7), len(dist)]

compound = ['PHG00686', 'SEW05414', 'HTS12268', 'BTB13496']
traj_dict = {i: glob('SUMO1_2uyz_{}_F*_5000ns.h5'.format(i)) for i in compound}
traj_files = sum(list(traj_dict.values()))

d = dict(zip(traj_files, Pool(48).map(map_fn1, traj_files)))

for cp in compound:
    bound_frames, total_frames = np.array(list(map(lambda i: d[i], traj_dict[cp]))).sum(0)
    fraction = bound_frames/total_frames
    print(cp, round(fraction, 3), total_frames//1000)
    

  
