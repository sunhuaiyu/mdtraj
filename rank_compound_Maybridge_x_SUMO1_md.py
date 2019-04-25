# Extract free energy score from docking results and rank the compounds.
# 
from glob import glob
import pandas as pd
from os import system
import numpy as np

result = pd.DataFrame()

# frames and respective frame_counts are NOT in chronological order wrt MD trajectory
# calculation refer to rmsd_based_clustering_SUMO3_2mp2.py
frames = ['frame543', 'frame3300', 'frame98', 'frame1809', 'frame3727']
frame_counts = [1334, 1531, 118, 566, 453] 

for frame in frames:
    G = dict()
    for name in glob('vina_' + frame + '/maybridge*_out.pdbqt'):
        f = open(name, 'rt')
        f.readline()
        G[name.split('/')[-1].split('_')[0]] = float(f.readline().split()[3])
 
    result[frame] = pd.Series(G)
    # result.to_csv('SUMO3_' + frame + '_docking_rank.txt', sep='\t')

result['weighted_mean_dG'] = np.average(result, axis=1, weights=frame_counts).round(2)
ranked_result = result.sort_values(by='weighted_mean_dG')

ranked_result.to_csv('SUMO1_dock_Maybridge_20171031.txt', sep='\t')


terms = ['MOLFORMULA', 'MOLWEIGHT', 'product_name', 'c_log_p',
         'h_bond_donors', 'h_bond_acceptors', 'ACD_Code', 'code']


#ranked_result = pd.read_csv('SUMO1_dock_Maybridge_20171031.txt', sep='\t', index_col=0)
#collect top250 docking result

system('mkdir ./top250_per_frame3300')
top250 = ranked_result[:250].copy()
for i in top250.index:
    top250.loc[i, 'SMILES'] = open('/media/data/hsun/Maybridge/{}.smi'.format(i), 'rt').readline().rstrip()
    system('cp /media/data/hsun/Maybridge/{}.sdf ./top250_per_frame3300/'.format(i))
    system('cp ./vina_frame3300/{}_out.pdbqt ./top250_per_frame3300/'.format(i))
    system('babel ./vina_frame3300/{0}_out.pdbqt -o sdf ./top250_per_frame3300/{0}_out.sdf'.format(i))
    system('babel ./top250_per_frame3300/{0}.sdf -o svg ./top250_per_frame3300/{0}.svg'.format(i))

    with open('./top250_per_frame3300/{}.sdf'.format(i)) as sdf:
        line = sdf.readline()
        while len(line) > 0:
            if line[0] == '>':
                name = line.split('<')[1].split('>')[0]
                if name in terms:
                    top250.loc[i, name] = sdf.readline().rstrip()
            line = sdf.readline()

top250['MOLWEIGHT'] = top250['MOLWEIGHT'].astype(float).round(2)
top250['c_log_p'] = top250['c_log_p'].astype(float).round(2)
top250.loc[:, ['h_bond_donors', 'h_bond_acceptors']] =\
    top250.loc[:, ['h_bond_donors', 'h_bond_acceptors']].astype(float).astype(int)

top250.to_csv('SUMO1_dock_Maybridge_20171031_top250.txt', sep='\t')

