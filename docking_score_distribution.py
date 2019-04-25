
import pandas as pd
from glob import glob

# list of all purchased compounds
g = [i.split('.')[0] for i in glob('*.svg')]


# SUMO1 docking result
m = pd.read_csv('SUMO1_dock_Maybridge_20171031.txt', sep='\t')

plot(m.weighted_mean_dG, '.', color='grey', markersize=1)
ylim(-10, 0)

for i in m.index: 
    if m.iloc[i, 0] in g: 
        plot(i, m.iloc[i, -1], marker='o', color='k')
        
savefig('SUMO1_all.jpg', dpi=600)


#SUMO2 docking result
m2 = pd.read_csv('SUMO2_4npn_1wm2_1wm3_dock_Maybridge_201807.txt', sep='\t')

plot(m2.weighted_mean_dG, '.', color='grey', markersize=1)
ylim(-10, 0)

for i in m2.index: 
    if m2.iloc[i, 0] in g: 
        plot(i, m2.iloc[i, -1], marker='o', color='k')

savefig('SUMO2_all.jpg', dpi=600)


#
x = [np.real(m.loc[m.iloc[:, 0]==i, 'weighted_mean_dG'])[0] for i in g]
y = [np.real(m2.loc[m2.iloc[:, 0]==i, 'weighted_mean_dG'])[0] for i in g]
scatter(x, y)


p1 = m.sort_values(axis=0, by='Unnamed: 0')
p2 = m2.sort_values(axis=0, by='Unnamed: 0')  
scatter(p1['weighted_mean_dG'], p2['weighted_mean_dG'], s=2, marker='.', color='k')
xlim(-10, 0); ylim(-10, 0)
savefig('SUMO1_vs_SUMO2_all.jpg', dpi=600)

