#2018-02-14 script to plot dose effect curve of compounds against 3xSu1 and 3xSu2

import pandas as pd
from scipy.optimize import curve_fit
import numpy as np

# concentration points (1:2 or 1:3 series from 500uM)
x_2 = np.array([500 / 2**i for i in arange(8)])[::-1]
x_3 = np.array([500 / 3**i for i in arange(8)])[::-1]

# points used to plot fitted curve
xp = arange(0, 1000, 0.1)

def dose_curve_fit(m, x):
    #F = lambda x, k, IC50: 1 / ( 1 + exp( k * (x - IC50) )) 
    F = lambda x, k, IC50: 1 / ( 1 + ( x / IC50 ) ** (-k) )

    IC50_all = []
    for y in array(m.iloc[:, 1:].T):
        p, q = curve_fit(F, x, y, p0=[-1, 50] )
        IC50_all.append(p[1])
        print(p[1])
        plot(xp, F(xp, *p), '0.8')
        #plot (p[1], 0.5, 'o')

    boxplot(array(m.iloc[:, 1:]).T, labels=x.round().astype(int), positions=x, 
            widths=0.2*x)
    xscale('log'); xlim(0.1, 1000)

    print('mean IC50: ', mean(IC50_all), '+/-', array(IC50_all).std())
    print('n: ', len(IC50_all))

    return IC50_all

# copy data from Excel 
m = pd.read_clipboard(header=None)

# choose the correct concentration series for fitting
dose_curve_fit(m, x_2) 


'''
# mean values fit
y = array(m.iloc[:, 1:].mean(1))
y_std = array(m.iloc[:, 1:].std(1))

p, q = curve_fit(F, x, y, p0=[-1, 50], sigma=y_std)
print('mean IC50: ', p[1])

xp = arange(0, 500, 0.1)
plot(xp, F(xp, p[0], p[1]))
'''





