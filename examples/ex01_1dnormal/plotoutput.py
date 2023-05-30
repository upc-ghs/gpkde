'''
Example 01 Gaussian particle distribution
'''
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm


basedir = os.getcwd()
fname   = 'gpkde.out'
log     = 'gpkde.log'

# load density
df = pd.read_csv( os.path.join( basedir, fname ) ,
            header=None,
            delim_whitespace=True,
            skiprows=0,
            index_col=None
        )
colx = 0
coly = 1
colz = 2
cold = 3
colh = 4
# one-base grid indexes
binx = df[colx].to_numpy().astype(np.int32)

# nparticles and grid
N      = 1e5
dx     = 0.01
xo     = -6
xmax   = 6
xg     = np.arange(xo+0.5*dx,xmax+0.5*dx,dx)
nx     = int(abs(xmax-xo)/dx)

# transfer data to arrays
density    = np.full( (nx,), np.nan )
density[binx-1] = df[ cold ].to_numpy().astype(np.float32)
histogram  = np.full( (nx), np.nan )
histogram[binx-1] = df[ colh ].to_numpy().astype(np.float32)
maxh = np.max(histogram[~np.isnan(histogram)])

# figure 
fig  = plt.figure(figsize=(8,8))
ax   = fig.add_subplot(111)
ax.plot( xg, density/N         , label='GPKDE', zorder=2)
ax.plot( xg, histogram/N       , label='HIST' , zorder=0)
ax.plot( xg, norm.pdf(xg, 0, 1), label='NORM' , zorder=1)
ax.set_xlabel('x')
ax.set_ylabel('rho')
ax.legend()

# print the figure
fig.savefig(os.path.join( basedir, 'ex01_figure.png') )
