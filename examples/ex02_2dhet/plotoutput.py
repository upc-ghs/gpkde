'''
Example 02 heterogeneous particle distribution
'''
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm

basedir = os.getcwd()
fname   = 'gpkde.out'

# load density
colx = 0
coly = 1
colz = 2
cold = 3
colh = 4
df   = pd.read_csv( os.path.join( basedir, fname ) ,
            header=None,
            delim_whitespace=True,
            skiprows=0,
            index_col=None
        )
# grid indexes in output file are one-based
binx = df[colx].to_numpy().astype(np.int32)
biny = df[coly].to_numpy().astype(np.int32)

# grid size
lx, ly = 400, 200
dx, dy = 1, 1
nx, ny = int(lx/dx), int(ly/dy)
X,Y    = np.mgrid[0:nx,0:ny]
X      = X*dx + 0.5*dx 
Y      = Y*dy + 0.5*dy

# transfer data to arrays
density    = np.full( (nx,ny), np.nan )
density[binx-1,biny-1] = df[ cold ].to_numpy().astype(np.float32)
histogram  = np.full( (nx,ny), np.nan )
histogram[binx-1,biny-1] = df[ colh ].to_numpy().astype(np.float32)
maxh = np.max(histogram[~np.isnan(histogram)])

# figure 
fig  = plt.figure(figsize=(16,12))
norm = colors.Normalize(vmin=0, vmax=1.0)
# gpkde
ax   = fig.add_subplot(211)
im   = ax.pcolormesh(X,Y,density/maxh, cmap=cm.turbo, shading='nearest')
ax.set_aspect('equal')
im.set_norm(norm)
fig.colorbar(im,ax=ax, label='rho')
ax.set_xlabel('x[m]')
ax.set_ylabel('y[m]')
ax.text(0.8,0.8,'gpkde',transform=ax.transAxes)
# histogram
ax   = fig.add_subplot(212)
im   = ax.pcolormesh(X,Y,histogram/maxh, cmap=cm.turbo, shading='nearest')
ax.set_aspect('equal')
im.set_norm(norm)
fig.colorbar(im,ax=ax, label='rho')
ax.set_xlabel('x[m]')
ax.set_ylabel('y[m]')
ax.text(0.8,0.8,'histogram',transform=ax.transAxes)

# print the figure
fig.savefig(os.path.join( basedir, 'ex02_figure.png') )
