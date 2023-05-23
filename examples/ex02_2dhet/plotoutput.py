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
binx = df[colx].to_numpy().astype(np.int32)
biny = df[coly].to_numpy().astype(np.int32)
binz = df[colz].to_numpy().astype(np.int32)

# grid size
nx, ny, nz = 400, 200, 1

# transfer data to arrays
density    = np.full( (nx,ny,nz), np.nan )
density[binx-1,biny-1,binz-1] = df[ cold ].to_numpy().astype(np.float32)
histogram  = np.full( (nx,ny,nz), np.nan )
histogram[binx-1,biny-1,binz-1] = df[ colh ].to_numpy().astype(np.float32)
maxh = np.max(histogram[~np.isnan(histogram)])


# figure 
fig  = plt.figure(figsize=(16,12))
norm = colors.Normalize(vmin=0, vmax=1.0)
# gpkde
ax   = fig.add_subplot(211)
im   = ax.imshow(density[:,:,0].transpose()/maxh, origin='lower', cmap=cm.turbo)
ax.set_aspect('equal')
im.set_norm(norm)
fig.colorbar(im,ax=ax, label='rho')
ax.set_xlabel('x[m]')
ax.set_ylabel('y[m]')
ax.text(0.8,0.8,'gpkde',transform=ax.transAxes)
# histogram
ax   = fig.add_subplot(212)
im   = ax.imshow(histogram[:,:,0].transpose()/maxh, origin='lower', cmap=cm.turbo)
ax.set_aspect('equal')
im.set_norm(norm)
fig.colorbar(im,ax=ax, label='rho')
ax.set_xlabel('x[m]')
ax.set_ylabel('y[m]')
ax.text(0.8,0.8,'histogram',transform=ax.transAxes)

fig.savefig(os.path.join( basedir, 'ex02_figure.png') )
