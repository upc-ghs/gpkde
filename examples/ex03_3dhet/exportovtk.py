'''
Example 03 heterogeneous 3D distribution
'''
import os
import numpy as np
import pandas as pd
import pyevtk
import argparse

# arguments
parser = argparse.ArgumentParser( description='Control plot options.' )
parser.add_argument( '--fname', type=str, help='load data from the given fname' )
args   = parser.parse_args()

# fname
basedir = os.getcwd()
if args.fname is None:
    fname = 'gpkde.out'
else:
    fname = args.fname

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
binz = df[colz].to_numpy().astype(np.int32)

# grid size
lx, ly, lz = 200, 150, 100
dx, dy, dz = 1, 1, 1
nx, ny, nz = int(lx/dx), int(ly/dy), int(lz/dz)
X,Y,Z      = np.mgrid[0:nx+1,0:ny+1,0:nz+1]
X          = X*dx + 0.5*dx 
Y          = Y*dy + 0.5*dy
Z          = Z*dz + 0.5*dz

# transfer data to arrays
density    = np.full( (nx,ny,nz), np.nan )
density[binx-1,biny-1,binz-1] = df[ cold ].to_numpy().astype(np.float32)
histogram  = np.full( (nx,ny,nz), np.nan )
histogram[binx-1,biny-1,binz-1] = df[ colh ].to_numpy().astype(np.float32)

# write vtk
pyevtk.hl.gridToVTK(
        os.path.join( basedir, 'gpkde' ),
        X,Y,Z,
        cellData={
            'histrho' : histogram,
            'gpkderho': density
        },
    )
