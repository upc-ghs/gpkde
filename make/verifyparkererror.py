import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd 
import scipy.special as sci
import os 

# PARKER
def function( x, v, t, D, co=1 ):
    first = 0.5*sci.erfc((x-v*t)/(2*np.sqrt(D*t))) 
    sec   = np.sqrt(v**2*t/(np.pi*D))*np.exp(-(x-v*t)**2/(4*D*t))
    third = 0.5*(1+v*x/D+v**2*t/D)*np.exp(v*x/D)*sci.erfc((x+v*t)/(2*np.sqrt(D*t)))
    return first + sec - third


DoverV = 0.1
#DoverV = 0.01
v      = 118.16967337525458  # WITH CONSTANT HEAD BOUNDARIES, STD < 1e-12
D      = DoverV*v
L      = 4
vTEND  = 0.5
TEND   = vTEND/v # END OF SOLUTE INJECTION
vTFIN  = 1 
TFIN   = vTFIN/v # END OF SIMULATION
deltax = 0.005
x      = np.arange( 0, L+deltax, deltax) 
delr   = deltax

vol    = 9.750000000000002e-08
co     = 1




fig = plt.figure(figsize=(10,10))

xlims =[0,2]

ax0  = fig.add_subplot(311)
ax1  = fig.add_subplot(312)
ax2  = fig.add_subplot(313)
#ax3  = fig.add_subplot(614)
#ax4  = fig.add_subplot(615)
#ax5  = fig.add_subplot(616)


# ANALYTICAL
#ax0.plot(x/(v*TFIN),function(x,v,TFIN,D) - function(x,v,TFIN-TEND,D), zorder=0, linewidth=2, label='analytical', color='k', alpha=0.5)

# LOAD FILE
# STATE ZERO
# DV 01
fname = 'gpkde_particles_dev.csv.np1.024.dv01'
#fname = 'gpkde_particles_dev.csv.np0.128.dv01'

# DV001
#fname = 'gpkde_particles_dev.csv.np0.0128'
#fname = 'gpkde_particles_dev.csv.np0.0016'
#fname = 'gpkde_particles_dev.csv.np0.8192'
#fname = 'gpkde_particles_dev.csv.np0.1024'
#fname = 'gpkde_particles_dev.csv.np1.6'
fname = fname + str(999)
zerodf = pd.read_csv( os.path.join( os.getcwd(), fname  ),
                     header=None,
                     delim_whitespace=True,
                     skiprows=0,
                     index_col=None
                 
                     )
ax0.plot( zerodf[ 0 ], zerodf[1].to_numpy(), zorder=1, linewidth=1, color='g', label='reldensity' )
ax0.plot( zerodf[ 0 ], zerodf[3].to_numpy(), color='r', zorder=2, linewidth=0.8, alpha=0.5, label='relroughness' )
ax0.plot( zerodf[ 0 ], zerodf[5].to_numpy(), color='b', zorder=2, linewidth=0.8, alpha=0.5, label='relsmoothing' )

ax0.set_xlabel('relerrorloop')
ax0.set_xlabel('nloop')
ax0.set_yscale('log')
ax0.set_ylim([1e-5,10])
ax0.axhline(1e-2)
ax0.axhline(5e-2, color='k')

ax1.plot( zerodf[ 0 ], zerodf[7].to_numpy(), zorder=1, linewidth=1, color='g', label='rmsehist' )
ax1.set_yscale('log')
ax1.set_xlabel('nloop')

ax2.plot( zerodf[ 0 ], zerodf[8].to_numpy(), color='r', zorder=2, linewidth=0.8, alpha=0.5, label='almiseproxy' )
ax2.set_yscale('log')
ax2.set_xlabel('nloop')







#colx           = 0
#coly           = 1
#colz           = 2
#coldensity     = 3
#colsmoothingx   = 4
#colsmoothingy   = 5
#colsmoothingz   = 6
#colsmoothingshapex = 7
#colsmoothingshapey = 8
#colsmoothingshapez = 9
#colsmoothingscale  = 10
#
#colroughnessxx   = 16
#colroughnessyy   = 17
#colroughnesszz   = 18
#colroughness     = 19
#
#cellsize = 0.005
#
#
#nloops  = 7
#loopids = np.arange(0,nloops+1,1).astype(np.int32)
#colors  =  ['b','c', 'lime', 'r', 'pink', 'y', 'gold', 'tab:purple', 'tab:orange', 'tab:brown'  ]
#for li in loopids:
#    # LOAD FILE
#    df = pd.read_csv( os.path.join( os.getcwd(), fname +str(li) ),
#                         header=None,
#                         delim_whitespace=True,
#                         skiprows=0,
#                         index_col=None
#                     )
#    ax0.plot( x[ df[0]-4 ]/(1.0*v*TFIN), 1.0*df[coldensity].to_numpy()*mass*delr/vol , linewidth=0.7, zorder=10+li, label=str(li), color=colors[li] )
#
#    # SMOOTHING
#    ax1.plot( x[ df[0]-4 ]/(1.0*v*TFIN), df[colsmoothingx].to_numpy()/cellsize , linewidth=0.7, zorder=10+li , color=colors[li])
#
#    # ROUGHNESS
#    ax2.plot( x[ df[0]-4 ]/(1.0*v*TFIN), df[colroughness].to_numpy() , linewidth=0.7, zorder=10+li, color=colors[li] )
#
#
#tx = 0.5
#ty = 0.5
#
#ax0.text(tx,ty,'density',transform=ax0.transAxes )
#ax1.text(tx,ty,'smoothing',transform=ax1.transAxes )
#ax2.text(tx,ty,'roughness',transform=ax2.transAxes )
##ax2.text(tx,ty,'smoothing',transform=ax2.transAxes )
##ax3.text(tx,ty,'smoothing',transform=ax3.transAxes )
##ax4.text(tx,ty,'smoothing',transform=ax4.transAxes )
##ax5.text(tx,ty,'smoothing',transform=ax5.transAxes )
##ax6.text(tx,ty,'smoothing',transform=ax6.transAxes )
#
#ax0.set_xlim(xlims)
#ax1.set_xlim(xlims)
#ax2.set_xlim(xlims)
##ax2.set_xlim(xlims)
##ax3.set_xlim(xlims)
##ax4.set_xlim(xlims)
##ax5.set_xlim(xlims)
#
#
##ax2.set_yscale('log')
#
##hlambda = 5 
##h       = deltax*hlambda/(v*TFIN)
##xo      = 0.71/(v*TFIN)
##print( xo )
##ax.axvline(xo    , linewidth=0.75, zorder=0, color='k', linestyle='--')
##ax.axvline(xo+h  , linewidth=0.75, zorder=0, color='k', linestyle='--')
##ax.axvline(xo-h  , linewidth=0.75, zorder=0, color='k', linestyle='--')
##ax.axvline(xo+3*h, linewidth=0.75, zorder=0, color='k', linestyle='--')
##ax.axvline(xo-3*h, linewidth=0.75, zorder=0, color='k', linestyle='--')



#ax.set_yscale('log')
#ax.set_ylim([0.25,1])
#ax.set_xlim([0.5,1])
fig.legend()
plt.savefig( 'figureopterror.png' )


plt.show( block=False) 
import pdb
pdb.set_trace()


