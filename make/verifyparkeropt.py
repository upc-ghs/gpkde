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


DoverV = 0.01
#v      = 118.46153846153845 # CLASSICAL WITH FLOW BOUNDARIES, SOME VARIABILITY
v      = 118.16967337525458  # WITH CONSTANT HEAD BOUNDARIES, STD < 1e-12
#v      = 59.230769230769226  # cm/hr
D      = DoverV*v
L      = 4
vTEND  = 0.5
TEND   = vTEND/v # END OF SOLUTE INJECTION
vTFIN  = 1 
TFIN   = vTFIN/v # END OF SIMULATION
deltax = 0.005
x      = np.arange( 0, L+deltax, deltax) 
delr   = 1000*deltax
#vol    = 1.9500000000000004e-07

vol    = 9.750000000000002e-08
#vol    = 1000*9.750000000000002e-08
facmas = 1.0005794666076087
mass  = 6.10880084019226e-12 #17
mass  = mass*facmas
#mass = 1.1931251641000507e-11           # 16
#mass = 4.887040672153808e-11           # 15
#mass   = 9.545001312800406e-11
#mass   = 2.2625188297008372e-10
#mass   = 7.636001050240325e-10
#mass   = 9.545001312800406e-11
#mass   = 1.1931251641000507e-11
#mass   = 6.10880084019226e-12
#mass   = 1.1931251641000507e-11
#mass   = 1.190185546875e-11
#mass   = 6.09375e-12
#mass   = 6.09375e-12  
#mass   = 1.190185546875e-11
print('MMASSS ', mass )
#ass   = 2.8211805555555564e-11 
#ass   = 9.521484375e-11
# mass = 6.093750000000001e-10
#mass   = 1.2187500000000002e-09
#mass   = 4.7607421875e-11

#vol    = 9.750000000000002e-08
#mass   = 1.50462962962963e-10

#mass   = 3.00925925925926e-10 # WITH NSUBS 3
#mass   = 1.5234375e-09 # WITH NSUBS 2
#mass   = 3.80859375e-10    # WITH NSUBS 4 
#mass   = 4.7607421875e-11  # WITH NSUBS 8 

co = 1

nparticles = 10 # NOMINAL
#nparticlestotal = 1600000 # REAL
##nparticlestotal = 1599988 # REAL
##nparticlestotal = 1599616 # REAL
##nparticles = 8 # NOMINAL
##nparticlestotal = 819007 # REAL
#total_flow_rate = 0.00231
#mass = co*total_flow_rate*TEND/nparticlestotal
##mass = co*total_flow_rate*dtrelease/(nparticles)
#print('MMASSS ', mass )




fig = plt.figure()
ax  = fig.add_subplot(111)


# ANALYTICAL
ax.plot(x/(v*TFIN),function(x,v,TFIN,D) - function(x,v,TFIN-TEND,D), zorder=0, linewidth=2)


# LOAD FILE
# STATE ZERO
zerodf = pd.read_csv( os.path.join( os.getcwd(), 'gpkde_density_output_' ),
                     header=None,
                     delim_whitespace=True,
                     skiprows=0,
                     index_col=None
                 )
ax.plot( x[ zerodf[0]-4 ]/(1.0*v*TFIN), 1.0*zerodf[3].to_numpy()*mass*delr/vol, zorder=1, linewidth=1, color='r' )
ax.plot( x[ zerodf[0]-4 ]/(1.0*v*TFIN), 1.0*zerodf[4].to_numpy()*mass/vol, color='k', zorder=2, linewidth=0.8, alpha=0.5 )

loopids = [0,1]
for li in loopids:
    # LOAD FILE
    df = pd.read_csv( os.path.join( os.getcwd(), 'gpkde_density_output_'+str(li) ),
                         header=None,
                         delim_whitespace=True,
                         skiprows=0,
                         index_col=None
                     )
    print( np.mean(df[3].to_numpy()*mass*delr/vol ) )

    ax.plot( x[ df[0]-4 ]/(1.0*v*TFIN), 1.0*df[3].to_numpy()*mass*delr/vol , linewidth=0.7, zorder=10-li, label=str(li) )



hlambda = 5 
h       = deltax*hlambda/(v*TFIN)
xo      = 0.71/(v*TFIN)
print( xo )
ax.axvline(xo    , linewidth=0.75, zorder=0, color='k', linestyle='--')
ax.axvline(xo+h  , linewidth=0.75, zorder=0, color='k', linestyle='--')
ax.axvline(xo-h  , linewidth=0.75, zorder=0, color='k', linestyle='--')
ax.axvline(xo+3*h, linewidth=0.75, zorder=0, color='k', linestyle='--')
ax.axvline(xo-3*h, linewidth=0.75, zorder=0, color='k', linestyle='--')



#ax.set_yscale('log')
ax.set_ylim([0.25,1])
#ax.set_xlim([0.5,1])
fig.legend()
plt.savefig( 'figureopt.png' )


plt.show( block=False) 
import pdb
pdb.set_trace()


