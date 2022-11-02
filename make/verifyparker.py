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
v      = 118.46153846153845
#v      = 59.230769230769226  # cm/hr
D      = DoverV*v
L      = 4
vTEND  = 0.5
TEND   = vTEND/v # END OF SOLUTE INJECTION
vTFIN  = 1 
TFIN   = vTFIN/v # END OF SIMULATION
deltax = 0.005
x      = np.arange( 0, L+deltax, deltax) 
delr   = deltax
#vol    = 1.9500000000000004e-07

vol    = 9.750000000000002e-08
mass   = 2.8211805555555564e-11 
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



fig = plt.figure()
ax  = fig.add_subplot(111)


# ANALYTICAL
ax.plot(x/(v*TFIN),function(x,v,TFIN,D) - function(x,v,TFIN-TEND,D))


# LOAD FILE
df = pd.read_csv( os.path.join( os.getcwd(), 'gpkde_density_output_' ),
                     header=None,
                     delim_whitespace=True,
                     skiprows=0,
                     index_col=None
                 )
ax.plot( x[ df[0]-4 ]/(1.01*v*TFIN), df[3].to_numpy()*mass*delr/vol, color='r' )
ax.plot( x[ df[0]-4 ]/(1.01*v*TFIN), df[4].to_numpy()*mass/vol )


plt.show( block=False) 
import pdb
pdb.set_trace()


