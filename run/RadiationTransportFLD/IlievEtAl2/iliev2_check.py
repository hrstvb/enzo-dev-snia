# matplotlib-based plotting script for Iliev et al. test #2
# Daniel R. Reynolds, reynolds@smu.edu

# imports
from pylab import *
import numpy as np
execfile("../utilities.py")

# set the total number of snapshots
te = 50

# set the solution tolerance
tol = 0.002

# load the reference solution
r_sol = np.load('r_sol.npy')

# set some constants
Ngammadot = 5.0e48     # ionization source strength [photons/sec]
aHII = 2.52e-13        # recombination rate coefficient
mp = 1.67262171e-24    # proton mass [g]
Myr = 3.15576e13       # duration of a Megayear [sec]
nH = 1.0e-3            # input hydrogen number density [cm^(-3)]
trec = 1.0/(aHII*nH)   # recombination time [sec]
rs0 = (3.0*Ngammadot/4/pi/aHII/nH/nH)**(1.0/3.0)   # Stromgren radius


# initialize time-history outputs
#    row 1: time (t)
#    row 2: computed i-front radius
rdata = zeros( (2, te+1), dtype=float);


# loop over snapshots, loading values and times
for tstep in range(0,te+1):
    
    # load relevant information
    t, vol, Eg, xHI, xHII, Temp = load_vals(tstep)
    
    # compute volume element
    nx, ny, nz = Eg.shape
    dV = vol/nx/ny/nz
    
    # compute I-front radius (assuming spherical)
    HIIvolume = sum(xHII)*dV*8.0
    radius = (3.0/4.0*HIIvolume/pi)**(1.0/3.0)
    
    # store data
    rdata[0][tstep] = t/trec
    rdata[1][tstep] = radius
    


# compute I-front radius comparison, error norm
r_err = (rdata[1] - r_sol)/rs0
r_err_norm = (np.sum(np.multiply(r_err,r_err))/te)**(0.5)
if (r_err_norm < tol):
    print 'Error of ',r_err_norm,' is below tolerance ',tol
    print 'PASS'
else:
    print 'Error of ',r_err_norm,' is above tolerance ',tol
    print 'FAIL'
