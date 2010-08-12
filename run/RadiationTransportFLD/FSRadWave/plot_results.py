# matplotlib-based plotting script for plane-wave streaming radiation test
# Daniel R. Reynolds, reynolds@smu.edu

import h5py
from pylab import *

# set the graphics output type
pictype = '.png'

# define some helpful functions
def get_params(file):
    """Returns dUnit, tUnit, lUnit from a given parameter file"""
    import shlex
    f = open(file)
    for line in f:
        text = shlex.split(line)
        if ("DensityUnits" in text):
            dUnit = float(text[len(text)-1])
        elif ("TimeUnits" in text):
            tUnit = float(text[len(text)-1])
        elif ("LengthUnits" in text):
            lUnit = float(text[len(text)-1])
    return [dUnit, tUnit, lUnit]

def load_snapshot(tdump):
    """Returns x, E for a given data output"""
    import h5py
    import numpy as np
    sdump = repr(tdump).zfill(4)
    pfile = 'DD' + sdump + '/data' + sdump
    hfile = pfile + '.cpu0000'
    f = h5py.File(hfile,'r')
    Ef3D = f.get('/Grid00000001/FS_Radiation')
    dUnit, tUnit, lUnit = get_params(pfile)
    rUnit = dUnit*lUnit*lUnit/tUnit/tUnit
    nx, ny, nz = Ef3D.shape
    if nx > ny*nz:
        E = sum(sum(Ef3D,axis=2),axis=1)/ny/nz*rUnit
    elif ny > nx*nz:
        E = sum(sum(Ef3D,axis=2),axis=0)/nx/nz*rUnit
    else:
        E = sum(sum(Ef3D,axis=0),axis=0)/nx/ny*rUnit
    N = E.size
    x = linspace(0.0, 1.0, N)
    return [x, E]
    

# plot first snapshot
figure(1)
[x,E] = load_snapshot(0)
xanal = [0.095, 0.105]
Eanal = [max(E), min(E)]
plot(x,E,'k-',xanal,Eanal,'k:')
xlabel('x')
ylabel('E')
title('Streaming radiation history')
figure(2)
semilogy(x,E,'k-',xanal,Eanal,'k:')
xlabel('x')
ylabel('E')
title('Streaming radiation history')

# repeat for next dataset
[x,E] = load_snapshot(1)
xanal = [0.195, 0.205]
figure(1)
plot(x,E,'b-',xanal,Eanal,'b:')
figure(2)
semilogy(x,E,'b-',xanal,Eanal,'b:')

# repeat for next dataset
[x,E] = load_snapshot(2)
xanal = [0.195, 0.205]
figure(1)
plot(x,E,'g-',xanal,Eanal,'g:')
figure(2)
semilogy(x,E,'g-',xanal,Eanal,'g:')

# repeat for next dataset
[x,E] = load_snapshot(3)
xanal = [0.295, 0.305]
figure(1)
plot(x,E,'r-',xanal,Eanal,'r:')
figure(2)
semilogy(x,E,'r-',xanal,Eanal,'r:')

# repeat for next dataset
[x,E] = load_snapshot(4)
xanal = [0.395, 0.405]
figure(1)
plot(x,E,'c-',xanal,Eanal,'c:')
figure(2)
semilogy(x,E,'c-',xanal,Eanal,'c:')

# repeat for next dataset
[x,E] = load_snapshot(5)
xanal = [0.495, 0.505]
figure(1)
plot(x,E,'m-',xanal,Eanal,'m:')
figure(2)
semilogy(x,E,'m-',xanal,Eanal,'m:')

# repeat for next dataset
[x,E] = load_snapshot(6)
xanal = [0.595, 0.605]
figure(1)
plot(x,E,'y-',xanal,Eanal,'y:')
figure(2)
semilogy(x,E,'y-',xanal,Eanal,'y:')

# repeat for next dataset
[x,E] = load_snapshot(7)
xanal = [0.695, 0.705]
figure(1)
plot(x,E,'k-',xanal,Eanal,'k:')
figure(2)
semilogy(x,E,'k-',xanal,Eanal,'k:')

# repeat for next dataset
[x,E] = load_snapshot(8)
xanal = [0.795, 0.805]
figure(1)
plot(x,E,'b--',xanal,Eanal,'b:')
figure(2)
semilogy(x,E,'b--',xanal,Eanal,'b:')

# repeat for next dataset
[x,E] = load_snapshot(9)
xanal = [0.895, 0.905]
figure(1)
plot(x,E,'g--',xanal,Eanal,'g:')
figure(2)
semilogy(x,E,'g--',xanal,Eanal,'g:')

# repeat for next dataset
[x,E] = load_snapshot(10)
xanal = [0.995, 1.0]
figure(1)
plot(x,E,'r--',xanal,Eanal,'r:')
figure(2)
semilogy(x,E,'r--',xanal,Eanal,'r:')

# finish off plot and save to file
figure(1)
#axis([0.0, 1.0, -30.0, -14.0])
savefig('rad_snapshots' + pictype)
figure(2)
#axis([0.0, 1.0, -60, -15])
savefig('rad_snapshots_log' + pictype)
