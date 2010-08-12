# matplotlib-based plotting script for point-source streaming radiation test
# Daniel R. Reynolds, reynolds@smu.edu

import h5py
from pylab import *
import numpy as np

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
    E = Ef3D[:][0][0]
    N = E.size
    x = linspace(0.01, 1.0, N)
    return [x, E]

def load_contour(tdump):
    """Returns E for a given data output"""
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
    E = Ef3D[:][:][0]
    return [E]
    

# plot first snapshot
figure(1)
[x,E] = load_snapshot(0)
plot(x,E,'k-')
xlabel('x')
ylabel('E')
title('Streaming radiation history')
figure(2)
semilogy(x,E,'k-')
xlabel('x')
ylabel('E')
title('Streaming radiation history (log)')
figure(3)
Escale = np.multiply(np.multiply(E,x),x)*4.0*pi
plot(x,Escale,'k-')
xlabel('x')
ylabel('$4 \pi r^2 E$')
title('Scaled streaming radiation history')
figure(4)
semilogy(x,Escale,'k-')
xlabel('x')
ylabel('$4 \pi r^2 E$')
title('Scaled streaming radiation history (log)')

# repeat for next dataset
[x,E] = load_snapshot(1)
figure(1)
plot(x,E,'b-')
figure(2)
semilogy(x,E,'b-')
Escale = np.multiply(np.multiply(E,x),x)*4.0*pi
figure(3)
plot(x,Escale,'b-')
figure(4)
semilogy(x,Escale,'b-')

# repeat for next dataset
[x,E] = load_snapshot(2)
figure(1)
plot(x,E,'g-')
figure(2)
semilogy(x,E,'g-')
Escale = np.multiply(np.multiply(E,x),x)*4.0*pi
figure(3)
plot(x,Escale,'g-')
figure(4)
semilogy(x,Escale,'g-')

# repeat for next dataset
[x,E] = load_snapshot(3)
figure(1)
plot(x,E,'r-')
figure(2)
semilogy(x,E,'r-')
Escale = np.multiply(np.multiply(E,x),x)*4.0*pi
figure(3)
plot(x,Escale,'r-')
figure(4)
semilogy(x,Escale,'r-')

# repeat for next dataset
[x,E] = load_snapshot(4)
figure(1)
plot(x,E,'c-')
figure(2)
semilogy(x,E,'c-')
Escale = np.multiply(np.multiply(E,x),x)*4.0*pi
figure(3)
plot(x,Escale,'c-')
figure(4)
semilogy(x,Escale,'c-')

# repeat for next dataset
[x,E] = load_snapshot(5)
figure(1)
plot(x,E,'m-')
figure(2)
semilogy(x,E,'m-')
Escale = np.multiply(np.multiply(E,x),x)*4.0*pi
figure(3)
plot(x,Escale,'m-')
figure(4)
semilogy(x,Escale,'m-')

# repeat for next dataset
[x,E] = load_snapshot(6)
figure(1)
plot(x,E,'y-')
figure(2)
semilogy(x,E,'y-')
Escale = np.multiply(np.multiply(E,x),x)*4.0*pi
figure(3)
plot(x,Escale,'y-')
figure(4)
semilogy(x,Escale,'y-')

# repeat for next dataset
[x,E] = load_snapshot(7)
figure(1)
plot(x,E,'k-')
figure(2)
semilogy(x,E,'k-')
Escale = np.multiply(np.multiply(E,x),x)*4.0*pi
figure(3)
plot(x,Escale,'k-')
figure(4)
semilogy(x,Escale,'k-')

# repeat for next dataset
[x,E] = load_snapshot(8)
figure(1)
plot(x,E,'b--')
figure(2)
semilogy(x,E,'b--')
Escale = np.multiply(np.multiply(E,x),x)*4.0*pi
figure(3)
plot(x,Escale,'b--')
figure(4)
semilogy(x,Escale,'b--')

# repeat for next dataset
[x,E] = load_snapshot(9)
figure(1)
plot(x,E,'g--')
figure(2)
semilogy(x,E,'g--')
Escale = np.multiply(np.multiply(E,x),x)*4.0*pi
figure(3)
plot(x,Escale,'g--')
figure(4)
semilogy(x,Escale,'g--')

# repeat for next dataset
[x,E] = load_snapshot(10)
figure(1)
plot(x,E,'r--')
figure(2)
semilogy(x,E,'r--')
Escale = np.multiply(np.multiply(E,x),x)*4.0*pi
figure(3)
plot(x,Escale,'r--')
figure(4)
semilogy(x,Escale,'r--')

# finish off plot and save to file
figure(1)
#axis([0.0, 1.0, -30.0, -14.0])
savefig('rad_snapshots' + pictype)
figure(2)
#axis([0.0, 1.0, 0.0, -15.0])
savefig('rad_snapshots_log' + pictype)
figure(3)
savefig('scaled_rad_snapshots' + pictype)
figure(4)
savefig('scaled_rad_snapshots_log' + pictype)


# plot contour at halfway point
[E] = load_contour(5)
nx, ny = E.shape
x = linspace(0.0,1.0,nx)
y = linspace(0.0,1.0,ny)
X, Y = meshgrid(x,y)
figure(5)
h = imshow(E, hold=False, extent=(0.0, 1.0, 0.0, 1.0), origin='lower')
colorbar(h)
title('E contour at T/2')
savefig('contour' + pictype)

figure(6)
h = imshow(log10(E), hold=False, extent=(0.0, 1.0, 0.0, 1.0), origin='lower')
colorbar(h)
title('log(E) contour at T/2')
savefig('contour_log' + pictype)


