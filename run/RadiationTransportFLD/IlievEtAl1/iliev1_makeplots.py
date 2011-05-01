# matplotlib-based plotting script for Iliev et al. test #1
# Daniel R. Reynolds, reynolds@smu.edu

# imports
from pylab import *
execfile("../utilities.py")

# set the total number of snapshots
te = 50

# set the graphics output type
pictype = '.png'

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
#    row 3: predicted i-front radius (rI)
#    row 4: stromgren sphere radius (rs)
rdata = zeros( (4, te+1), dtype=float);


# loop over snapshots, loading values and times
for tstep in range(te+1):
    
    # load relevant information
    t, vol, Eg, xHI, xHII, Temp = load_vals(tstep)
    
    # compute volume element
    nx, ny, nz = Eg.shape
    dV = vol/nx/ny/nz
    
    # compute I-front radius (assuming spherical)
    HIIvolume = sum(xHII)*dV*8.0
    radius = (3.0/4.0*HIIvolume/pi)**(1.0/3.0)
    
    # compute analytical solution
    ranal = rs0*(1.0 - exp(-t/trec))**(1.0/3.0)
    
    # store data
    rdata[0][tstep] = t/trec
    rdata[1][tstep] = radius
    rdata[2][tstep] = ranal
    rdata[3][tstep] = rs0
    
    # generate 2D plots at certain times
    if (tstep == 1) or (tstep == 10) or (tstep == 50):
        
        # set time label
        if (tstep == 1):
            Myr = '10'
        elif (tstep == 10):
            Myr = '100'
        else:
            Myr = '500'
        
        # set mesh
        x = linspace(0.0,1.0,nx)
        y = linspace(0.0,1.0,ny)
        X, Y = meshgrid(x,y)
        
        # xHI slice through z=0
        figure()
        sl = log10(xHI[:][:][0])
        h = imshow(sl, hold=False, extent=(0.0, 1.0, 0.0, 1.0), origin='lower')
        colorbar(h)
        title('log HI fraction, t =' + Myr + ' Myr')
        savefig('HIcontour_' + Myr + 'Myr' + pictype)
        
        # Eg slice through z=0
        figure()
        sl = log10(Eg[:][:][0])
        h = imshow(sl, hold=False, extent=(0.0, 1.0, 0.0, 1.0), origin='lower')
        colorbar(h)
        title('log radiation density, t =' + Myr + ' Myr')
        savefig('Econtour_' + Myr + 'Myr' + pictype)
        
        # spherically-averaged profiles for xHI, xHII
        Nradii = nx*3/2
        Hradii = linspace(0.0,sqrt(3.0),Nradii)
        rad_idx = zeros( (nx,ny,nz), dtype=float)
        for k in range(nz):
            zloc = (k+0.5)/nz
            for j in range(ny):
                yloc = (j+0.5)/ny
                for i in range(nx):
                    xloc = (i+0.5)/nx
                    rad_idx[i][j][k] = max(0,floor(sqrt(xloc*xloc + yloc*yloc + zloc*zloc)/sqrt(3.0)*Nradii))
        Hcount = 1.0e-16*ones(Nradii)
        HIprof = zeros(Nradii, dtype=float)
        HIIprof = zeros(Nradii, dtype=float)
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    idx = rad_idx[i][j][k]
                    HIprof[idx] += xHI[i][j][k]
                    HIIprof[idx] += xHII[i][j][k]
                    Hcount[idx] += 1
        HIprof = log10(HIprof/Hcount)
        HIIprof = log10(HIIprof/Hcount)
        figure()
        plot(Hradii,HIprof,'b-',Hradii,HIIprof,'r--')
        grid()
        xlabel('$r/L_{box}$')
        ylabel('log(xHI), log(xHII)')
        title('HI, HII Profiles, t =' + Myr + ' Myr')
        legend( ('xHI','xHII') )
        axis([ 0.0, 1.2, -7.0, 1.0 ])
        savefig('profiles_' + Myr + 'Myr' + pictype)


# I-front radius/velocity plots
figure()
plot(rdata[0],rdata[1]/rs0,'b-',rdata[0],rdata[2]/rs0,'r--')
xlabel('$t/t_{rec}$')
ylabel('$r_I/r_S$')
title('Propagation of HII Region')
legend( ('computed', 'analytical'), loc=4 )
grid()
axis([ 0.0, 4.5, 0.0, 1.1 ])
savefig('rad_vs_time' + pictype)

figure()
times = (rdata[0][1:te] + rdata[0][2:te+1])/2.0
velocity = (rdata[1][2:te+1] - rdata[1][1:te])/(rdata[0][2:te+1]-rdata[0][1:te])/rs0
vel_anal = (rdata[2][2:te+1] - rdata[2][1:te])/(rdata[0][2:te+1]-rdata[0][1:te])/rs0
plot(times,velocity,'b-',times,vel_anal,'r--')
xlabel('$t/t_{rec}$')
ylabel('$v_I/(r_S/t_{rec})$')
title('Velocity of HII Region')
legend( ('computed', 'analytical') )
grid()
axis([ 0.0, 4.5, -0.2, 1.2 ])
savefig('vel_vs_time' + pictype)
