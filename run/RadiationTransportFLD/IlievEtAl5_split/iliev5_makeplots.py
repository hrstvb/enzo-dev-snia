# matplotlib-based plotting script for Iliev et al. test #5
# Daniel R. Reynolds, reynolds@smu.edu

# imports
from pylab import *
import numpy as np
execfile("../utilities.py")

# set the total number of snapshots and graphics output type
te = 50
pictype = '.png'

# set some constants
Ngammadot = 5.0e48     # ionization source strength [photons/sec]
aHII = 2.52e-13        # recombination rate coefficient
mp = 1.67262171e-24    # proton mass [g]
kb = 1.3806504e-16     # Boltzmann constant
Myr = 3.15576e13       # duration of a Megayear [sec]
nH = 1.0e-3            # input hydrogen number density [cm^(-3)]
trec = 1.0/(aHII*nH)   # recombination time [sec]
tS = 0

# initialize time-history outputs
#    row 1: time (t)
#    row 2: computed i-front radius
rdata = zeros( (5, te+1), dtype=float);
HIfraction = zeros( (2, te+1), dtype=float);

# overlay the profile plots, so generate storage arrays here
t, vol, gamma, xR, yR, zR, rho, HI, HII, Eg, tE, vx, vy, vz = load_vals(0)
nx, ny, nz = Eg.shape
Nradii = nx*3/2
T_multiprof   = zeros((Nradii, 5), dtype=float)
Eg_multiprof  = zeros((Nradii, 5), dtype=float)
P_multiprof   = zeros((Nradii, 5), dtype=float)
nH_multiprof  = zeros((Nradii, 5), dtype=float)
    

# loop over snapshots, loading values and times
for tstep in range(0,te+1):
    
    # load relevant information
    t, vol, gamma, xR, yR, zR, rho, xHI, xHII, Eg, tE, vx, vy, vz = load_vals_hydro(tstep)

    # compute some derived quantities
    kE = 0.5*(np.multiply(vx,vx) + np.multiply(vy,vy) + np.multiply(vz,vz))
    Press = (gamma-1.0)*np.multiply(rho,tE-kE)
    Temp = mp/kb*np.divide(Press,2.0*rho - HI)
    iE = np.subtract(tE,kE)
    nDens = np.divide(rho,mp)
    
    # compute volume element
    nx, ny, nz = Eg.shape
    dV = vol/nx/ny/nz
    
    # compute I-front radius (assuming spherical)
    HIIvolume = sum(xHII)*dV*8.0
    radius = (3.0/4.0*HIIvolume/pi)**(1.0/3.0)
    
    # store data
    HIvolume = sum(xHI)*dV
    HIfraction[0][tstep] = t/trec
    HIfraction[1][tstep] = HIvolume/vol
    
    # compute Stromgren radius, predicted HII region radius
    #   pressure/density/temperature behind I-front
    Icells = radius/xR*nx       # number of cells behind front
    PressIn = 0.0
    DensIn  = 0.0
    TempIn  = 0.0
    TempOut = 0.0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                if (sqrt(1.0*i*i+1.0*j*j+1.0*k*k) <= Icells):
                    PressIn += tE[i][j][k]*rho[i][j][k]*(gamma-1.0)
                    DensIn += rho[i][j][k]
    rs0 = (3.0*Ngammadot/4.0/pi/aHII/nH/nH)**(1.0/3.0)
    TempIn  = 1.0e4                            # temperature inside HII region
    TempOut = 1.0e2                            # temperature outside HII region
    rf = rs0*(2.0*TempIn/TempOut)**(2.0/3.0)   # predicted final HII region radius
    DensIn = max(DensIn, 1.0e-50)              # stability when no ionized region
    cspeed = sqrt(PressIn/DensIn)              # average sound speed behind front


    # if radius <= rs0, use static solution, otherwise use dynamic (with 't'
    # set as the time since reaching the Stromgren radius.
    # (in both of these, accomodate for the fact that Matlab screws up
    #  fractional powers of negative numbers)
    if (radius < rs0):     # use static solution
        arg = 1.0-exp(-t/trec)
        if (arg < 0.0):
            ranal = -rs0*abs(arg)**(1.0/3.0)
        else:
            ranal = rs0*arg**(1.0/3.0)
    else:                  # use dynamic solution
       if (tS == 0):
           tS = t          # store time at which rs0 is reached
       tdiff = t - tS
       arg = 1.0 + 7.0*cspeed*tdiff/4.0/rs0
       if (arg < 0): 
	  ranal = -rs0*abs(arg)**(4.0/7.0)
       else:
          ranal = rs0*arg**(4.0/7.0)

    # fill radius/time array 
    rdata[0][tstep] = t/trec
    rdata[1][tstep] = radius
    rdata[2][tstep] = ranal
    rdata[3][tstep] = rs0
    rdata[4][tstep] = rf

    # generate 2D plots at certain times
    if ((tstep == 0) or (tstep == 1) or (tstep == 3) or (tstep == 10) or (tstep == 20) or (tstep == 50) or (tstep == 100)):
        
        # set time label
        if (tstep == 0):
            Myr = '0'
            mpentry = -1
        elif (tstep == 1):
            Myr = '10'
            mpentry = 0
        elif (tstep == 3):
            Myr = '30'
            mpentry = 1
        elif (tstep == 10):
            Myr = '100'
            mpentry = 2
        elif (tstep == 20):
            Myr = '200'
            mpentry = 3
        elif (tstep == 50):
            Myr = '500'
            mpentry = 4
        else:
            Myr = '1000'
            mpentry = -1
        
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
        
        # Pressure slice through z=0
        figure()
        sl = log10(Press[:][:][0])
        h = imshow(sl, hold=False, extent=(0.0, 1.0, 0.0, 1.0), origin='lower')
        colorbar(h)
        title('log Pressure, t =' + Myr + ' Myr')
        savefig('PressContour_' + Myr + 'Myr' + pictype)
        
        # Temp slice through z=0
        figure()
        sl = log10(Temp[:][:][0])
        h = imshow(sl, hold=False, extent=(0.0, 1.0, 0.0, 1.0), origin='lower')
        colorbar(h)
        title('log Temperature, t =' + Myr + ' Myr')
        savefig('TempContour_' + Myr + 'Myr' + pictype)
        
        # spherically-averaged profiles for xHI, xHII, Temp, Eg, Press, kE, iE, tE, nH
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
        Tprof = zeros(Nradii, dtype=float)
        Egprof = zeros(Nradii, dtype=float)
        kEprof = zeros(Nradii, dtype=float)
        tEprof = zeros(Nradii, dtype=float)
        iEprof = zeros(Nradii, dtype=float)
        Pprof = zeros(Nradii, dtype=float)
        nHprof = zeros(Nradii, dtype=float)
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    idx = rad_idx[i][j][k]
                    HIprof[idx]  += xHI[i][j][k]
                    HIIprof[idx] += xHII[i][j][k]
                    Tprof[idx]   += Temp[i][j][k]
                    Egprof[idx]  += Eg[i][j][k]
                    Pprof[idx]   += Press[i][j][k]
                    kEprof[idx]  += kE[i][j][k]
                    tEprof[idx]  += tE[i][j][k]
                    iEprof[idx]  += iE[i][j][k]
                    nHprof[idx]  += nDens[i][j][k]
                    Hcount[idx]  += 1
        HIprof  = log10(HIprof/Hcount)
        HIIprof = log10(HIIprof/Hcount)
        Tprof   = log10(Tprof/Hcount)
        Egprof  = log10(Egprof/Hcount)
        Pprof   = log10(Pprof/Hcount)
        kEprof  = log10(kEprof/Hcount)
        tEprof  = log10(tEprof/Hcount)
        iEprof  = log10(iEprof/Hcount)
        nHprof  = log10(nHprof/Hcount)
        for i in range(Nradii):
            HIprof[i]  = min(max(HIprof[i],-100),100)
            HIIprof[i] = min(max(HIIprof[i],-100),100)
            Tprof[i]   = min(max(Tprof[i],-100),100)
            Egprof[i]  = min(max(Egprof[i],-100),100)
            Pprof[i]   = min(max(Pprof[i],-100),100)
            kEprof[i]  = min(max(kEprof[i],-100),100)
            tEprof[i]  = min(max(tEprof[i],-100),100)
            iEprof[i]  = min(max(iEprof[i],-100),100)
            nHprof[i]  = min(max(nHprof[i],-100),100)

        # store profiles for later use
        if (mpentry != -1):
            T_multiprof[:,mpentry]   = Tprof
            Eg_multiprof[:,mpentry]  = Egprof
            P_multiprof[:,mpentry]   = iEprof
            nH_multiprof[:,mpentry]  = nHprof
        
        # chemistry profiles
        figure()
        plot(Hradii,HIprof,'b-',Hradii,HIIprof,'r--')
        grid()
        xlabel('$r/L_{box}$')
        ylabel('log(xHI), log(xHII)')
        title('HI, HII Profiles, t =' + Myr + ' Myr')
        legend( ('xHI','xHII') )
        axis([ 0.0, 1.2, -7.0, 1.0 ])
        savefig('profiles_' + Myr + 'Myr' + pictype)
        
        # Temperature profile
        figure()
        plot(Hradii,Tprof)
        grid()
        xlabel('$r/L_{box}$')
        ylabel('log(T) [K]')
        title('Temperature Profile, t =' + Myr + ' Myr')
        axis([ 0.0, 1.2, 2.0, 5.0 ])
        savefig('TempProfile_' + Myr + 'Myr' + pictype)
        
        # Radiation profile
        figure()
        plot(Hradii,Egprof)
        grid()
        xlabel('$r/L_{box}$')
        ylabel('log(E)')
        title('Radiation Profile, t =' + Myr + ' Myr')
        axis([ 0.0, 1.2, -35.0, -10.0 ])
        savefig('EgProfile_' + Myr + 'Myr' + pictype)
        
        # Pressure profile
        figure()
        plot(Hradii,Pprof)
        grid()
        xlabel('$r/L_{box}$')
        ylabel('log(P)')
        title('Pressure Profile, t =' + Myr + ' Myr')
        savefig('PressProfile_' + Myr + 'Myr' + pictype)
        
        # Energy profiles
        figure()
        plot(Hradii,tEprof,'b-',Hradii,kEprof,'r--',Hradii,iEprof,'k:')
        grid()
        xlabel('$r/L_{box}$')
        ylabel('log(energy)')
        title('Energy Profiles, t =' + Myr + ' Myr')
        legend( ('Total','Kinetic','Internal') )
        axis([ 0.0, 1.2, 7.0, 14.0 ])
        savefig('EnergyProfile_' + Myr + 'Myr' + pictype)
        
        # number density profile
        figure()
        plot(Hradii,nHprof)
        grid()
        xlabel('$r/L_{box}$')
        ylabel('log(n)')
        title('Number Density Profile, t =' + Myr + ' Myr')
        axis([ 0.0, 1.2, -4.0, -2.0 ])
        savefig('nProfile_' + Myr + 'Myr' + pictype)
        

# I-front radius
figure()
plot(rdata[0],rdata[1]/rs0,'b-',rdata[0],rdata[2]/rs0,'r--')
xlabel('$t/t_{rec}$')
ylabel('$r_I/r_S$')
legend( ('Computed','Analytical') )
title('Propagation of HII Region')
axis([ 0.0, 4.5, 0.0, 3.0 ])
grid()
savefig('rad_vs_time' + pictype)

# Neutral fraction history
figure()
plot(HIfraction[0],HIfraction[1],'b-')
xlabel('$t/t_{rec}$')
ylabel('1-x')
title('Total Neutral Fraction')
axis([ 0.0, 4.5, 0.5, 1.1 ])
grid()
savefig('frac_vs_time' + pictype)

# overlaid Temperature profiles
figure()
plot(Hradii,T_multiprof[:,0], Hradii,T_multiprof[:,1], Hradii,T_multiprof[:,2],
     Hradii,T_multiprof[:,3], Hradii,T_multiprof[:,4])
grid()
xlabel('$r/L_{box}$')
ylabel('log(T) [K]')
title('Temperature Profiles')
legend( ('10 Myr','30 Myr','100 Myr','200 Myr','500 Myr') )
axis([ 0.0, 1.2, 2.0, 5.0 ])
savefig('TempProfiles' + pictype)

# overlaid Radiation profiles
figure()
plot(Hradii,Eg_multiprof[:,0], Hradii,Eg_multiprof[:,1], Hradii,Eg_multiprof[:,2],
     Hradii,Eg_multiprof[:,3], Hradii,Eg_multiprof[:,4])
grid()
xlabel('$r/L_{box}$')
ylabel('log(E)')
title('Radiation Profiles')
legend( ('10 Myr','30 Myr','100 Myr','200 Myr','500 Myr') )
axis([ 0.0, 1.2, -35.0, -10.0 ])
savefig('EgProfiles' + pictype)

# overlaid Pressure profiles
figure()
plot(Hradii,P_multiprof[:,0], Hradii,P_multiprof[:,1], Hradii,P_multiprof[:,2], 
     Hradii,P_multiprof[:,3], Hradii,P_multiprof[:,4])
grid()
xlabel('$r/L_{box}$')
ylabel('log(P)')
title('Pressure Profiles')
legend( ('10 Myr','30 Myr','100 Myr','200 Myr','500 Myr') )
savefig('PressProfiles' + pictype)

# overlaid number density profiles
figure()
plot(Hradii,nH_multiprof[:,0], Hradii,nH_multiprof[:,1], Hradii,nH_multiprof[:,2],
     Hradii,nH_multiprof[:,3], Hradii,nH_multiprof[:,4])
grid()
xlabel('$r/L_{box}$')
ylabel('log(n)')
title('Number Density Profiles')
legend( ('10 Myr','30 Myr','100 Myr','200 Myr','500 Myr') )
axis([ 0.0, 1.2, -4.0, -2.0 ])
savefig('nProfiles' + pictype)
