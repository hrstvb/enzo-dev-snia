# matplotlib-based plotting script for Iliev et al. test #5
# Daniel R. Reynolds, reynolds@smu.edu

# imports
from pylab import *
import numpy as np

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

# define some helpful functions
def get_params(file):
    """Returns t, vol, gamma, dUnit, tUnit, lUnit from a given parameter file"""
    import shlex
    f = open(file)
    for line in f:
        text = shlex.split(line)
        if ("InitialTime" in text):
            tval = float(text[len(text)-1])
        elif ("DensityUnits" in text):
            dUnit = float(text[len(text)-1])
        elif ("TimeUnits" in text):
            tUnit = float(text[len(text)-1])
        elif ("LengthUnits" in text):
            lUnit = float(text[len(text)-1])
        elif ("Gamma" in text):
            gamma = float(text[len(text)-1])
        elif ("DomainLeftEdge" in text):
            xL = float(text[len(text)-3])
            yL = float(text[len(text)-2])
            zL = float(text[len(text)-1])
        elif ("DomainRightEdge" in text):
            xR = float(text[len(text)-3])
            yR = float(text[len(text)-2])
            zR = float(text[len(text)-1])
    file2 = file + '.rtmodule'
    f = open(file2)
    for line in f:
        text = shlex.split(line)
        if ("E1Units" in text):
            E1Unit = float(text[len(text)-1])
        elif ("E2Units" in text):
            E2Unit = float(text[len(text)-1])
        elif ("E3Units" in text):
            E3Unit = float(text[len(text)-1])
    xL = xL*lUnit
    xR = xR*lUnit
    yL = yL*lUnit
    yR = yR*lUnit
    zL = zL*lUnit
    zR = zR*lUnit
    vol = (xR-xL)*(yR-yL)*(zR-zL)
    tval = tval*tUnit
    return [tval, vol, xR, yR, zR, gamma, dUnit, tUnit, lUnit, E1Unit, E2Unit, E3Unit]

def load_vals(tdump):
    """Returns t, vol, gamma, rho, HI, HII, Ef, E1, E2, E3, tE, vx, vy, vz from a given data dump"""
    import h5py
    import numpy as np
    sdump = repr(tdump).zfill(4)
    pfile = 'DD' + sdump + '/data' + sdump
    hfile = pfile + '.cpu0000'
    tval, vol, xR, yR, zR, gamma, dUnit, tUnit, lUnit, E1Unit, E2Unit, E3Unit = get_params(pfile)
    f = h5py.File(hfile,'r')
    Ef = f.get('/Grid00000001/FS_Radiation_Energy')
    E1 = f.get('/Grid00000001/Radiation_Energy1')
    E2 = f.get('/Grid00000001/Radiation_Energy2')
    E3 = f.get('/Grid00000001/Radiation_Energy3')
    tE = f.get('/Grid00000001/Total_Energy')
    HI = f.get('/Grid00000001/HI_Density')
    HII = f.get('/Grid00000001/HII_Density')
    rho = f.get('/Grid00000001/Density')
    vx = f.get('/Grid00000001/x-velocity')
    vy = f.get('/Grid00000001/y-velocity')
    vz = f.get('/Grid00000001/z-velocity')
    HI = np.multiply(HI,dUnit)
    HII = np.multiply(HII,dUnit)
    rho = np.multiply(rho,dUnit)
    Ef = np.multiply(Ef,dUnit*lUnit*lUnit/tUnit/tUnit)
    E1 = np.multiply(E1,E1Unit*dUnit*lUnit*lUnit/tUnit/tUnit)
    E2 = np.multiply(E2,E2Unit*dUnit*lUnit*lUnit/tUnit/tUnit)
    E3 = np.multiply(E3,E3Unit*dUnit*lUnit*lUnit/tUnit/tUnit)
    tE = np.multiply(tE,lUnit*lUnit/tUnit/tUnit)
    vx = np.multiply(vx,lUnit/tUnit)
    vy = np.multiply(vy,lUnit/tUnit)
    vz = np.multiply(vz,lUnit/tUnit)
    return [tval, vol, gamma, xR, yR, zR, rho, HI, HII, Ef, E1, E2, E3, tE, vx, vy, vz]



# loop over snapshots, loading values and times
for tstep in range(0,te+1):
    
    # load relevant information
    t, vol, gamma, xR, yR, zR, rho, HI, HII, Ef, E1, E2, E3, tE, vx, vy, vz = load_vals(tstep)

    # compute some derived quantities
    xHI  = np.divide(HI, rho)
    xHII = np.divide(HII,rho)
    kE = 0.5*(np.multiply(vx,vx) + np.multiply(vy,vy) + np.multiply(vz,vz))
    Press = (gamma-1.0)*np.multiply(rho,tE-kE)
    Temp = mp/kb*np.divide(Press,2.0*rho - HI)
    iE = np.subtract(tE,kE)
    nDens = np.divide(rho,mp)
    
    # compute volume element
    nx, ny, nz = Ef.shape
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
    #TempIn  = 0.0
    #TempOut = 0.0
    #ncells = 0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                if (sqrt(i*i+j*j+k*k) <= Icells):
                    #PressIn += Te[i][j][k]*rho[i][j][k]*(gamma-1) - 0.5*(
                    #    vx[i][j][k]**2 + vy[i][j][k]**2 + vz[i][j][k]**2) 
                    PressIn += tE[i][j][k]*rho[i][j][k]*(gamma-1.0)
                    DensIn += rho[i][j][k]
                    #TempIn += Te[i][j][k]*rho[i][j][k]*(gamma-1.0)/(2.0*rho[i][j][k] 
                    #                                                - HI[i][j][k])*mp/kb
                    #ncells += 1
                #else:
                #    TempOut += Te[i][j][k]*rho[i][j][k]*(gamma-1.0)/(2.0*rho[i][j][k] 
                #                                                     - HI[i][j][k])*mp/kb
    #Ti = Ti/ncells                            # average temperature behind i-front
    #Te = Te/(nx*ny*nz - ncells)               # average temperature beyond i-front
    rs0 = (3.0*Ngammadot/4.0/pi/aHII/nH/nH)**(1.0/3.0)
    TempIn = 1.0e4                             # temperature inside HII region
    TempOut = 1.0e2                            # temperature outside HII region
    rf  = (2.0*TempIn/TempOut)**(2.0/3.0)*rs0  # predicted final HII region radius
    DensIn = max(DensIn, rho[0][0][0])
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
#    if (tstep > -1):
        
        # set time label
        if (tstep == 0):
            Myr = '0'
        elif (tstep == 1):
            Myr = '10'
        elif (tstep == 3):
            Myr = '30'
        elif (tstep == 10):
            Myr = '100'
        elif (tstep == 20):
            Myr = '200'
        elif (tstep == 50):
            Myr = '500'
        else:
            Myr = '1000'
#        Myr = repr(tstep)
        
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
        
        # Ef slice through z=0
        figure()
        sl = log10(Ef[:][:][0])
        h = imshow(sl, hold=False, extent=(0.0, 1.0, 0.0, 1.0), origin='lower')
        colorbar(h)
        title('log FS radiation density, t =' + Myr + ' Myr')
        savefig('Efcontour_' + Myr + 'Myr' + pictype)
        
        # E1 slice through z=0
        figure()
        sl = log10(E1[:][:][0])
        h = imshow(sl, hold=False, extent=(0.0, 1.0, 0.0, 1.0), origin='lower')
        colorbar(h)
        title('log radiation 1 density, t =' + Myr + ' Myr')
        savefig('E1contour_' + Myr + 'Myr' + pictype)
        
        # E2 slice through z=0
        figure()
        sl = log10(E2[:][:][0])
        h = imshow(sl, hold=False, extent=(0.0, 1.0, 0.0, 1.0), origin='lower')
        colorbar(h)
        title('log radiation 2 density, t =' + Myr + ' Myr')
        savefig('E2contour_' + Myr + 'Myr' + pictype)
        
        # E3 slice through z=0
        figure()
        sl = log10(E3[:][:][0])
        h = imshow(sl, hold=False, extent=(0.0, 1.0, 0.0, 1.0), origin='lower')
        colorbar(h)
        title('log radiation 3 density, t =' + Myr + ' Myr')
        savefig('E3contour_' + Myr + 'Myr' + pictype)
        
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
        Efprof = zeros(Nradii, dtype=float)
        E1prof = zeros(Nradii, dtype=float)
        E2prof = zeros(Nradii, dtype=float)
        E3prof = zeros(Nradii, dtype=float)
        kEprof = zeros(Nradii, dtype=float)
        tEprof = zeros(Nradii, dtype=float)
        iEprof = zeros(Nradii, dtype=float)
        Pprof = zeros(Nradii, dtype=float)
        nHprof = zeros(Nradii, dtype=float)
        rsquared = linspace(1,Nradii,Nradii)
        rsquared = np.multiply(rsquared,rsquared)
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    idx = rad_idx[i][j][k]
                    HIprof[idx]  += xHI[i][j][k]
                    HIIprof[idx] += xHII[i][j][k]
                    Tprof[idx]   += Temp[i][j][k]
                    Efprof[idx]  += Ef[i][j][k]
                    E1prof[idx]  += E1[i][j][k]
                    E2prof[idx]  += E2[i][j][k]
                    E3prof[idx]  += E3[i][j][k]
                    Pprof[idx]   += Press[i][j][k]
                    kEprof[idx]  += kE[i][j][k]
                    tEprof[idx]  += tE[i][j][k]
                    iEprof[idx]  += iE[i][j][k]
                    nHprof[idx]  += nDens[i][j][k]
                    Hcount[idx]  += 1
        HIprof  = log10(HIprof/Hcount)
        HIIprof = log10(HIIprof/Hcount)
        Tprof   = log10(Tprof/Hcount)
        Efprof  = log10(Efprof/Hcount)
        E1prof  = log10(E1prof/Hcount)
        E2prof  = log10(E2prof/Hcount)
        E3prof  = log10(E3prof/Hcount)
        Pprof   = log10(Pprof/Hcount)
        kEprof  = log10(kEprof/Hcount)
        tEprof  = log10(tEprof/Hcount)
        iEprof  = log10(iEprof/Hcount)
        nHprof  = log10(nHprof/Hcount)
        for i in range(Nradii):
            HIprof[i]  = min(max(HIprof[i],-100),100)
            HIIprof[i] = min(max(HIIprof[i],-100),100)
            Tprof[i]   = min(max(Tprof[i],-100),100)
            Efprof[i]  = min(max(Efprof[i],-100),100)
            E1prof[i]  = min(max(E1prof[i],-100),100)
            E2prof[i]  = min(max(E2prof[i],-100),100)
            E3prof[i]  = min(max(E3prof[i],-100),100)
            Pprof[i]   = min(max(Pprof[i],-100),100)
            kEprof[i]  = min(max(kEprof[i],-100),100)
            tEprof[i]  = min(max(tEprof[i],-100),100)
            iEprof[i]  = min(max(iEprof[i],-100),100)
            nHprof[i]  = min(max(nHprof[i],-100),100)
        
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
        
        # FS Radiation profile
        figure()
        plot(Hradii,Efprof)
        grid()
        xlabel('$r/L_{box}$')
        ylabel('log(E_f)')
        title('FS Radiation Profile, t =' + Myr + ' Myr')
        #axis([ 0.0, 1.2, -35.0, -10.0 ])
        savefig('EfProfile_' + Myr + 'Myr' + pictype)
        
        # Radiation 1 profile
        figure()
        plot(Hradii,E1prof)
        grid()
        xlabel('$r/L_{box}$')
        ylabel('log(E_1)')
        title('Radiation 1 Profile, t =' + Myr + ' Myr')
        #axis([ 0.0, 1.2, -35.0, -10.0 ])
        savefig('E1Profile_' + Myr + 'Myr' + pictype)
        
        # Radiation 2 profile
        figure()
        plot(Hradii,E2prof)
        grid()
        xlabel('$r/L_{box}$')
        ylabel('log(E_2)')
        title('Radiation 2 Profile, t =' + Myr + ' Myr')
        #axis([ 0.0, 1.2, -35.0, -10.0 ])
        savefig('E2Profile_' + Myr + 'Myr' + pictype)
        
        # Radiation 3 profile
        figure()
        plot(Hradii,E3prof)
        grid()
        xlabel('$r/L_{box}$')
        ylabel('log(E_3)')
        title('Radiation 3 Profile, t =' + Myr + ' Myr')
        #axis([ 0.0, 1.2, -35.0, -10.0 ])
        savefig('E3Profile_' + Myr + 'Myr' + pictype)
        
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
