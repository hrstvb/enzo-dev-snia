!=======================================================================
!
! Copyright 2009 Daniel R. Reynolds
!
! This software is released under the terms of the "Enzo Public License"
! in the accompanying LICENSE file.
!
!=======================================================================
subroutine MFProb_LocResid(res_ec, res_HI, res_HeI, res_HeII, src_ec,   &
     src_HI, src_HeI, src_HeII, vx, vy, vz, rho, ec, eh, HI, HI0, HeI,  &
     HeI0, HeII, HeII0, piHI, piHeI, piHeII, GHI, GHeI, GHeII, dt,      &
     theta, DualEnergy, a, a0, adot, adot0, gamma, HFrac, Model, CompA, &
     Comp_xray, Comp_temp, NTempBins, TempStart, TempEnd, k1Tb, k2Tb,   &
     k3Tb, k4Tb, k5Tb, k6Tb, ceHITb, ceHeITb, ceHeIITb, ciHITb,         &
     ciHeITb, ciHeISTb, ciHeIITb, reHIITb, reHeII1Tb, reHeII2Tb,        &
     reHeIIITb, bremTb, aUn, dUn, dUn0, vUn, lUn, lUn0, rUn, rUn0, eUn, &
     nUn, nUn0, ecScale, Nchem, dx, dy, dz, Nx, Ny, Nz, NGxl, NGxr,     &
     NGyl, NGyr, NGzl, NGzr, ier)
  !=======================================================================
  !  written by: Daniel R. Reynolds
  !  date:       March, 2009
  !
  !  PURPOSE: Computes the theta-method residuals for the gas energy and 
  !           chemistry equations.
  !
  !  INPUTS:
  !     src_ec     - source function values for gas energy correction eq.
  !     src_HI     - source function values for HI eq.
  !     src_HeI    - source function values for HeI eq.
  !     src_HeII   - source function values for HeII eq.
  !     v*         - velocity arrays in each direction
  !     rho        - density array
  !     ec         - fluid energy correction array (new time)
  !     eh         - total fluid energy array (old time)
  !     HI         - Hydrogen I number density array
  !     HI0        - Hydrogen I number density array (old time)
  !     HeI        - Helium I number density array
  !     HeI0       - Helium I number density array (old time)
  !     HeII       - Helium II number density array
  !     HeII0      - Helium II number density array (old time)
  !     piHI       - Hydrogen I photo-ionization rates
  !     piHeI      - Helium I photo-ionization rates
  !     piHeII     - Helium II photo-ionization rates
  !     GHI        - Hydrogen I photo-heating coefficients
  !     GHeI       - Helium I photo-heating coefficients
  !     GHeII      - Helium II photo-heating coefficients
  !     dt         - time step size
  !     theta      - parameter in time stepping method
  !     DualEnergy - 1 if eh is internal energy, 0 if it is total energy
  !     a          - cosmological expansion parameter
  !     a0         - cosmological expansion parameter (old time)
  !     adot       - da/dt
  !     adot0      - da/dt (old time)
  !     gamma      - constant in ideal gas law
  !     HFrac      - percentage of mass composed of Hydrogen
  !     Model      - flag denoting physical model to use
  !     CompA      - Compton cooling coefficient 1 (multiplier)
  !     Comp_xray  - X-ray Compton heating coefficient
  !     Comp_temp  - X-ray Compton heating temperature 
  !     *Temp*     - bins for temperature-dependent rates
  !     *Tb        - temperature-dependent rate tables
  !     *Un,*Un0   - variable scaling constants (new and old times)
  !     ecScale    - for model 4, this holds the constant temperature
  !     Nchem      - number of chemical species in simulation
  !     dx,dy,dz   - mesh spacing (comoving) in each direction
  !     Nx,Ny,Nz   - active mesh size in each direction
  !     NG*l/NG*r  - left/right ghost cells in each direction
  !
  !     Note: the vector inputs are of size (Nx + NGxl + NGxr) in 
  !     the x-direction, others are similar.
  !
  !  OUTPUT ARGUMENTS: 
  !     res_ec     - nonlinear residual for the fluid energy equation
  !     res_HI     - nonlinear residual for the HI chemistry equation
  !     res_HeI    - nonlinear residual for the HeI chemistry eq.
  !     res_HeII   - nonlinear residual for the HeII chemistry eq.
  !     ier        - success/failure flag (0->failure, 1->success)
  !
  !  EXTERNALS: 
  !
  !  LOCALS:
  !
  !=======================================================================
#include "fortran.def"
  implicit none
    
  !--------------
  ! argument declarations
  integer, intent(in)  :: DualEnergy, Model, Nchem, NTempBins
  integer, intent(in)  :: Nx, NGxl, NGxr
  integer, intent(in)  :: Ny, NGyl, NGyr
  integer, intent(in)  :: Nz, NGzl, NGzr
  integer, intent(out) :: ier
  REALSUB, intent(in)  :: dt, theta, a, a0, adot, adot0
  real, intent(in) :: dx, dy, dz, gamma, HFrac
  real, intent(in) :: CompA, Comp_xray, Comp_temp
  real, intent(in) :: TempStart, TempEnd
  real, intent(in) :: aUn, dUn, dUn0, vUn, lUn, lUn0, rUn, rUn0, &
       eUn, nUn, nUn0, ecScale
  real, intent(in),                                             &
       dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr) ::   &
       src_ec, src_HI, src_HeI, src_HeII, vx, vy, vz, rho, ec, HI,  &
       HI0, HeI, HeI0, HeII, HeII0, eh, piHI, piHeI, piHeII, GHI, GHeI, GHeII
  real, intent(in), dimension(NTempBins) :: k1Tb, k2Tb, k3Tb, k4Tb,  &
       k5Tb, k6Tb, ceHITb, ceHeITb, ceHeIITb, ciHITb, ciHeITb, ciHeISTb, &
       ciHeIITb, reHIITb, reHeII1Tb, reHeII2Tb, reHeIIITb, bremTb
  real, intent(out),                                          &
       dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr) :: &
       res_ec, res_HI, res_HeI, res_HeII

!--------------
! locals
  integer :: i, j, k, Tidx, Tidxp
  real :: thetam, lTempS, lTempE, dlTemp, lTemp, Tl, Tr, Tfac, KEconst
  real :: k1, k2, k3, k4, k5, k6, eint
  real :: k1o, k2o, k3o, k4o, k5o, k6o
  real :: mol_weight, aval, afac, c, hp, kb, mp, zr, gam_1
  real :: dxi2, dyi2, dzi2, DivV, GradRhoDotV
  real :: rhoval, ecval, ehval, nH, nHI, nHII, nHe, nHeI, nHeII, nHeIII, ne
  real :: rhoval0, nH0, nHI0, nHII0, nHe0, nHeI0, nHeII0, nHeIII0, ne0
  real :: min_temp, T, T0, lamT, lamT0, G, G0, Lambda, Lambda0
  real :: ceHI, ceHeI, ceHeII, ciHI, ciHeI, ciHeIS, ciHeII
  real :: ceHI0, ceHeI0, ceHeII0, ciHI0, ciHeI0, ciHeIS0, ciHeII0
  real :: reHII, reHeII1, reHeII2, reHeIII, brem, Comp1, Comp2
  real :: reHII0, reHeII10, reHeII20, reHeIII0, brem0
  
!=======================================================================

  ! initialize outputs to have all zero values, flag to success
  ier = 1
  res_ec = 0.d0
  res_HI = 0.d0
  if (Nchem == 3) then
     res_HeI  = 0.d0
     res_HeII = 0.d0
  endif
  
  ! check that chemistry constants make sense
  if ((Nchem /= 1) .and. (Nchem /= 3)) then
     write(*,*) 'Chemistry ERROR: illegal value, Nchem = ',Nchem, &
          ',  Nchem must be one of {1, 3}.  Returning!'
     ier = 0
     return
  endif
  if ((HFrac < 0.d0) .or. (HFrac > 1.d0)) then
     write(*,*) 'Chemistry ERROR: illegal value, HFrac = ',HFrac, &
          ',  value must be in the interval [0,1].  Returning!'
     ier = 0
     return     
  endif


!!$  print *, 'Entering MFProb_LocResid routine:'
!!$  print *, '    DualEnergy =',DualEnergy
!!$  print *, '    Model =',Model
!!$  print *, '    Nchem =',Nchem
!!$  print *, '    NTempBins =',NTempBins
!!$  print *, '    Nx*',Nx,NGxl,NGxr
!!$  print *, '    Ny*',Ny,NGyl,NGyr
!!$  print *, '    Nz*',Nz,NGzl,NGzr
!!$  print *, '    dt,theta =',dt,theta
!!$  print *, '    a*',a,a0,adot,adot0
!!$  print *, '    dx*',dx,dy,dz
!!$  print *, '    gamma =',gamma
!!$  print *, '    HFrac =',HFrac
!!$  print *, '    Comp* =',CompA,Comp_xray,Comp_temp
!!$  print *, '    Temp* =',TempStart,TempEnd
!!$  print *, '    *Un =',aUn,dUn,dUn0,vUn,lUn,lUn0,rUn,rUn0,eUn,nUn,nUn0,ecScale
!!$  print *, '    src* =',sum(src_ec),sum(src_HI),sum(src_HeI),sum(src_HeII)
!!$  print *, '    v* =',sum(vx),sum(vy),sum(vz)
!!$  print *, '    rho =',sum(rho)
!!$  print *, '    ec,eh =',sum(ec),sum(eh)
!!$  print *, '    H* =',sum(HI),sum(HI0),sum(HeI),sum(HeI0),sum(HeII),sum(HeII0)
!!$  print *, '    pi* =',sum(piHI),sum(piHeI),sum(piHeII)
!!$  print *, '    G* =',sum(GHI),sum(GHeI),sum(GHeII)
!!$  print *, '    k*Tb =',sum(k1Tb),sum(k2Tb),sum(k3Tb),sum(k4Tb),sum(k5Tb),sum(k6Tb)
!!$  print *, '    ce* =',sum(ceHITb),sum(ceHeITb),sum(ceHeIITb)
!!$  print *, '    ci* =',sum(ciHITb),sum(ciHeITb),sum(ciHeISTb),sum(ciHeIITb)
!!$  print *, '    re* =',sum(reHIITb),sum(reHeII1Tb),sum(reHeII2Tb),sum(reHeIIITb)
!!$  print *, '    brem =',sum(bremTb)



  ! initialize constants
  thetam = 1.d0-theta
  aval = a*aUn
  gam_1 = gamma-1.d0
  c  = 2.99792458d10           ! speed of light [cm/s]
  hp = 6.6260693d-27           ! Planck's constant [ergs*s]
  kb = 1.3806504d-16           ! boltzmann constant [erg/K]
  mp = 1.67262171d-24          ! Mass of a proton [g]
  zr = 1.d0/(aval) - 1.d0      ! cosmological redshift
  afac  = adot/a               ! adot/a
  min_temp = 1.d-1             ! minimum temperature [K]
  KEconst = 0.5d0              ! conversion to internal energy
  if (DualEnergy == 1)  KEconst = 0.d0

  !   lookup table constants
  lTempS = log(TempStart)
  lTempE = log(TempEnd)
  dlTemp = (lTempE - lTempS)/(1.d0*NTempBins - 1.d0)

  ! compute shortcuts
  dxi2 = 0.5d0*a/(dx*lUn)  ! convert to proper units during division
  dyi2 = 0.5d0*a/(dy*lUn)
  dzi2 = 0.5d0*a/(dz*lUn)
  if (aval .ne. 1.d0) then        ! Compton cooling coefficients
     Comp1 = CompA*(1.d0 + zr)**4
     Comp2 = 2.73d0*(1.d0 + zr)
  else
     Comp1 = 0.d0
     Comp2 = 0.d0
  endif


  ! compute right-hand sides depending on physical model

  !  First compute gas energy correction-based adjustments alone
  !  (as these involve derivatives, use loop for appropriate dimension)
!!$  ! 1D model
!!$  if (Ny == 1) then
!!$     do k=1,Nz,1
!!$        do j=1,Ny,1
!!$           do i=1,Nx,1
!!$
!!$              ! shortcuts
!!$              ecval = ec(i,j,k) * eUn
!!$              
!!$              ! velocity divergence
!!$              DivV = dxi2*(vx(i+1,j,k)-vx(i-1,j,k))*vUn
!!$              
!!$              ! (grad density).dot.(velocity)/density
!!$              GradRhoDotV = dxi2*(rho(i+1,j,k)-rho(i-1,j,k)) &
!!$                   * vx(i,j,k) / rho(i,j,k) * vUn
!!$              
!!$              ! put it together
!!$              res_ec(i,j,k) = ecval/aval*(DivV - gam_1*GradRhoDotV)
!!$              
!!$           enddo
!!$        enddo
!!$     enddo
!!$
!!$  ! 2D model
!!$  else if (Nz == 1) then
!!$     do k=1,Nz,1
!!$        do j=1,Ny,1
!!$           do i=1,Nx,1
!!$              
!!$              ! shortcuts
!!$              ecval = ec(i,j,k) * eUn
!!$              
!!$              ! velocity divergence
!!$              DivV = (dxi2*(vx(i+1,j,k)-vx(i-1,j,k))  &
!!$                    + dyi2*(vy(i,j+1,k)-vy(i,j-1,k))) * vUn
!!$              
!!$              ! (grad density).dot.(velocity)/density
!!$              GradRhoDotV = (dxi2*(rho(i+1,j,k)-rho(i-1,j,k))*vx(i,j,k)  &
!!$                           + dyi2*(rho(i,j+1,k)-rho(i,j-1,k))*vy(i,j,k)) &
!!$                           / rho(i,j,k)*vUn
!!$              
!!$              ! put it together
!!$              res_ec(i,j,k) = ecval/aval*(DivV - gam_1*GradRhoDotV)
!!$              
!!$           enddo
!!$        enddo
!!$     enddo
!!$     
!!$  ! 3D model
!!$  else
!!$     do k=1,Nz,1
!!$        do j=1,Ny,1
!!$           do i=1,Nx,1
!!$              
!!$              ! shortcuts
!!$              ecval = ec(i,j,k) * eUn
!!$              
!!$              ! velocity divergence
!!$              DivV = (dxi2*(vx(i+1,j,k)-vx(i-1,j,k))  &
!!$                    + dyi2*(vy(i,j+1,k)-vy(i,j-1,k))  &
!!$                    + dzi2*(vz(i,j,k+1)-vz(i,j,k-1))) * vUn
!!$              
!!$              ! (grad density).dot.(velocity)/density
!!$              GradRhoDotV = (dxi2*(rho(i+1,j,k)-rho(i-1,j,k))*vx(i,j,k)  &
!!$                           + dyi2*(rho(i,j+1,k)-rho(i,j-1,k))*vy(i,j,k)  &
!!$                           + dzi2*(rho(i,j,k+1)-rho(i,j,k-1))*vz(i,j,k)) &
!!$                          / rho(i,j,k) * vUn
!!$              
!!$              ! put it together
!!$              res_ec(i,j,k) = ecval/aval*(DivV - gam_1*GradRhoDotV)
!!$              
!!$           enddo
!!$        enddo
!!$     enddo
!!$  endif



  !==================================================
  !  Now incorporate Model-specific equations

  !    isothermal case-B Hydrogen ionization test case
  if (Model == 4) then
     
     ! set molecular weight for temperature computations
     mol_weight = 0.6d0

     ! iterate over the domain
     do k=1,Nz,1
        do j=1,Ny,1
           do i=1,Nx,1
              
              ! set shortcut values for this spatial location,
              ! converting densities from comoving to proper, and 
              ! put all shortcuts in CGS units
              rhoval  = rho(i,j,k)*dUn
              rhoval0 = rho(i,j,k)*dUn0
              nHI  = HI(i,j,k)*nUn
              nHI0 = HI0(i,j,k)*nUn0
              nH  = Hfrac*rhoval/mp
              nH0 = Hfrac*rhoval0/mp
              nHII  = max(1.d0*(nH - nHI), 0.d0)
              nHII0 = max(1.d0*(nH0 - nHI0), 0.d0)
              ne  = nHII
              ne0 = nHII0

              ! temperature shortcuts
              !*** For Model 4 with cosmology, the temperature is held in ecScale ***!
              if (adot == 0.d0) then
                 eint = vUn*vUn*(eh(i,j,k) &
                      - KEconst*(vx(i,j,k)**2 + vy(i,j,k)**2 + vz(i,j,k)**2))
                 T  = gam_1*mol_weight*mp*(ec(i,j,k)*eUn + eint)/kb
                 T  = max(1.d0*T,1.d0*min_temp)
                 T0 = gam_1*mol_weight*mp*eint/kb
                 T0 = max(1.d0*T0,1.d0*min_temp)
              else
                 T  = ecScale
                 T0 = ecScale
              endif
              lamT  = 3.15614d5/T
              lamT0 = 3.15614d5/T0

              ! look up rates
              lTemp = min(max(log(T), lTempS), lTempE)
              Tidx = min(NTempBins-1, max(1, int((lTemp-lTempS)/dlTemp)+1))
              Tidxp = Tidx+1
              Tl = lTempS + (Tidx-1)*dlTemp
              Tr = lTempS +  Tidx*dlTemp
              Tfac = (lTemp - Tl)/(Tr - Tl)
              k1 = k1Tb(Tidx) + (k1Tb(Tidxp) - k1Tb(Tidx))*Tfac

              lTemp = min(max(log(T0), lTempS), lTempE)
              Tidx = min(NTempBins-1, max(1, int((lTemp-lTempS)/dlTemp)+1))
              Tidxp = Tidx+1
              Tl = lTempS + (Tidx-1)*dlTemp
              Tr = lTempS +  Tidx*dlTemp
              Tfac = (lTemp - Tl)/(Tr - Tl)
              k1o = k1Tb(Tidx) + (k1Tb(Tidxp) - k1Tb(Tidx))*Tfac

              ! compute case B Hydrogen recombination coefficient 
              ! [Hui & Gnedin, 1997: RI^B_{HII}]
              k2 = 2.753d-14*lamT**(1.5d0) *                  &
                   (1.d0+(lamT/2.74d0)**(0.407d0))**(-2.242d0)
              k2o = 2.753d-14*lamT0**(1.5d0) *                &
                   (1.d0+(lamT0/2.74d0)**(0.407d0))**(-2.242d0)

              ! compute (comoving, scaled) residual for HI species
              res_HI(i,j,k) = HI(i,j,k) - HI0(i,j,k)          &
                   - theta*dt*(src_HI(i,j,k) + k2*ne*nHII     &
                      - nHI*(k1*ne + piHI(i,j,k)))/nUn        &
                   - thetam*dt*(src_HI(i,j,k) + k2o*ne0*nHII0 &
                      - nHI0*(k1o*ne0 + piHI(i,j,k)))/nUn0
              
              ! decouple gas energy from all other variables
              res_ec(i,j,k) = 0.d0
              
           enddo
        enddo
     enddo

  !==================================================
  !     case B HII recombination rate
  else if (Model == 1) then

     ! Hydrogen only case
     if (Nchem == 1) then
        ! iterate over the domain
        do k=1,Nz,1
           do j=1,Ny,1
              do i=1,Nx,1
                 
                 ! set shortcut values for this spatial location,
                 ! converting densities from comoving to proper, and 
                 ! put all shortcuts in ACTUAL units
                 ecval = ec(i,j,k)*eUn
                 rhoval  = rho(i,j,k)*dUn
                 rhoval0 = rho(i,j,k)*dUn0
                 nHI  = HI(i,j,k)*nUn
                 nHI0 = HI0(i,j,k)*nUn0
                 nH  = Hfrac*rhoval/mp
                 nH0 = Hfrac*rhoval0/mp
                 nHII  = max(1.d0*(nH - nHI), 0.d0)
                 nHII0 = max(1.d0*(nH0 - nHI0), 0.d0)
                 ne  = nHII
                 ne0 = nHII0
                 
                 ! compute temperatures
                 !!!! This assumes ALL density is H, otherwise mol_weight it wrong !!!!
                 eint = vUn*vUn*(eh(i,j,k) &
                      - KEconst*(vx(i,j,k)**2 + vy(i,j,k)**2 + vz(i,j,k)**2))
                 mol_weight = rhoval/mp/(nHI + nHII + ne)
                 T = gam_1*mol_weight*mp*(ec(i,j,k)*eUn + eint)/kb
                 T = max(1.d0*T,1.d0*min_temp)
                 lamT = 3.15614d5/T
                 mol_weight = rhoval0/mp/(nHI0 + nHII0 + ne0)
                 T0 = gam_1*mol_weight*mp*eint/kb
                 T0 = max(1.d0*T0,1.d0*min_temp)
                 lamT0 = 3.15614d5/T0

                 ! look up rates
                 lTemp = min(max(log(T), lTempS), lTempE)
                 Tidx = min(NTempBins-1, max(1, int((lTemp-lTempS)/dlTemp)+1))
                 Tidxp = Tidx+1
                 Tl = lTempS + (Tidx-1)*dlTemp
                 Tr = lTempS +  Tidx*dlTemp
                 Tfac = (lTemp - Tl)/(Tr - Tl)
                 k1 = k1Tb(Tidx) + (k1Tb(Tidxp) - k1Tb(Tidx))*Tfac
                 ceHI = ceHITb(Tidx) + (ceHITb(Tidxp) - ceHITb(Tidx))*Tfac
                 ciHI = ciHITb(Tidx) + (ciHITb(Tidxp) - ciHITb(Tidx))*Tfac
                 reHII = reHIITb(Tidx) + (reHIITb(Tidxp) - reHIITb(Tidx))*Tfac
                 brem = bremTb(Tidx) + (bremTb(Tidxp) - bremTb(Tidx))*Tfac

                 lTemp = min(max(log(T0), lTempS), lTempE)
                 Tidx = min(NTempBins-1, max(1, int((lTemp-lTempS)/dlTemp)+1))
                 Tidxp = Tidx+1
                 Tl = lTempS + (Tidx-1)*dlTemp
                 Tr = lTempS +  Tidx*dlTemp
                 Tfac = (lTemp - Tl)/(Tr - Tl)
                 k1o = k1Tb(Tidx) + (k1Tb(Tidxp) - k1Tb(Tidx))*Tfac
                 ceHI0 = ceHITb(Tidx) + (ceHITb(Tidxp) - ceHITb(Tidx))*Tfac
                 ciHI0 = ciHITb(Tidx) + (ciHITb(Tidxp) - ciHITb(Tidx))*Tfac
                 reHII0 = reHIITb(Tidx) + (reHIITb(Tidxp) - reHIITb(Tidx))*Tfac
                 brem0 = bremTb(Tidx) + (bremTb(Tidxp) - bremTb(Tidx))*Tfac

                 ! compute case B Hydrogen recombination coefficient 
                 ! [Hui & Gnedin, 1997: RI^B_{HII}]
                 k2 = 2.753d-14*lamT**(1.5d0) *                  &
                      (1.d0+(lamT/2.74d0)**(0.407d0))**(-2.242d0)
                 k2o = 2.753d-14*lamT0**(1.5d0) *                &
                      (1.d0+(lamT0/2.74d0)**(0.407d0))**(-2.242d0)
                 
                 ! compute fluid cooling rate.  Terms (in order):
                 !    Collisional Excitations
                 !    Collisional Ionizations
                 !    Recombinations
                 !    Compton cooling or heating 
                 !    X-ray Compton heating
                 !    Bremsstrahlung
                 Lambda = ne/rhoval               &
                      *(ceHI*nHI                  &
                      + ciHI*nHI                  &
                      + reHII*nHII                &
                      + Comp1*(T-Comp2)           &
                      + Comp_xray*(T-Comp_temp)   &
                      + brem*nHII                 &
                      )
                 Lambda0 = ne0/rhoval0            &
                      *(ceHI0*nHI0                &
                      + ciHI0*nHI0                &
                      + reHII0*nHII0              &
                      + Comp1*(T0-Comp2)          &
                      + Comp_xray*(T0-Comp_temp)  &
                      + brem0*nHII0               &
                      )
                 
                 ! compute fluid heating rates
                 G  = nHI/rhoval*GHI(i,j,k)
                 G0 = nHI0/rhoval0*GHI(i,j,k)
                 
                 ! put it all together
                 res_ec(i,j,k) = ec(i,j,k)                       &
                      - theta*dt*(src_ec(i,j,k) + res_ec(i,j,k)  &
                         + G - Lambda - 2.d0*afac*ecval)/eUn     &
                      - thetam*dt*(src_ec(i,j,k) + res_ec(i,j,k) &
                         + G0 - Lambda0)/eUn
                 res_HI(i,j,k) = HI(i,j,k) - HI0(i,j,k)          &
                      - theta*dt*(src_HI(i,j,k) + k2*ne*nHII     &
                         - nHI*(k1*ne + piHI(i,j,k)))/nUn        &
                      - thetam*dt*(src_HI(i,j,k) + k2o*ne0*nHII0 &
                         - nHI0*(k1o*ne0 + piHI(i,j,k)))/nUn0
                 
              enddo
           enddo
        enddo

     ! Hydrogen + Helium case        
     else 
        ! iterate over the domain
        do k=1,Nz,1
           do j=1,Ny,1
              do i=1,Nx,1
                 
                 ! set shortcut values for this spatial location,
                 ! converting densities from comoving to proper, and 
                 ! put all shortcuts in ACTUAL units
                 ecval   = ec(i,j,k)*eUn
                 rhoval  = rho(i,j,k)*dUn
                 rhoval0 = rho(i,j,k)*dUn0
                 nHI     = HI(i,j,k)*nUn
                 nHI0    = HI0(i,j,k)*nUn0
                 nH      = Hfrac*rhoval/mp
                 nH0     = Hfrac*rhoval0/mp
                 nHII    = max(1.d0*(nH - nHI), 0.d0)
                 nHII0   = max(1.d0*(nH0 - nHI0), 0.d0)
                 nHe     = (1.d0-HFrac)*rhoval/4.d0/mp
                 nHe0    = (1.d0-HFrac)*rhoval0/4.d0/mp
                 nHeI    = HeI(i,j,k)*nUn/4.d0
                 nHeI0   = HeI0(i,j,k)*nUn0/4.d0
                 nHeII   = HeII(i,j,k)*nUn/4.d0
                 nHeII0  = HeII0(i,j,k)*nUn0/4.d0
                 nHeIII  = max(1.d0*(nHe - nHeI - nHeII), 0.d0)
                 nHeIII0 = max(1.d0*(nHe0 - nHeI0 - nHeII0), 0.d0)
                 ne      = nHII + 0.25d0*nHeII + 0.5d0*nHeIII
                 ne0     = nHII0 + 0.25d0*nHeII0 + 0.5d0*nHeIII0
                 
                 ! compute temperatures
                 eint = vUn*vUn*(eh(i,j,k) &
                      - KEconst*(vx(i,j,k)**2 + vy(i,j,k)**2 + vz(i,j,k)**2))
                 mol_weight = rhoval/mp/(0.25d0*(nHeI+nHeII+nHeIII)+nHI+nHII+ne)
                 T = gam_1*mol_weight*mp*(ec(i,j,k)*eUn + eint)/kb
                 T = max(1.d0*T,1.d0*min_temp)
                 lamT = 3.15614d5/T
                 mol_weight = rhoval0/mp/(0.25d0*(nHeI0+nHeII0+nHeIII0)+nHI0+nHII0+ne0)
                 T0 = gam_1*mol_weight*mp*eint/kb
                 T0 = max(1.d0*T0,1.d0*min_temp)
                 lamT0 = 3.15614d5/T0

                 ! look up rates
                 lTemp = min(max(log(T), lTempS), lTempE)
                 Tidx = min(NTempBins-1, max(1, int((lTemp-lTempS)/dlTemp)+1))
                 Tidxp = Tidx+1
                 Tl = lTempS + (Tidx-1)*dlTemp
                 Tr = lTempS +  Tidx*dlTemp
                 Tfac = (lTemp - Tl)/(Tr - Tl)
                 k1 = k1Tb(Tidx) + (k1Tb(Tidxp) - k1Tb(Tidx))*Tfac
                 k3 = k3Tb(Tidx) + (k3Tb(Tidxp) - k3Tb(Tidx))*Tfac
                 k4 = k4Tb(Tidx) + (k4Tb(Tidxp) - k4Tb(Tidx))*Tfac
                 k5 = k5Tb(Tidx) + (k5Tb(Tidxp) - k5Tb(Tidx))*Tfac
                 k6 = k6Tb(Tidx) + (k6Tb(Tidxp) - k6Tb(Tidx))*Tfac
                 ceHI = ceHITb(Tidx) + (ceHITb(Tidxp) - ceHITb(Tidx))*Tfac
                 ceHeI = ceHeITb(Tidx) + (ceHeITb(Tidxp) - ceHeITb(Tidx))*Tfac
                 ceHeII = ceHeIITb(Tidx) + (ceHeIITb(Tidxp) - ceHeIITb(Tidx))*Tfac
                 ciHI = ciHITb(Tidx) + (ciHITb(Tidxp) - ciHITb(Tidx))*Tfac
                 ciHeI = ciHeITb(Tidx) + (ciHeITb(Tidxp) - ciHeITb(Tidx))*Tfac
                 ciHeII = ciHeIITb(Tidx) + (ciHeIITb(Tidxp) - ciHeIITb(Tidx))*Tfac
                 ciHeIS = ciHeISTb(Tidx) + (ciHeISTb(Tidxp) - ciHeISTb(Tidx))*Tfac
                 reHII = reHIITb(Tidx) + (reHIITb(Tidxp) - reHIITb(Tidx))*Tfac
                 reHeII1 = reHeII1Tb(Tidx) + (reHeII1Tb(Tidxp) - reHeII1Tb(Tidx))*Tfac
                 reHeII2 = reHeII2Tb(Tidx) + (reHeII2Tb(Tidxp) - reHeII2Tb(Tidx))*Tfac
                 reHeIII = reHeIIITb(Tidx) + (reHeIIITb(Tidxp) - reHeIIITb(Tidx))*Tfac
                 brem = bremTb(Tidx) + (bremTb(Tidxp) - bremTb(Tidx))*Tfac

                 lTemp = min(max(log(T0), lTempS), lTempE)
                 Tidx = min(NTempBins-1, max(1, int((lTemp-lTempS)/dlTemp)+1))
                 Tidxp = Tidx+1
                 Tl = lTempS + (Tidx-1)*dlTemp
                 Tr = lTempS +  Tidx*dlTemp
                 Tfac = (lTemp - Tl)/(Tr - Tl)
                 k1o = k1Tb(Tidx) + (k1Tb(Tidxp) - k1Tb(Tidx))*Tfac
                 k3o = k3Tb(Tidx) + (k3Tb(Tidxp) - k3Tb(Tidx))*Tfac
                 k4o = k4Tb(Tidx) + (k4Tb(Tidxp) - k4Tb(Tidx))*Tfac
                 k5o = k5Tb(Tidx) + (k5Tb(Tidxp) - k5Tb(Tidx))*Tfac
                 k6o = k6Tb(Tidx) + (k6Tb(Tidxp) - k6Tb(Tidx))*Tfac
                 ceHI0 = ceHITb(Tidx) + (ceHITb(Tidxp) - ceHITb(Tidx))*Tfac
                 ceHeI0 = ceHeITb(Tidx) + (ceHeITb(Tidxp) - ceHeITb(Tidx))*Tfac
                 ceHeII0 = ceHeIITb(Tidx) + (ceHeIITb(Tidxp) - ceHeIITb(Tidx))*Tfac
                 ciHI0 = ciHITb(Tidx) + (ciHITb(Tidxp) - ciHITb(Tidx))*Tfac
                 ciHeI0 = ciHeITb(Tidx) + (ciHeITb(Tidxp) - ciHeITb(Tidx))*Tfac
                 ciHeII0 = ciHeIITb(Tidx) + (ciHeIITb(Tidxp) - ciHeIITb(Tidx))*Tfac
                 ciHeIS0 = ciHeISTb(Tidx) + (ciHeISTb(Tidxp) - ciHeISTb(Tidx))*Tfac
                 reHII0 = reHIITb(Tidx) + (reHIITb(Tidxp) - reHIITb(Tidx))*Tfac
                 reHeII10 = reHeII1Tb(Tidx) + (reHeII1Tb(Tidxp) - reHeII1Tb(Tidx))*Tfac
                 reHeII20 = reHeII2Tb(Tidx) + (reHeII2Tb(Tidxp) - reHeII2Tb(Tidx))*Tfac
                 reHeIII0 = reHeIIITb(Tidx) + (reHeIIITb(Tidxp) - reHeIIITb(Tidx))*Tfac
                 brem0 = bremTb(Tidx) + (bremTb(Tidxp) - bremTb(Tidx))*Tfac

                 ! compute case B Hydrogen recombination coefficient 
                 ! [Hui & Gnedin, 1997: RI^B_{HII}]
                 k2 = 2.753d-14*lamT**(1.5d0) /                 &
                      (1.d0+(lamT/2.74d0)**(0.407d0))**(2.242d0)
                 k2o = 2.753d-14*lamT0**(1.5d0) /               &
                      (1.d0+(lamT0/2.74d0)**(0.407d0))**(2.242d0)
                 
                 ! compute fluid cooling rate.  Terms (in order):
                 !    Collisional Excitations (3)
                 !    Collisional Ionizations (4)
                 !    Recombinations (4)
                 !    Compton cooling or heating (1)
                 !    X-ray Compton heating (1)
                 !    Bremsstrahlung (1)
                 Lambda = ne/rhoval                    &
                      *(ceHI*nHI                       &
                      + ceHeI*nHeII*ne/4.d0            &
                      + ceHeII*nHeII/4.d0              &
                      + ciHI*nHI                       &
                      + ciHeI*nHeI/4.d0                &
                      + ciHeII*nHeII/4.d0              &
                      + ciHeIS*nHeII*ne/4.d0           &
                      + reHII*nHII                     &
                      + reHeII1*nHeII/4.d0             &
                      + reHeII2*nHeII/4.d0             &
                      + reHeIII*nHeIII/4.d0            &
                      + Comp1*(T-Comp2)                &
                      + Comp_xray*(T-Comp_temp)        &
                      + brem*(nHII+nHeII/4.d0+nHeIII)  &
                      )
                 Lambda0 = ne0/rhoval0                     &
                      *(ceHI0*nHI0                         &
                      + ceHeI0*nHeII0*ne0/4.d0             &
                      + ceHeII0*nHeII0/4.d0                &
                      + ciHI0*nHI0                         &
                      + ciHeI0*nHeI0/4.d0                  &
                      + ciHeII0*nHeII0/4.d0                &
                      + ciHeIS0*nHeII0*ne0/4.d0            &
                      + reHII0*nHII0                       &
                      + reHeII10*nHeII0/4.d0               &
                      + reHeII20*nHeII0/4.d0               &
                      + reHeIII0*nHeIII0/4.d0              &
                      + Comp1*(T0-Comp2)                   &
                      + Comp_xray*(T0-Comp_temp)           &
                      + brem0*(nHII0+nHeII0/4.d0+nHeIII0)  &
                      )
                 
                 ! compute fluid heating rates
                 G  = (nHI*GHI(i,j,k) + nHeI*GHeI(i,j,k)   &
                      + nHeII*GHeII(i,j,k))/rhoval
                 G0 = (nHI0*GHI(i,j,k) + nHeI0*GHeI(i,j,k) &
                      + nHeII0*GHeII(i,j,k))/rhoval0
                 
                 ! put it all together
                 res_ec(i,j,k) = ec(i,j,k)                       &
                      - theta*dt*(src_ec(i,j,k) + res_ec(i,j,k)  &
                         + G - Lambda - 2.d0*afac*ecval)/eUn     &
                      - thetam*dt*(src_ec(i,j,k) + res_ec(i,j,k) &
                         + G0 - Lambda0)/eUn
                 res_HI(i,j,k) = HI(i,j,k) - HI0(i,j,k)          &
                      - theta*dt*(src_HI(i,j,k) + k2*ne*nHII     &
                         - nHI*(k1*ne + piHI(i,j,k)))/nUn        &
                      - thetam*dt*(src_HI(i,j,k) + k2o*ne0*nHII0 &
                         - nHI0*(k1o*ne0 + piHI(i,j,k)))/nUn0
                 res_HeI(i,j,k) = HeI(i,j,k) - HeI0(i,j,k)         &
                      - theta*dt*(src_HeI(i,j,k) + k4*ne*nHeII     &
                         - nHeI*(k3*ne + piHeI(i,j,k)))/nUn        &
                      - thetam*dt*(src_HeI(i,j,k) + k4o*ne0*nHeII0 &
                         - nHeI0*(k3o*ne0 + piHeI(i,j,k)))/nUn0
                 res_HeII(i,j,k) = HeII(i,j,k) - HeII0(i,j,k)        &
                      - theta*dt*(src_HeII(i,j,k) + ne*k6*nHeIII     &
                         + nHeI*(ne*k3 + piHeI(i,j,k))               &
                         - nHeII*(ne*(k4+k5) + piHeII(i,j,k)))/nUn   &
                      - thetam*dt*(src_HeII(i,j,k) + ne0*k6o*nHeIII0 &
                         + nHeI0*(ne0*k3o + piHeI(i,j,k))            &
                         - nHeII0*(ne0*(k4o+k5o) + piHeII(i,j,k)))/nUn0

              enddo
           enddo
        enddo
     endif  ! Nchem

  !==================================================
  else

     write(0,*) 'MFProb_LocResid: Model =',Model,' undefined!'

  endif ! Model


  return
end subroutine MFProb_LocResid
!=======================================================================
