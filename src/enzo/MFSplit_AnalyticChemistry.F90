!=======================================================================
!
! Copyright 2009 Daniel R. Reynolds
! Copyright 2009 Laboratory for Computational Astrophysics
! Copyright 2009 Regents of the University of California
!
! This software is released under the terms of the "Enzo Public License"
! in the accompanying LICENSE file.
!
!=======================================================================
subroutine MFSplit_AnalyticChemistry(ec, HI, HeI, HeII, HI0, HeI0, HeII0, &
     dt, vx, vy, vz, rho, eh, src_ec, src_HI, src_HeI, src_HeII, piHI,    &
     piHeI, piHeII, GHI, GHeI, GHeII, gamma, HFrac, Model, DualEnergy, a, &
     a0, adot, adot0, CompA, CompX, CompT, NTempBins, TStart, TEnd, k1Tb, &
     k2Tb, k3Tb, k4Tb, k5Tb, k6Tb, ceHITb, ceHeITb, ceHeIITb, ciHITb,     &
     ciHeITb, ciHeISTb, ciHeIITb, reHIITb, reHeII1Tb, reHeII2Tb,          &
     reHeIIITb, bremTb, aUn, dUn, dUn0, vUn, lUn, lUn0, eUn, nUn, nUn0,   &
     ecScale, Nchem, Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr, NGzl, NGzr, ier)
  !=======================================================================
  !  written by: Daniel R. Reynolds
  !  date:       December 2009
  !
  !  PURPOSE: Computes the solutions to the gas energy correction and 
  !           chemistry equations, using a quasi-steady-state 
  !           approximation and the resulting analytical solution of the 
  !           relevant ODEs.
  !
  !  INPUTS:
  !     ec         - fluid energy correction array
  !     HI         - Hydrogen I number density array
  !     HeI        - Helium I number density array
  !     HeII       - Helium II number density array
  !     HI0        - Hydrogen I number density array (old)
  !     HeI0       - Helium I number density array (old)
  !     HeII0      - Helium II number density array (old)
  !     dt         - time step size
  !     vx,vy,vz   - velocity arrays in each direction
  !     rho        - density array
  !     eh         - total fluid energy array
  !     src_ec     - source function values for gas energy correction eq.
  !     src_HI     - source function values for HI eq.
  !     src_HeI    - source function values for HeI eq.
  !     src_HeII   - source function values for HeII eq.
  !     piH*       - chemistry photo-ionization rates
  !     GH*        - chemistry photo-heating coefficients
  !     gamma      - constant in ideal gas law
  !     HFrac      - percentage of mass composed of Hydrogen
  !     Model      - flag denoting physical model to use
  !     DualEnergy - flag denoting dual energy formalism
  !     a,a0       - cosmological expansion parameter (new and old)
  !     adot,adot0 - da/dt (new and old)
  !     CompA      - Compton cooling coefficient 1 (multiplier)
  !     CompX      - X-ray Compton heating coefficient
  !     CompT      - X-ray Compton heating temperature 
  !     NTempBins  - number of temperature bins for rate tables
  !     TStart     - starting temperature for rate table
  !     TEnd       - ending temperature for rate table
  !     k*Tb       - chemistry rate tables
  !     c*Tb, r*Tb - heating/cooling rate tables
  !     *Un,*Un0   - variable scaling constants (new and old)
  !     ecScale    - for model 4, this holds the constant temperature
  !     Nchem      - number of chemical species in simulation
  !     Nx,Ny,Nz   - active mesh size in each direction
  !     NG*l/NG*r  - left/right ghost cells in each direction
  !
  !     Note: the vector inputs are of size (Nx + NGxl + NGxr) in 
  !     the x-direction, others are similar.
  !
  !  OUTPUT ARGUMENTS: 
  !     ec         - updated fluid energy correction array
  !     HI         - updated Hydrogen I number density array
  !     HeI        - updated Helium I number density array
  !     HeII       - updated Helium II number density array
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
  integer, intent(in)  :: Model, Nchem, NTempBins, DualEnergy
  integer, intent(in)  :: Nx, NGxl, NGxr
  integer, intent(in)  :: Ny, NGyl, NGyr
  integer, intent(in)  :: Nz, NGzl, NGzr
  integer, intent(out) :: ier
  REALSUB, intent(in)  :: a, a0, adot, adot0
  real, intent(in) :: dt, gamma, HFrac, TStart, TEnd
  real, intent(in) :: CompA, CompX, CompT
  real, intent(in) :: dUn, dUn0, vUn, lUn, lUn0, &
       eUn, nUn, nUn0, aUn, ecScale
  real, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr) :: &
       ec, HI, HeI, HeII
  real, intent(in),                                             &
       dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr) ::   &
       HI0, HeI0, HeII0, vx, vy, vz, rho, eh, src_ec, src_HI,       &
       src_HeI, src_HeII, piHI, piHeI, piHeII, GHI, GHeI, GHeII
  real, intent(in), dimension(NTempBins) :: k1Tb, k2Tb, k3Tb,    &
       k4Tb, k5Tb, k6Tb, ceHITb, ceHeITb, ceHeIITb, ciHITb, ciHeITb, &
       ciHeIITb, ciHeISTb, reHIITb, reHeII1Tb, reHeII2Tb, reHeIIITb, &
       bremTb

  !--------------
  ! locals
  integer :: i, j, k, l, sweeps, sweeps2
  real :: Temp, mp, kb, mol_weight, ERes
  real :: nH, nHI, nHII, ne, nHe, nHeI, nHeII, nHeIII, rhoval
  real :: aval, dadt, dUnit, lUnit, nUnit, gam_1
  real :: zr, Comp1, Comp2, KEconst, eint, dx
  real :: res_Er, res_ec, res_HI, res_HeI, res_HeII
  real :: ecnew, HInew, HeInew, HeIInew, change, lam, lam2, FPtol
  
  !=======================================================================

  ! initialize success/fail flag to success
  ier = 1

  ! we only have this enabled for Models 1 & 4 (case B HII recomb rate), 
  ! with chemistry 
  if ((Nchem < 1) .or. ((Model /= 1) .and. (Model /= 4))) then
     write(*,*) 'AnalyticResid ERROR: only implemented for Models 1 & 4', &
          ' w/ chem, called with Model =',Model,' and Nchem = ',Nchem
     ier = 0
     return
  endif

  ! compute shortcuts
  dx = 1.d0       ! send in a dummy value, since residual does not use Er
  mp = 1.67262171d-24          ! Mass of a proton [g]
  kb = 1.3806504d-16           ! boltzmann constant [erg/K]
  gam_1 = gamma-1.d0
  aval = (a+a0)*0.5d0
  dadt = (adot+adot0)*0.5d0
  dUnit = (dUn+dUn0)*0.5d0
  lUnit = (lUn+lUn0)*0.5d0
  nUnit = (nUn+nUn0)*0.5d0
  Comp1 = 0.d0
  Comp2 = 0.d0
  zr = 1.d0/(a*aUn) - 1.d0
  if (a*aUn /= 1.d0) then
     Comp1 = CompA*(1.d0 + zr)**4
     Comp2 = 2.73d0*(1.d0 + zr)
  endif
  KEconst = 0.5d0
  if (DualEnergy == 1)  KEconst = 0.d0

  ! set fixed-point iteration parameters
  FPtol = 1.0d-8
  sweeps = 50
  sweeps2 = 10000
  lam = 1.d0
  lam2 = 0.1d0


!!$  ! skip chemistry and heating for testing
!!$  return


  ! perform iteration based on chemistry

  !    Hydrogen chemistry
  if (Nchem == 1) then

     ! iterate over the domain
     do k=1,Nz
        do j=1,Ny
           do i=1,Nx
              
              ! get shortcut values at this location
              ecnew = ec(i,j,k)
              HInew = HI(i,j,k)
              eint = vUn*vUn*(eh(i,j,k)                                  &
                   - KEconst*(vx(i,j,k)**2 + vy(i,j,k)**2 + vz(i,j,k)**2))
              
              ! perform coarse fixed-point iteration to find analytical solution
              do l=0,sweeps
                 
                 ! compute local temperature (based on Model)
                 if (Model == 4) then
                    if (adot == 0.d0) then
                       mol_weight = 0.6d0
                       Temp = (gamma-1.d0)*0.6d0*mp*eint/kb
                       Temp = max(1.d0*Temp,1.d-1)
                    else
                       Temp = ecScale
                    endif
                 else 
                    rhoval = rho(i,j,k)*dUnit
                    nHI  = HI(i,j,k)*nUnit
                    nH  = Hfrac*rhoval/mp
                    nHII  = max(1.d0*(nH - nHI), 0.d0)
                    ne  = nHII
                    mol_weight = rhoval/mp/(nHI + nHII + ne)
                    Temp = gam_1*mol_weight*mp*(ec(i,j,k)*eUn*0.5d0 + eint)/kb
                    Temp = max(1.d0*Temp,1.d-1)
                 endif
                    
                 ! call the local residual routine
                 call MFSplit_AnalyticChemResid(res_HI, res_HeI, res_HeII,       &
                      HInew, 0.d0, 0.d0, HI0(i,j,k), 0.d0, 0.d0, dt, rho(i,j,k), &
                      Temp, src_HI(i,j,k), 0.d0, 0.d0, piHI(i,j,k), 0.d0, 0.d0,  &
                      NTempBins, TStart, TEnd, k1Tb, k2Tb, k3Tb, k4Tb, k5Tb,     &
                      k6Tb, dUnit, nUnit, Nchem, HFrac, ier)
                 call MFSplit_AnalyticGasResid(res_ec, ecnew, HInew, 0.d0, 0.d0, &
                      HI0(i,j,k), 0.d0, 0.d0, dt, rho(i,j,k), eint,              &
                      src_ec(i,j,k), gamma, HFrac, Model, a, adot, Comp1, Comp2, &
                      CompX, CompT, GHI(i,j,k), 0.d0, 0.d0, NTempBins, TStart,   &
                      TEnd, ceHITb, ceHeITb, ceHeIITb, ciHITb, ciHeITb,          &
                      ciHeISTb, ciHeIITb, reHIITb, reHeII1Tb, reHeII2Tb,         &
                      reHeIIITb, bremTb, dUnit, eUn, nUnit, Nchem, ier)
                 
                 ! check the error flag
                 if (ier /= 1)  return
                 
                 ! update the current guesses
                 ecnew = ecnew - lam*res_ec
                 HInew = HInew - lam*res_HI
                 
                 ! compute the FP change in iteration, check for convergence
                 change = max(abs(res_ec),abs(res_HI))
                 if (change < FPtol)  exit
                 
              end do
              
              ! if the coarse iteration unsuccessful, repeat with a finer update
              if (change >= FPtol) then
                 
                 ! perform fine fixed-point iteration
                 do l=0,sweeps2
                    
                    ! compute local temperature (based on Model)
                    if (Model == 4) then
                       if (adot == 0.d0) then
                          mol_weight = 0.6d0
                          Temp = (gamma-1.d0)*0.6d0*mp*eint/kb
                          Temp = max(1.d0*Temp,1.d-1)
                       else
                          Temp = ecScale
                       endif
                    else 
                       rhoval = rho(i,j,k)*dUnit
                       nHI  = HI(i,j,k)*nUnit
                       nH  = Hfrac*rhoval/mp
                       nHII  = max(1.d0*(nH - nHI), 0.d0)
                       ne  = nHII
                       mol_weight = rhoval/mp/(nHI + nHII + ne)
                       Temp = gam_1*mol_weight*mp*(ec(i,j,k)*eUn*0.5d0 + eint)/kb
                       Temp = max(1.d0*Temp,1.d-1)
                    endif
                    
                    ! call the local residual routine
                    call MFSplit_AnalyticChemResid(res_HI, res_HeI, res_HeII,    &
                         HInew, 0.d0, 0.d0, HI0(i,j,k), 0.d0, 0.d0, dt,          &
                         rho(i,j,k), Temp, src_HI(i,j,k), 0.d0, 0.d0,            &
                         piHI(i,j,k), 0.d0, 0.d0, NTempBins, TStart, TEnd, k1Tb, &
                         k2Tb, k3Tb, k4Tb, k5Tb, k6Tb, dUnit, nUnit, Nchem,      &
                         HFrac, ier)
                    call MFSplit_AnalyticGasResid(res_ec, ecnew, HInew, 0.d0,    &
                         0.d0, HI0(i,j,k), 0.d0, 0.d0, dt, rho(i,j,k), eint,     &
                         src_ec(i,j,k), gamma, HFrac, Model, a, adot, Comp1,     &
                         Comp2, CompX, CompT, GHI(i,j,k), 0.d0, 0.d0, NTempBins, &
                         TStart, TEnd, ceHITb, ceHeITb, ceHeIITb, ciHITb,        &
                         ciHeITb, ciHeISTb, ciHeIITb, reHIITb, reHeII1Tb,        &
                         reHeII2Tb, reHeIIITb, bremTb, dUnit, eUn, nUnit, Nchem, &
                         ier)
                    
                    ! check the error flag
                    if (ier /= 1)  return
                    
                    ! update the current guesses
                    ecnew = ecnew - lam2*res_ec
                    HInew = HInew - lam2*res_HI
                    
                    ! compute the FP change in iteration, check for convergence
                    change = max(abs(res_ec),abs(res_HI))
                    if (change < FPtol)  exit
                    
                 end do
                 
              end if
              
              ! fill the relevant residuals
              ec(i,j,k) = ecnew
              HI(i,j,k) = HInew
              
           end do
        end do
     end do


  !    Hydrogen + Helium chemistry
  else

     ! iterate over the domain
     do k=1,Nz
        do j=1,Ny
           do i=1,Nx
              
              ! get shortcut values at this location
              ecnew   = ec(i,j,k)
              HInew   = HI(i,j,k)
              HeInew  = HeI(i,j,k)
              HeIInew = HeII(i,j,k)
              eint = vUn*vUn*(eh(i,j,k)                                  &
                   - KEconst*(vx(i,j,k)**2 + vy(i,j,k)**2 + vz(i,j,k)**2))
              
              ! perform coarse fixed-point iteration to find analytical solution
              do l=0,sweeps
                 
                 ! compute local temperature (based on Model)
                 if (Model == 4) then
                    if (adot == 0.d0) then
                       mol_weight = 0.6d0
                       Temp = (gamma-1.d0)*0.6d0*mp*eint/kb
                       Temp = max(1.d0*Temp,1.d-1)
                    else
                       Temp = ecScale
                    endif
                 else 
                    rhoval = rho(i,j,k)*dUnit
                    nHI  = HI(i,j,k)*nUnit
                    nH  = Hfrac*rhoval/mp
                    nHII = max(1.d0*(nH - nHI), 0.d0)
                    nHe = (1.d0-HFrac)*rhoval/4.d0/mp
                    nHeI = HeI(i,j,k)*nUnit/4.d0
                    nHeII = HeII(i,j,k)*nUnit/4.d0
                    nHeIII = max(1.d0*(nHe - nHeI - nHeII), 0.d0)
                    ne = nHII + 0.25d0*nHeII + 0.5d0*nHeIII
                    mol_weight = rhoval/mp/(0.25d0*(nHeI+nHeII+nHeIII)+nHI+nHII+ne)
                    Temp = gam_1*mol_weight*mp*(ec(i,j,k)*eUn*0.5d0 + eint)/kb
                    Temp = max(1.d0*Temp,1.d-1)
                 endif

                 ! call the local residual routines
                 call MFSplit_AnalyticChemResid(res_HI, res_HeI, res_HeII,       &
                      HInew, HeInew, HeIInew, HI0(i,j,k), HeI0(i,j,k),           &
                      HeII0(i,j,k), dt, rho(i,j,k), Temp, src_HI(i,j,k),         &
                      src_HeI(i,j,k), src_HeII(i,j,k), piHI(i,j,k),              &
                      piHeI(i,j,k), piHeII(i,j,k), NTempBins, TStart, TEnd,      &
                      k1Tb, k2Tb, k3Tb, k4Tb, k5Tb, k6Tb, dUnit, nUnit, Nchem,   &
                      HFrac, ier)
                 call MFSplit_AnalyticGasResid(res_ec, ecnew, HInew, HeInew,     &
                      HeIInew, HI0(i,j,k), HeI0(i,j,k), HeII0(i,j,k), dt,        &
                      rho(i,j,k), eint, src_ec(i,j,k), gamma, HFrac, Model, a,   &
                      adot, Comp1, Comp2, CompX, CompT, GHI(i,j,k), GHeI(i,j,k), &
                      GHeII(i,j,k), NTempBins, TStart, TEnd, ceHITb, ceHeITb,    &
                      ceHeIITb, ciHITb, ciHeITb, ciHeISTb, ciHeIITb, reHIITb,    &
                      reHeII1Tb, reHeII2Tb, reHeIIITb, bremTb, dUnit, eUn,       &
                      nUnit, Nchem, ier)
                 
                 ! check the error flag
                 if (ier /= 1)  return
                 
                 ! update the current guesses
                 ecnew   = ecnew   - lam*res_ec
                 HInew   = HInew   - lam*res_HI
                 HeInew  = HeInew  - lam*res_HeI
                 HeIInew = HeIInew - lam*res_HeII

                 ! compute the FP change in iteration, check for convergence
                 change = max(abs(res_ec),abs(res_HI))
                 change = max(change,abs(res_HeI))
                 change = max(change,abs(res_HeII))
                 if (change < FPtol)  exit
                 
              end do
              
              ! if the coarse iteration unsuccessful, repeat with a finer update
              if (change >= FPtol) then
                 
                 ! perform fine fixed-point iteration
                 do l=0,sweeps2
                    
                    ! compute local temperature (based on Model)
                    if (Model == 4) then
                       if (adot == 0.d0) then
                          mol_weight = 0.6d0
                          Temp = (gamma-1.d0)*0.6d0*mp*eint/kb
                          Temp = max(1.d0*Temp,1.d-1)
                       else
                          Temp = ecScale
                       endif
                    else 
                       rhoval = rho(i,j,k)*dUnit
                       nHI  = HI(i,j,k)*nUnit
                       nH  = Hfrac*rhoval/mp
                       nHII = max(1.d0*(nH - nHI), 0.d0)
                       nHe = (1.d0-HFrac)*rhoval/4.d0/mp
                       nHeI = HeI(i,j,k)*nUnit/4.d0
                       nHeII = HeII(i,j,k)*nUnit/4.d0
                       nHeIII = max(1.d0*(nHe - nHeI - nHeII), 0.d0)
                       ne = nHII + 0.25d0*nHeII + 0.5d0*nHeIII
                       mol_weight = rhoval/mp/(0.25d0*(nHeI+nHeII+nHeIII)+nHI+nHII+ne)
                       Temp = gam_1*mol_weight*mp*(ec(i,j,k)*eUn*0.5d0 + eint)/kb
                       Temp = max(1.d0*Temp,1.d-1)
                    endif

                    ! call the local residual routines
                    call MFSplit_AnalyticChemResid(res_HI, res_HeI, res_HeII,    &
                         HInew, HeInew, HeIInew, HI0(i,j,k), HeI0(i,j,k),        &
                         HeII0(i,j,k), dt, rho(i,j,k), Temp, src_HI(i,j,k),      &
                         src_HeI(i,j,k), src_HeII(i,j,k), piHI(i,j,k),           &
                         piHeI(i,j,k), piHeII(i,j,k), NTempBins, TStart, TEnd,   &
                         k1Tb, k2Tb, k3Tb, k4Tb, k5Tb, k6Tb, dUnit, nUnit,       &
                         Nchem, HFrac, ier)
                    call MFSplit_AnalyticGasResid(res_ec, ecnew, HInew, HeInew,  &
                         HeIInew, HI0(i,j,k), HeI0(i,j,k), HeII0(i,j,k), dt,     &
                         rho(i,j,k), eint, src_ec(i,j,k), gamma, HFrac, Model,   &
                         a, adot, Comp1, Comp2, CompX, CompT, GHI(i,j,k),        &
                         GHeI(i,j,k), GHeII(i,j,k), NTempBins, TStart, TEnd,     &
                         ceHITb, ceHeITb, ceHeIITb, ciHITb, ciHeITb, ciHeISTb,   &
                         ciHeIITb, reHIITb, reHeII1Tb, reHeII2Tb, reHeIIITb,     &
                         bremTb, dUnit, eUn, nUnit, Nchem, ier)
                    
                    ! check the error flag
                    if (ier /= 1)  return
                    
                    ! update the current guesses
                    ecnew   = ecnew   - lam2*res_ec
                    HInew   = HInew   - lam2*res_HI
                    HeInew  = HeInew  - lam2*res_HeI
                    HeIInew = HeIInew - lam2*res_HeII
                    
                    ! compute the FP change in iteration, check for convergence
                    change = max(abs(res_ec),abs(res_HI))
                    change = max(change,abs(res_HeI))
                    change = max(change,abs(res_HeII))
                    if (change < FPtol)  exit
                    
                 end do
                 
              end if
              
              ! fill the relevant residuals
              ec(i,j,k)   = ecnew
              HI(i,j,k)   = HInew
              HeI(i,j,k)  = HeInew
              HeII(i,j,k) = HeIInew
              
           end do
        end do
     end do

  end if

  ! exit subroutine
  return

end subroutine MFSplit_AnalyticChemistry
!=======================================================================






subroutine MFSplit_AnalyticInitGuess(Ef, E1, E2, E3, ec, HI, HeI,         &
     HeII, dt, vx, vy, vz, rho, eh, src_E1, src_E2, src_E3, src_ec,       &
     src_HI, src_HeI, src_HeII, gamma, HFrac, Model, DualEnergy,          &
     ESpectrum, a, a0, adot, adot0, CompA, CompX, CompT, piHI, piHeI,     &
     piHeII, GHI, GHeI, GHeII, NTempBins, TStart, TEnd, k1Tb, k2Tb, k3Tb, &
     k4Tb, k5Tb, k6Tb, ceHITb, ceHeITb, ceHeIITb, ciHITb, ciHeITb,        &
     ciHeISTb, ciHeIITb, reHIITb, reHeII1Tb, reHeII2Tb, reHeIIITb,        &
     bremTb, chibar, aUn, dUn, dUn0, vUn, lUn, lUn0, fsUn, eUn, nUn,      &
     nUn0, ecScale, Nchem, dx, dy, dz, Nx, Ny, Nz, NGxl, NGxr, NGyl,      &
     NGyr, NGzl, NGzr, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       August 2009
!
!  PURPOSE: Computes an initial guess to the time-evolved PDE system 
!     using a Gauss-Seidel iteration on the local ODE reaction system.
!
!  INPUTS:
!     Ef         - free-streaming radiation energy (new & old)
!     E1         - radiation energy freq 1 (old time on input, new on output)
!     E1         - radiation energy freq 1 (old time on input, new on output)
!     E2         - radiation energy freq 2
!     E3         - radiation energy freq 3
!     ec         - fluid energy correction
!     HI         - Hydrogen I number density
!     HeI        - Helium I number density
!     HeII       - Helium II number density
!     dt         - time step size
!     vx,vy,vz   - velocity arrays in each direction
!     rho        - density array
!     eh         - total fluid energy array
!     src_*      - source function values for equations
!     gamma      - constant in ideal gas law
!     HFrac      - percentage of mass composed of Hydrogen
!     Model      - flag denoting physical model to use
!     DualEnergy - flag denoting dual energy formalism
!     ESpectrum  - flag denoting type of free-streaming radiation spectrum
!     a,a0       - cosmological expansion parameter (new & old)
!     adot,adot0 - da/dt (new & old)
!     CompA      - Compton cooling coefficient 1 (multiplier)
!     CompX      - X-ray Compton heating coefficient
!     CompT      - X-ray Compton heating temperature 
!     piH*       - chemistry photo-ionization rates
!     GH*        - chemistry photo-heating coefficients
!     NTempBins  - number of temperature bins for rate tables
!     TStart     - starting temperature for rate table
!     TEnd       - ending temperature for rate table
!     k*Tb       - chemistry rate tables
!     c*Tb, r*Tb - heating/cooling rate tables
!     *Un,*Un0   - variable scaling constants (new & old)
!     ecScale    - holds the temperature for model 4
!     Nchem      - number of chemical species in simulation
!     Nx,Ny,Nz   - active mesh size in each direction
!     dx,dy,dz   - spatial mesh cell size in each direction
!     NG*l/NG*r  - left/right ghost cells in each direction
!
!     Note: the vector inputs are of size (Nx + NGxl + NGxr) in 
!     the x-direction, others are similar.
!
!  OUTPUT ARGUMENTS: 
!     E1   - guess at radiation energy density, freq 1
!     E2   - guess at radiation energy density, freq 2
!     E3   - guess at radiation energy density, freq 3
!     ec   - guess at gas energy correction 
!     HI   - guess at HI density (unfilled if Nchem < 1)
!     HeI  - guess at HeI density (unfilled if Nchem < 3)
!     HeII - guess at HeII density (unfilled if Nchem < 3)
!     ier  - success/failure flag (0->failure, 1->success)
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
  integer, intent(in)  :: Model, Nchem, NTempBins, DualEnergy, ESpectrum
  integer, intent(in)  :: Nx, NGxl, NGxr
  integer, intent(in)  :: Ny, NGyl, NGyr
  integer, intent(in)  :: Nz, NGzl, NGzr
  integer, intent(out) :: ier
  REALSUB, intent(in)  :: a, a0, adot, adot0
  real, intent(in) :: dt, dx, dy, dz, gamma, HFrac, TStart, TEnd
  real, intent(in) :: CompA, CompX, CompT, chibar
  real, intent(in) :: dUn, dUn0, vUn, lUn, lUn0, fsUn, &
       eUn, nUn, nUn0, aUn, ecScale
  real, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr), &
       intent(in) :: vx, vy, vz, rho, eh, src_E1, src_E2, src_E3,    &
       src_ec, src_HI, src_HeI, src_HeII
  real, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr) :: &
       ec, Ef, E1, E2, E3, HI, HeI, HeII, piHI, piHeI, piHeII,    &
       GHI, GHeI, GHeII
  real, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr) &
       :: ecres, HIres, HeIres, HeIIres
  real, intent(in), dimension(NTempBins) :: k1Tb, k2Tb, k3Tb,    &
       k4Tb, k5Tb, k6Tb, ceHITb, ceHeITb, ceHeIITb, ciHITb, ciHeITb, &
       ciHeIITb, ciHeISTb, reHIITb, reHeII1Tb, reHeII2Tb, reHeIIITb, &
       bremTb

!--------------
! locals
  integer  :: i, j, k, l, l2, sweeps, sweeps2
  real :: eint, mp, kb, mol_weight, Temp, ERes, dtmp
  real :: nH, nHI, nHII, ne, nHe, nHeI, nHeII, nHeIII
  real :: FPtol, lam, lam2, change, zr, Comp1, Comp2, KEconst
  real :: deltax, gam_1
  real :: rhoval, aval, dadt, dUnit, lUnit, nUnit
  real :: res_E1, res_E2, res_E3, res_ec, res_HI, res_HeI, res_HeII
  real :: E1new, E2new, E3new, ecnew, HInew, HeInew, HeIInew
  
!=======================================================================

  ! initialize success/fail flag to success
  ier = 1

  ! we only have this enabled for Model 1 (case B HII recomb rate), 
  ! with chemistry 
  if ((Nchem == 0) .or. (Model /= 1)) then
     write(*,*) 'AnalyticInitGuess ERROR: only implemented for Model 1, w/ chem', &
          'called with Model =',Model,' and Nchem = ',Nchem
     ier = 0
     return
  endif

  ! compute shortcuts
  deltax = (dx*dy*dz)**(1.d0/3.d0)
  aval = (a+a0)*0.5d0
  dadt = (adot+adot0)*0.5d0
  dUnit = (dUn+dUn0)*0.5d0
  lUnit = (lUn+lUn0)*0.5d0
  nUnit = (nUn+nUn0)*0.5d0
  Comp1 = 0.d0
  Comp2 = 0.d0
  gam_1 = gamma-1.d0
  mp = 1.67262171d-24          ! Mass of a proton [g]
  kb = 1.3806504d-16           ! boltzmann constant [erg/K]
  zr = 1.d0/(a*aUn) - 1.d0
  if (a*aUn /= 1.d0) then
     Comp1 = CompA*(1.d0 + zr)**4
     Comp2 = 2.73d0*(1.d0 + zr)
  endif
  KEconst = 0.d0
  if (DualEnergy == 1)  KEconst = 0.d0

  ! set some parameters for the Gauss-Seidel iteration
  FPtol = 1.0d-8
  sweeps = 50
  sweeps2 = 10000
  lam = 1.d0
  lam2 = 0.1d0


  ! perform method based on chemistry
  !    Hydrogen chemistry
  if (Nchem == 1) then
  
     ! iterate over the domain
     do k=1,Nz
        do j=1,Ny
           do i=1,Nx
              
              ! get shortcut values at this location
              E1new = E1(i,j,k)
              E2new = E2(i,j,k)
              E3new = E3(i,j,k)
              ecnew = ec(i,j,k)
              HInew = HI(i,j,k)
              eint = vUn*vUn*(eh(i,j,k)                                 &
                   - KEconst*(vx(i,j,k)**2 + vy(i,j,k)**2 + vz(i,j,k)**2))
              
              ! initialize iteration counters
              l = 0
              l2 = 0
              ERes = 1.d0

              ! perform coarse fixed-point iteration
              do l=0,sweeps
                 
                 if (ERes > FPTol) then
                    ! compute the photo-ionization and photo-heating coeffs
                    call MFSplit_ComputeRadiationIntegrals(piHI(i,j,k), dtmp,    &
                         dtmp, GHI(i,j,k), dtmp, dtmp, Ef(i,j,k), E1(i,j,k),     &
                         E2(i,j,k), E3(i,j,k), E1(i,j,k), E2(i,j,k), E3(i,j,k),  &
                         Nchem, ESpectrum, chibar, fsUn, 1, 1, 1, 0, 0, 0, 0, 0, &
                         0, ier)
                 endif

                 ! compute local temperature (based on Model)
                 if (Model == 4) then
                    if (adot == 0.d0) then
                       mol_weight = 0.6d0
                       Temp = (gamma-1.d0)*0.6d0*mp*eint/kb
                       Temp = max(1.d0*Temp,1.d-1)
                    else
                       Temp = ecScale
                    endif
                 else 
                    rhoval = rho(i,j,k)*dUnit
                    nHI  = HI(i,j,k)*nUnit
                    nH  = Hfrac*rhoval/mp
                    nHII  = max(1.d0*(nH - nHI), 0.d0)
                    ne  = nHII
                    mol_weight = rhoval/mp/(nHI + nHII + ne)
                    Temp = gam_1*mol_weight*mp*(ec(i,j,k)*eUn*0.5d0 + eint)/kb
                    Temp = max(1.d0*Temp,1.d-1)
                 endif

                 ! call the local residual routines
                 call MFSplit_AnalyticRadResid(res_E1, res_E2, res_E3, E1new,    &
                      E2new, E3new, HInew, 0.d0, 0.d0, E1(i,j,k), E2(i,j,k),     &
                      E3(i,j,k), HI(i,j,k), 0.d0, 0.d0, dt, dx, src_E1(i,j,k),   &
                      src_E2(i,j,k), src_E3(i,j,k), a, lUnit, nUnit, Nchem, ier)
                 call MFSplit_AnalyticChemResid(res_HI, res_HeI, res_HeII,       &
                      HInew, 0.d0, 0.d0, HI(i,j,k), 0.d0, 0.d0, dt, rho(i,j,k),  &
                      Temp, src_HI(i,j,k), 0.d0, 0.d0, piHI(i,j,k), 0.d0, 0.d0,  &
                      NTempBins, TStart, TEnd, k1Tb, k2Tb, k3Tb, k4Tb, k5Tb,     &
                      k6Tb, dUnit, nUnit, Nchem, HFrac, ier)
                 call MFSplit_AnalyticGasResid(res_ec, ecnew, HInew, 0.d0, 0.d0, &
                      HI(i,j,k), 0.d0, 0.d0, dt, rho(i,j,k), eint, src_ec(i,j,k),&
                      gamma, HFrac, Model, a, adot, Comp1, Comp2, CompX, CompT,  &
                      GHI(i,j,k), 0.d0, 0.d0, NTempBins, TStart, TEnd, ceHITb,   &
                      ceHeITb, ceHeIITb, ciHITb, ciHeITb, ciHeISTb, ciHeIITb,    &
                      reHIITb, reHeII1Tb, reHeII2Tb, reHeIIITb, bremTb, dUnit,   &
                      eUn, nUnit, Nchem, ier)
                 
                 ! check the error flag
                 if (ier /= 1)  return
                 
                 ! compute FP change in this iteration and check for convergence
                 ERes = max(abs(res_E1),abs(res_E2))
                 ERes = max(ERes,abs(res_E3))
                 change = max(ERes,abs(res_ec))
                 change = max(change,abs(res_HI))
                 if (change < FPtol)  exit
                 
                 ! update the current guesses
                 E1new = E1new - lam*res_E1
                 E2new = E2new - lam*res_E2
                 E3new = E3new - lam*res_E3
                 ecnew = ecnew - lam*res_ec
                 HInew = HInew - lam*res_HI
                 
              end do  ! coarse fixed-point iteration
              
              ! if the coarse iteration unsuccessful, repeat with a finer update
              if (change >= FPtol) then
                 
                 ! perform fine fixed-point iteration
                 do l2=0,sweeps2
                    
                    if (ERes > FPTol) then
                       ! compute the photo-ionization and photo-heating coeffs
                       call MFSplit_ComputeRadiationIntegrals(piHI(i,j,k), dtmp, &
                            dtmp, GHI(i,j,k), dtmp, dtmp, Ef(i,j,k), E1(i,j,k),  &
                            E2(i,j,k), E3(i,j,k), E1(i,j,k), E2(i,j,k),          &
                            E3(i,j,k), Nchem, ESpectrum, chibar, fsUn, 1, 1, 1,  &
                            0, 0, 0, 0, 0, 0, ier)
                    endif
                    
                    ! compute local temperature (based on Model)
                    if (Model == 4) then
                       if (adot == 0.d0) then
                          mol_weight = 0.6d0
                          Temp = (gamma-1.d0)*0.6d0*mp*eint/kb
                          Temp = max(1.d0*Temp,1.d-1)
                       else
                          Temp = ecScale
                       endif
                    else 
                       rhoval = rho(i,j,k)*dUnit
                       nHI  = HI(i,j,k)*nUnit
                       nH  = Hfrac*rhoval/mp
                       nHII  = max(1.d0*(nH - nHI), 0.d0)
                       ne  = nHII
                       mol_weight = rhoval/mp/(nHI + nHII + ne)
                       Temp = gam_1*mol_weight*mp*(ec(i,j,k)*eUn*0.5d0 + eint)/kb
                       Temp = max(1.d0*Temp,1.d-1)
                    endif
                    
                    ! call the local residual routines
                    call MFSplit_AnalyticRadResid(res_E1, res_E2, res_E3, E1new, &
                         E2new, E3new, HInew, 0.d0, 0.d0, E1(i,j,k), E2(i,j,k),  &
                         E3(i,j,k), HI(i,j,k), 0.d0, 0.d0, dt, dx, src_E1(i,j,k),&
                         src_E2(i,j,k), src_E3(i,j,k), a, lUnit, nUnit, Nchem, ier)
                    call MFSplit_AnalyticChemResid(res_HI, res_HeI, res_HeII,    &
                         HInew, 0.d0, 0.d0, HI(i,j,k), 0.d0, 0.d0, dt,           &
                         rho(i,j,k), Temp, src_HI(i,j,k), 0.d0, 0.d0,            &
                         piHI(i,j,k), 0.d0, 0.d0, NTempBins, TStart, TEnd, k1Tb, &
                         k2Tb, k3Tb, k4Tb, k5Tb, k6Tb, dUnit, nUnit, Nchem,      &
                         HFrac, ier)
                    call MFSplit_AnalyticGasResid(res_ec, ecnew, HInew, 0.d0,    &
                         0.d0, HI(i,j,k), 0.d0, 0.d0, dt, rho(i,j,k), eint,      &
                         src_ec(i,j,k), gamma, HFrac, Model, a, adot, Comp1,     &
                         Comp2, CompX, CompT, GHI(i,j,k), 0.d0, 0.d0, NTempBins, &
                         TStart, TEnd, ceHITb, ceHeITb, ceHeIITb, ciHITb,        &
                         ciHeITb, ciHeISTb, ciHeIITb, reHIITb, reHeII1Tb,        &
                         reHeII2Tb, reHeIIITb, bremTb, dUnit, eUn, nUnit, Nchem, &
                         ier)
                    
                    ! check the error flag
                    if (ier /= 1)  return
                    
                    ! compute FP change in iteration and check for convergence
                    ERes = max(abs(res_E1),abs(res_E2))
                    ERes = max(ERes,abs(res_E3))
                    change = max(ERes,abs(res_ec))
                    change = max(change,abs(res_HI))
                    if (change < FPtol)  exit
                    
                    ! update the current guesses
                    E1new = E1new - lam2*res_E1
                    E2new = E2new - lam2*res_E2
                    E3new = E3new - lam2*res_E3
                    ecnew = ecnew - lam2*res_ec
                    HInew = HInew - lam2*res_HI
                 
                 end do  ! fine fixed-point iteration
                 
              end if

              ! update the outputs with the current values
              E1(i,j,k) = E1new
              E2(i,j,k) = E2new
              E3(i,j,k) = E3new
              ec(i,j,k) = ecnew
              HI(i,j,k) = HInew

           end do
        end do
     end do
     
  !    Hydrogen + Helium chemistry
  else
  
     ! iterate over the domain
     do k=1,Nz
        do j=1,Ny
           do i=1,Nx
              
              ! get shortcut values at this location
              E1new   = E1(i,j,k)
              E2new   = E2(i,j,k)
              E3new   = E3(i,j,k)
              ecnew   = ec(i,j,k)
              HInew   = HI(i,j,k)
              HeInew  = HeI(i,j,k)
              HeIInew = HeII(i,j,k)
              eint = vUn*vUn*(eh(i,j,k)                                 &
                   - KEconst*(vx(i,j,k)**2 + vy(i,j,k)**2 + vz(i,j,k)**2))
              
              ! initialize iteration counters
              l = 0
              l2 = 0
              
              ! perform coarse fixed-point iteration
              do l=0,sweeps

                 if (ERes > FPTol) then
                    ! compute the photo-ionization and photo-heating coeffs
                    call MFSplit_ComputeRadiationIntegrals(piHI(i,j,k),          &
                         piHeI(i,j,k), piHeII(i,j,k), GHI(i,j,k), GHeI(i,j,k),   &
                         GHeII(i,j,k), Ef(i,j,k), E1(i,j,k), E2(i,j,k),          &
                         E3(i,j,k), E1(i,j,k), E2(i,j,k), E3(i,j,k), Nchem,      &
                         ESpectrum, chibar, fsUn, 1, 1, 1, 0, 0, 0, 0, 0, 0, ier)
                 endif

                 ! compute local temperature (based on Model)
                 if (Model == 4) then
                    if (adot == 0.d0) then
                       mol_weight = 0.6d0
                       Temp = (gamma-1.d0)*0.6d0*mp*eint/kb
                       Temp = max(1.d0*Temp,1.d-1)
                    else
                       Temp = ecScale
                    endif
                 else 
                    rhoval = rho(i,j,k)*dUnit
                    nHI  = HI(i,j,k)*nUnit
                    nH  = Hfrac*rhoval/mp
                    nHII = max(1.d0*(nH - nHI), 0.d0)
                    nHe = (1.d0-HFrac)*rhoval/4.d0/mp
                    nHeI = HeI(i,j,k)*nUnit/4.d0
                    nHeII = HeII(i,j,k)*nUnit/4.d0
                    nHeIII = max(1.d0*(nHe - nHeI - nHeII), 0.d0)
                    ne = nHII + 0.25d0*nHeII + 0.5d0*nHeIII
                    mol_weight = rhoval/mp/(0.25d0*(nHeI+nHeII+nHeIII)+nHI+nHII+ne)
                    Temp = gam_1*mol_weight*mp*(ec(i,j,k)*eUn*0.5d0 + eint)/kb
                    Temp = max(1.d0*Temp,1.d-1)
                 endif

                 ! call the local residual routines
                 call MFSplit_AnalyticRadResid(res_E1, res_E2, res_E3, E1new,    &
                      E2new, E3new, HInew, HeInew, HeIInew, E1(i,j,k),           &
                      E2(i,j,k), E3(i,j,k), HI(i,j,k), HeI(i,j,k), HeII(i,j,k),  &
                      dt, dx, src_E1(i,j,k), src_E2(i,j,k), src_E3(i,j,k), a,    &
                      lUnit, nUnit, Nchem, ier)
                 call MFSplit_AnalyticChemResid(res_HI, res_HeI, res_HeII,       &
                      HInew, HeInew, HeIInew, HI(i,j,k), HeI(i,j,k),             &
                      HeII(i,j,k), dt, rho(i,j,k), Temp, src_HI(i,j,k),          &
                      src_HeI(i,j,k), src_HeII(i,j,k), piHI(i,j,k),              &
                      piHeI(i,j,k), piHeII(i,j,k), NTempBins, TStart, TEnd,      &
                      k1Tb, k2Tb, k3Tb, k4Tb, k5Tb, k6Tb, dUnit, nUnit, Nchem,   &
                      HFrac, ier)
                 call MFSplit_AnalyticGasResid(res_ec, ecnew, HInew, HeInew,     &
                      HeIInew, HI(i,j,k), HeI(i,j,k), HeII(i,j,k), dt,           &
                      rho(i,j,k), eint, src_ec(i,j,k), gamma, HFrac, Model, a,   &
                      adot, Comp1, Comp2, CompX, CompT, GHI(i,j,k), GHeI(i,j,k), &
                      GHeII(i,j,k), NTempBins, TStart, TEnd, ceHITb, ceHeITb,    &
                      ceHeIITb, ciHITb, ciHeITb, ciHeISTb, ciHeIITb, reHIITb,    &
                      reHeII1Tb, reHeII2Tb, reHeIIITb, bremTb, dUnit, eUn,       &
                      nUnit, Nchem, ier)
                 
                 ! check the error flag
                 if (ier /= 1)  return
                 
                 ! compute FP change in this iteration and check for convergence
                 ERes = max(abs(res_E1),abs(res_E2))
                 ERes = max(ERes,abs(res_E3))
                 change = max(ERes,abs(res_ec))
                 change = max(change,abs(res_HI))
                 change = max(change,abs(res_HeI))
                 change = max(change,abs(res_HeII))
                 if (change < FPtol)  exit
              
                 ! update the current guesses
                 E1new   = E1new   - lam*res_E1
                 E2new   = E2new   - lam*res_E2
                 E3new   = E3new   - lam*res_E3
                 ecnew   = ecnew   - lam*res_ec
                 HInew   = HInew   - lam*res_HI
                 HeInew  = HeInew  - lam*res_HeI
                 HeIInew = HeIInew - lam*res_HeII
                 
              end do  ! coarse fixed-point iteration
              
              ! if the coarse iteration unsuccessful, repeat with a finer update
              if (change >= FPtol) then
                 
                 ! perform fine fixed-point iteration
                 do l2=0,sweeps2
                    
                    if (ERes > FPTol) then
                       ! compute the photo-ionization and photo-heating coeffs
                       call MFSplit_ComputeRadiationIntegrals(piHI(i,j,k),       &
                            piHeI(i,j,k), piHeII(i,j,k), GHI(i,j,k),             &
                            GHeI(i,j,k), GHeII(i,j,k), Ef(i,j,k), E1(i,j,k),     &
                            E2(i,j,k), E3(i,j,k), E1(i,j,k), E2(i,j,k),          &
                            E3(i,j,k), Nchem, ESpectrum, chibar, fsUn, 1, 1, 1,  &
                            0, 0, 0, 0, 0, 0, ier)
                    endif
                    
                    ! compute local temperature (based on Model)
                    if (Model == 4) then
                       if (adot == 0.d0) then
                          mol_weight = 0.6d0
                          Temp = (gamma-1.d0)*0.6d0*mp*eint/kb
                          Temp = max(1.d0*Temp,1.d-1)
                       else
                          Temp = ecScale
                       endif
                    else 
                       rhoval = rho(i,j,k)*dUnit
                       nHI  = HI(i,j,k)*nUnit
                       nH  = Hfrac*rhoval/mp
                       nHII = max(1.d0*(nH - nHI), 0.d0)
                       nHe = (1.d0-HFrac)*rhoval/4.d0/mp
                       nHeI = HeI(i,j,k)*nUnit/4.d0
                       nHeII = HeII(i,j,k)*nUnit/4.d0
                       nHeIII = max(1.d0*(nHe - nHeI - nHeII), 0.d0)
                       ne = nHII + 0.25d0*nHeII + 0.5d0*nHeIII
                       mol_weight = rhoval/mp/(0.25d0*(nHeI+nHeII+nHeIII)+nHI+nHII+ne)
                       Temp = gam_1*mol_weight*mp*(ec(i,j,k)*eUn*0.5d0 + eint)/kb
                       Temp = max(1.d0*Temp,1.d-1)
                    endif
                    
                    ! call the local residual routines
                    call MFSplit_AnalyticRadResid(res_E1, res_E2, res_E3, E1new, &
                         E2new, E3new, HInew, HeInew, HeIInew, E1(i,j,k),        &
                         E2(i,j,k), E3(i,j,k), HI(i,j,k), HeI(i,j,k),            &
                         HeII(i,j,k), dt, dx, src_E1(i,j,k), src_E2(i,j,k),      &
                         src_E3(i,j,k), a, lUnit, nUnit, Nchem, ier)
                    call MFSplit_AnalyticChemResid(res_HI, res_HeI, res_HeII,    &
                         HInew, HeInew, HeIInew, HI(i,j,k), HeI(i,j,k),          &
                         HeII(i,j,k), dt, rho(i,j,k), Temp, src_HI(i,j,k),       &
                         src_HeI(i,j,k), src_HeII(i,j,k), piHI(i,j,k),           &
                         piHeI(i,j,k), piHeII(i,j,k), NTempBins, TStart, TEnd,   &
                         k1Tb, k2Tb, k3Tb, k4Tb, k5Tb, k6Tb, dUnit, nUnit,       &
                         Nchem, HFrac, ier)
                    call MFSplit_AnalyticGasResid(res_ec, ecnew, HInew, HeInew,  &
                         HeIInew, HI(i,j,k), HeI(i,j,k), HeII(i,j,k), dt,        &
                         rho(i,j,k), eint, src_ec(i,j,k), gamma, HFrac, Model,   &
                         a, adot, Comp1, Comp2, CompX, CompT, GHI(i,j,k),        &
                         GHeI(i,j,k), GHeII(i,j,k), NTempBins, TStart, TEnd,     &
                         ceHITb, ceHeITb, ceHeIITb, ciHITb, ciHeITb, ciHeISTb,   &
                         ciHeIITb, reHIITb, reHeII1Tb, reHeII2Tb, reHeIIITb,     &
                         bremTb, dUnit, eUn, nUnit, Nchem, ier)
                    
                    ! check the error flag
                    if (ier /= 1)  return
                    
                    ! compute FP change in this iteration and check for convergence
                    ERes = max(abs(res_E1),abs(res_E2))
                    ERes = max(ERes,abs(res_E3))
                    change = max(ERes,abs(res_ec))
                    change = max(change,abs(res_HI))
                    change = max(change,abs(res_HeI))
                    change = max(change,abs(res_HeII))
                    if (change < FPtol)  exit
                    
                    ! update the current guesses
                    E1new   = E1new   - lam2*res_E1
                    E2new   = E2new   - lam2*res_E2
                    E3new   = E3new   - lam2*res_E3
                    ecnew   = ecnew   - lam2*res_ec
                    HInew   = HInew   - lam2*res_HI
                    HeInew  = HeInew  - lam2*res_HeI
                    HeIInew = HeIInew - lam2*res_HeII
                    
                 end do  ! fine fixed-point iteration
                 
              end if

              ! update the outputs with the current values
              E1(i,j,k)   = E1new
              E2(i,j,k)   = E2new
              E3(i,j,k)   = E3new
              ec(i,j,k)   = ecnew
              HI(i,j,k)   = HInew
              HeI(i,j,k)  = HeInew
              HeII(i,j,k) = HeIInew
              
           end do
        end do
     end do
     
  end if

  ! exit subroutine
  return

end subroutine MFSplit_AnalyticInitGuess
!=======================================================================






subroutine MFSplit_AnalyticGasResid(ecres, ec, HI, HeI, HeII, HI0, HeI0, &
     HeII0, dt, rho, eint, src_ec, gamma, HFrac, Model, a, adot, Comp1,  &
     Comp2, CompX, CompT, GHI, GHeI, GHeII, NTempBins, TStart, TEnd,     &
     ceHITb, ceHeITb, ceHeIITb, ciHITb, ciHeITb, ciHeISTb, ciHeIITb,     &
     reHIITb, reHeII1Tb, reHeII2Tb, reHeIIITb, bremTb, dUn, eUn, nUn,    &
     Nchem, ier)
  !=======================================================================
  !  written by: Daniel R. Reynolds
  !  date:       August 2009
  !
  !  PURPOSE: 
  !
  !  INPUTS:
  !     ec         - fluid energy correction array
  !     HI         - Hydrogen I number density array
  !     HeI        - Helium I number density array
  !     HeII       - Helium II number density array
  !     HI0        - old Hydrogen I number density array
  !     HeI0       - old Helium I number density array
  !     HeII0      - old Helium II number density array
  !     dt         - time step size
  !     rho        - density array
  !     eint       - fluid internal energy
  !     src_ec     - source function values for gas energy correction eq.
  !     gamma      - constant in ideal gas law
  !     HFrac      - percentage of mass composed of Hydrogen
  !     Model      - flag denoting physical model to use
  !     a          - cosmological expansion parameter
  !     adot       - da/dt
  !     Comp1      - Compton cooling coefficient 1
  !     Comp2      - Compton cooling coefficient 2
  !     CompX      - X-ray Compton heating coefficient
  !     CompT      - X-ray Compton heating temperature 
  !     piH*       - chemistry photo-ionization rates
  !     GH*        - chemistry photo-heating rates
  !     NTempBins  - number of temperature bins for rate tables
  !     TStart     - starting temperature for rate table
  !     TEnd       - ending temperature for rate table
  !     c*Tb, r*Tb - heating/cooling rate tables
  !     *Un        - variable scaling constants
  !     Nchem      - number of chemical species in simulation
  !
  !  OUTPUT ARGUMENTS: 
  !     ecres   - gas energy correction residual
  !     ier     - success/failure flag (0->failure, 1->success)
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
  integer, intent(in)  :: Model, Nchem, NTempBins
  integer, intent(out) :: ier
  REALSUB,  intent(in) :: a, adot
  real, intent(in) :: dt, gamma, HFrac, TStart, TEnd
  real, intent(in) :: Comp1, Comp2, CompX, CompT
  real, intent(in) :: dUn, eUn, nUn
  real, intent(in) :: ec, HI, HeI, HeII, HI0, HeI0, HeII0, rho, eint, &
       src_ec, GHI, GHeI, GHeII
  real, intent(out) :: ecres
  real, intent(in), dimension(NTempBins) :: ceHITb, ceHeITb, ceHeIITb, &
       ciHITb, ciHeITb, ciHeIITb, ciHeISTb, reHIITb, reHeII1Tb, reHeII2Tb, &
       reHeIIITb, bremTb

  !--------------
  ! locals
  integer :: Tidx, Tidxp
  real :: afac, c, kb, mp, lTempS, lTempE, dlTemp
  real :: min_temp, min_ni
  real :: HIval, HeIval, HeIIval, ecval, rhoval
  real :: nH, nHI, nHII, nHe, nHeI, nHeII, nHeIII, ne
  real :: HIanal, HeIanal, HeIIanal, ecanal
  real :: T, lamT, lTemp, Tl, Tr, Tfac
  real :: ceHI, ceHeI, ceHeII, ciHI, ciHeI, ciHeII, ciHeIS
  real :: reHII, reHeII1, reHeII2, reHeIII, brem, G, Lambda
  real :: P2, Q2, k1, k2, k3, k4, k5, k6
  real :: aHI, bHI, cHI, dHI, kHI, expArg, expV, sqD, rt1, rt2
  real :: aHeI, bHeI, cHeI, dHeI, kHeI
  real :: aHeII, bHeII, cHeII, dHeII, kHeII

  !=======================================================================

  ! initialize outputs to have all zero values, flag to success
  ier = 1
  ecres = 0.d0

  ! for isothermal models, just return
  if (Model == 4)  return

  ! initialize constants
  afac = adot/a
  kb = 1.3806504d-16           ! boltzmann constant [erg/K]
  mp = 1.67262171d-24          ! Mass of a proton [g]
  min_temp = 1.d0              ! minimum temperature [K]

  ! lookup table constants
  lTempS = log(TStart)
  lTempE = log(TEnd)
  dlTemp = (lTempE - lTempS)/(1.d0*NTempBins - 1.d0)

  ! get shortcut values time-centered variables
  ecval = ec*0.5d0
  HIval = (HI0 + HI)*0.5d0
  HeIval = (HeI0 + HeI)*0.5d0
  HeIIval = (HeII0 + HeII)*0.5d0
  rhoval = rho*dUn
  nH = Hfrac*rhoval/mp
  nHI = HIval*nUn
  nHII = max(1.d0*(nH - nHI), 0.d0)
  if (Nchem == 1) then
     nHe = 0.d0
     nHeI = 0.d0
     nHeII = 0.d0
     nHeIII = 0.d0
  else
     nHe = (1.d0-Hfrac)*rhoval/mp
     nHeI = HeIval*nUn
     nHeII = HeIIval*nUn
     nHeIII = max(1.d0*(nHe - nHeI - nHeII), 0.d0)
  endif
  ne = nHII + 0.25d0*nHeII + 0.5d0*nHeIII
  
  ! compute temperature and ODE terms
  T = (gamma-1.d0)*rhoval/(0.25d0*(nHeI+nHeII+nHeIII)+nHI+nHII+ne) &
       * (eint+ecval*eUn)/kb
  T = max(1.d0*T,1.d0*min_temp)
  lamT = 3.15614d5/T
  lTemp = min(max(log(T), lTempS), lTempE)
  Tidx = min(NTempBins-1, max(1, int((lTemp-lTempS)/dlTemp)+1))
  Tidxp = Tidx+1
  Tl = lTempS + (Tidx-1)*dlTemp
  Tr = lTempS +  Tidx*dlTemp
  Tfac = (lTemp - Tl)/(Tr - Tl)

  ! compute gas ODE rates
  ceHI = ceHITb(Tidx) + (ceHITb(Tidxp) - ceHITb(Tidx))*Tfac
  ciHI = ciHITb(Tidx) + (ciHITb(Tidxp) - ciHITb(Tidx))*Tfac
  reHII = reHIITb(Tidx) + (reHIITb(Tidxp) - reHIITb(Tidx))*Tfac
  brem = bremTb(Tidx) + (bremTb(Tidxp) - bremTb(Tidx))*Tfac
  Lambda = ne/rhoval*(ceHI*nHI + ciHI*nHI + reHII*nHII        &
       + Comp1*(T-Comp2) + CompX*(T-CompT) + brem*nHII)
  if (Nchem == 3) then
     ceHeI = ceHeITb(Tidx) + (ceHeITb(Tidxp) - ceHeITb(Tidx))*Tfac
     ceHeII = ceHeIITb(Tidx) + (ceHeIITb(Tidxp) - ceHeIITb(Tidx))*Tfac
     ciHeI = ciHeITb(Tidx) + (ciHeITb(Tidxp) - ciHeITb(Tidx))*Tfac
     ciHeII = ciHeIITb(Tidx) + (ciHeIITb(Tidxp) - ciHeIITb(Tidx))*Tfac
     ciHeIS = ciHeISTb(Tidx) + (ciHeISTb(Tidxp) - ciHeISTb(Tidx))*Tfac
     reHeII1 = reHeII1Tb(Tidx) + (reHeII1Tb(Tidxp) - reHeII1Tb(Tidx))*Tfac
     reHeII2 = reHeII2Tb(Tidx) + (reHeII2Tb(Tidxp) - reHeII2Tb(Tidx))*Tfac
     reHeIII = reHeIIITb(Tidx) + (reHeIIITb(Tidxp) - reHeIIITb(Tidx))*Tfac
     Lambda = Lambda + ne/rhoval*(brem*(nHeII/4.d0+nHeIII)       &
          + 0.25d0*(ceHeI*nHeII*ne + ceHeII*nHeII + ciHeI*nHeI   &
          + ciHeII*nHeII + ciHeIS*nHeII*ne + reHeII1*nHeII       &
          + reHeII2*nHeII + reHeIII*nHeIII))
  endif
  G = (nHI*GHI + nHeI*GHeI + nHeII*GHeII)/rhoval
  P2 = 2.d0*afac
  Q2 = (G - Lambda + src_ec)/eUn

  ! compute quasi-steady-state solution for ec, place in ecanal
  if (abs(P2) < 1.0d-14) then
     ecanal = Q2*dt
  else
     if (P2*dt > 7.0d2) then
        ecanal = Q2/P2
     else
        ecanal = (-Q2/P2)*exp(-P2*dt) + Q2/P2
     end if
  end if
  
  ! enforce bounds
  T = (gamma-1.d0)*rhoval/(0.25d0*(nHeI+nHeII+nHeIII)+nHI+nHII+ne) &
       *(eint+ecanal*eUn)/kb
  if (T < min_temp)  &
       ecanal = (min_temp/(gamma-1.d0)/rhoval &
       *(0.25d0*(nHeI+nHeII+nHeIII)+nHI+nHII+ne)*kb-eint)/eUn
  
  ! compute residual
  ecres = ec - ecanal
  
  ! check some things
  !    this statement checks if ecanal = NaN
  if (ecanal /= ecanal) then
     print *,'NaN encountered in AnalyticGasResid (ec)!!'
     ier = 0
  end if
  
  ! exit subroutine
  return

end subroutine MFSplit_AnalyticGasResid
!=======================================================================





subroutine MFSplit_AnalyticChemResid(HIres, HeIres, HeIIres, HI, HeI,   &
     HeII, HI0, HeI0, HeII0, dt, rho, Temp, src_HI, src_HeI, src_HeII,  &
     piHI, piHeI, piHeII, NTempBins, TStart, TEnd, k1Tb, k2Tb, k3Tb,    &
     k4Tb, k5Tb, k6Tb, dUn, nUn, Nchem, HFrac, ier)
  !=======================================================================
  !  written by: Daniel R. Reynolds
  !  date:       August 2009
  !
  !  PURPOSE: 
  !
  !  INPUTS:
  !     HI         - Hydrogen I number density array
  !     HeI        - Helium I number density array
  !     HeII       - Helium II number density array
  !     HI0        - old Hydrogen I number density array
  !     HeI0       - old Helium I number density array
  !     HeII0      - old Helium II number density array
  !     dt         - time step size
  !     rho        - density array
  !     Temp       - temperature
  !     src_HI     - source function values for HI eq.
  !     src_HeI    - source function values for HeI eq.
  !     src_HeII   - source function values for HeII eq.
  !     Model      - flag denoting physical model to use
  !     piH*       - chemistry photo-ionization rates
  !     NTempBins  - number of temperature bins for rate tables
  !     TStart     - starting temperature for rate table
  !     TEnd       - ending temperature for rate table
  !     k*Tb       - chemistry rate tables
  !     c*Tb, r*Tb - heating/cooling rate tables
  !     *Un        - variable scaling constants
  !     Nchem      - number of chemical species in simulation
  !
  !  OUTPUT ARGUMENTS: 
  !     HIres   - Hydrogen I number density residual
  !     HeIres  - Helium I number density residual
  !     HeIIres - Helium II number density residual
  !     ier     - success/failure flag (0->failure, 1->success)
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
  integer, intent(in)  :: Nchem, NTempBins
  integer, intent(out) :: ier
  real, intent(in) :: dt, TStart, TEnd, Temp, HFrac
  real, intent(in) :: dUn, nUn
  real, intent(in) :: HI, HeI, HeII, HI0, HeI0, HeII0, rho,      &
       src_HI, src_HeI, src_HeII, piHI, piHeI, piHeII
  real, intent(out) :: HIres, HeIres, HeIIres
  real, intent(in), dimension(NTempBins) :: k1Tb, k2Tb, k3Tb,    &
       k4Tb, k5Tb, k6Tb

  !--------------
  ! locals
  integer :: Tidx, Tidxp
  real :: mp, lTempS, lTempE, dlTemp, min_ni, HIval, HeIval, HeIIval
  real :: rhoval, nH, nHI, nHII, nHe, nHeI, nHeII, nHeIII, ne
  real :: HIanal, HeIanal, HeIIanal, lamT, lTemp, Tl, Tr, Tfac
  real :: k1, k2, k3, k4, k5, k6
  real :: aHI, bHI, cHI, dHI, kHI, expArg, expV, sqD, rt1, rt2
  real :: aHeI, bHeI, cHeI, dHeI, kHeI
  real :: aHeII, bHeII, cHeII, dHeII, kHeII

  !=======================================================================

  ! initialize outputs to have all zero values, flag to success
  ier = 1
  HIres = 0.d0
  if (Nchem == 3) then
     HeIres = 0.d0
     HeIIres = 0.d0
  end if
  HIanal = 0.d0
  HeIanal = 0.d0
  HeIIanal = 0.d0

  ! initialize constants
  mp = 1.67262171d-24          ! Mass of a proton [g]
  min_ni   = 0.d0              ! minimum chem number density [cm^{-3}]

  ! lookup table constants
  lTempS = log(TStart)
  lTempE = log(TEnd)
  dlTemp = (lTempE - lTempS)/(1.d0*NTempBins - 1.d0)


  ! get shortcut values time-centered variables
  HIval = (HI0 + HI)*0.5d0
  HeIval = (HeI0 + HeI)*0.5d0
  HeIIval = (HeII0 + HeII)*0.5d0
  rhoval = rho*dUn
  nH = Hfrac*rhoval/mp
  nHI = HIval*nUn
  nHII = max(1.d0*(nH - nHI), 0.d0)
  if (Nchem == 1) then
     nHe = 0.d0
     nHeI = 0.d0
     nHeII = 0.d0
     nHeIII = 0.d0
  else
     nHe = (1.d0-Hfrac)*rhoval/mp
     nHeI = HeIval*nUn
     nHeII = HeIIval*nUn
     nHeIII = max(1.d0*(nHe - nHeI - nHeII), 0.d0)
  endif
  ne = nHII + 0.25d0*nHeII + 0.5d0*nHeIII
  
  ! compute temperature-based ODE terms
  lamT = 3.15614d5/Temp
  lTemp = min(max(log(Temp), lTempS), lTempE)
  Tidx = min(NTempBins-1, max(1, int((lTemp-lTempS)/dlTemp)+1))
  Tidxp = Tidx+1
  Tl = lTempS + (Tidx-1)*dlTemp
  Tr = lTempS +  Tidx*dlTemp
  Tfac = (lTemp - Tl)/(Tr - Tl)

  ! compute chemistry ODE rates
  k1 = k1Tb(Tidx) + (k1Tb(Tidxp) - k1Tb(Tidx))*Tfac
  k2 = 2.753d-14*lamT**(1.5d0) *                 &
       (1.d0+(lamT/2.74d0)**(0.407d0))**(-2.242d0)
  aHI = (k1+k2)*nUn
  bHI = -(k1 + 2.d0*k2)*nH - (k1+k2)*(nHeII + 2.d0*nHeIII)/4.d0 - piHI
  cHI = (k2*nH*nH + k2*nH*(nHeII + 2.d0*nHeIII)/4.d0 + src_HI)/nUn
  if (Nchem == 3) then
     k3 = k3Tb(Tidx) + (k3Tb(Tidxp) - k3Tb(Tidx))*Tfac
     k4 = k4Tb(Tidx) + (k4Tb(Tidxp) - k4Tb(Tidx))*Tfac
     k5 = k5Tb(Tidx) + (k5Tb(Tidxp) - k5Tb(Tidx))*Tfac
     k6 = k6Tb(Tidx) + (k6Tb(Tidxp) - k6Tb(Tidx))*Tfac
     aHeI = k3*0.5d0*nUn
     bHeI = -piHeI - k3*(nHII + 0.5d0*nHe - 0.25d0*nHeII) - k4*0.5d0*nHeII
     cHeI = (k4*nHeII*(nHII + 0.5d0*nHe - 0.25d0*nHeII) + src_HeI)/nUn
     aHeII = 0.25d0*(k4+k5+k6)*nUn
     bHeII = -piHeII - (k4+k5+k6)*(nHII + 0.5d0*(nHe-nHeI)) &
          - 0.25d0*(k3*nHeI + k6*(nHe-nHeI))
     cHeII = ((k3*nHeI + k6*(nHe-nHeI))*(nHII + 0.5d0*(nHe-nHeI)) &
          + nHeI*piHeI + src_HeII)/nUn
  endif

  ! compute quasi-steady-state solution for HI, place in HIres
  dHI = bHI*bHI - 4.d0*aHI*cHI
  
  !    if the second-order terms are very small, treat it as a linear ODE
  if (abs(aHI*HIval**2/(cHI-bHI*HIval)) < 1.0d-8) then
     if (bHI*dt < -7.0d2) then
        HIanal = -cHI/bHI
     else
        HIanal = (HI0 + cHI/bHI)*exp(bHI*dt) - cHI/bHI
     end if
     
     !    otherwise, go ahead with the quadratic ODE solution
  else
     
     ! no roots to quadratic
     if (dHI/bHI/bHI < -1.0d-8) then
        sqD = sqrt(-dHI)
        kHI = 2.d0/sqD*atan((2.d0*aHI*HI0+bHI)/sqD)
        HIanal = (sqD*tan((dt+kHI)*0.5d0*sqD)-bHI)*0.5d0/aHI
        
        ! two roots to quadratic
     elseif (dHI/bHI/bHI > 1.0d-8) then
        rt1 = (-bHI+sqrt(dHI))*0.5d0/aHI
        rt2 = (-bHI-sqrt(dHI))*0.5d0/aHI 
        kHI = log(abs((HI0-rt2)/(HI0-rt1)))/aHI/(rt2-rt1)
        expArg = aHI*(dt+kHI)*(rt2-rt1)
        if (expArg < -7.0d2) then
           HIanal = rt2
        elseif (expArg > 7.0d2) then
           HIanal = rt1
        else
           expV = exp(expArg)
           if ((HI0-rt1)*(HI0-rt2) > 0.d0) then
              HIanal = (rt2-rt1*expV)/(1.d0-expV)
           else
              HIanal = (rt2+rt1*expV)/(1.d0+expV)
           end if
        end if
        
        ! one double root to quadratic
     else
        rt1 = -bHI*0.5d0/aHI
        kHI = 1.d0/aHI/(rt1-HI0)
        HIanal = rt1 - 1.d0/aHI/(dt+kHI)
        
     end if  ! roots
  end if  ! quadratic


  ! Helium chemistry
  if (Nchem == 3) then

     ! compute quasi-steady-state solution for HeI, place in HeIres
     dHeI = bHeI*bHeI - 4.d0*aHeI*cHeI

     !    if the second-order terms are very small, treat it as a linear ODE
     if (abs(aHeI*HeIval**2/(cHeI-bHeI*HeIval)) < 1.0d-8) then
        if (bHeI*dt < -7.0d2) then
           HeIanal = -cHeI/bHeI
        else
           HeIanal = (HeI0 + cHeI/bHeI)*exp(bHeI*dt) - cHeI/bHeI
        end if

     !    otherwise, go ahead with the quadratic ODE solution
     else

        ! no roots to quadratic
        if (dHeI/bHeI/bHeI < -1.0d-8) then
           sqD = sqrt(-dHeI)
           kHeI = 2.d0/sqD*atan((2.d0*aHeI*HeI0+bHeI)/sqD)
           HeIanal = (sqD*tan((dt+kHeI)*0.5d0*sqD)-bHeI)*0.5d0/aHeI

        ! two roots to quadratic
        elseif (dHeI/bHeI/bHeI > 1.0d-8) then
           rt1 = (-bHeI+sqrt(dHeI))*0.5d0/aHeI
           rt2 = (-bHeI-sqrt(dHeI))*0.5d0/aHeI 
           kHeI = log(abs((HeI0-rt2)/(HeI0-rt1)))/aHeI/(rt2-rt1)
           expArg = aHeI*(dt+kHeI)*(rt2-rt1)
           if (expArg < -7.0d2) then
              HeIanal = rt2
           elseif (expArg > 7.0d2) then
              HeIanal = rt1
           else
              expV = exp(expArg)
              if ((HeI0-rt1)*(HeI0-rt2) > 0.d0) then
                 HeIanal = (rt2-rt1*expV)/(1.d0-expV)
              else
                 HeIanal = (rt2+rt1*expV)/(1.d0+expV)
              end if
           end if

        ! one double root to quadratic
        else
           rt1 = -bHeI*0.5d0/aHeI
           kHeI = 1.d0/aHeI/(rt1-HeI0)
           HeIanal = rt1 - 1.d0/aHeI/(dt+kHeI)

        end if  ! roots
     end if  ! quadratic


     ! compute quasi-steady-state solution for HeII, place in HeIIres
     dHeII = bHeII*bHeII - 4.d0*aHeII*cHeII

     !    if the second-order terms are very small, treat it as a linear ODE
     if (abs(aHeII*HeIIval**2/(cHeII-bHeII*HeIIval)) < 1.0d-8) then
        if (bHeII*dt < -7.0d2) then
           HeIIanal = -cHeII/bHeII
        else
           HeIIanal = (HeII0 + cHeII/bHeII)*exp(bHeII*dt) - cHeII/bHeII
        end if

     !    otherwise, go ahead with the quadratic ODE solution
     else

        ! no roots to quadratic
        if (dHeII/bHeII/bHeII < -1.0d-8) then
           sqD = sqrt(-dHeII)
           kHeII = 2.d0/sqD*atan((2.d0*aHeII*HeII0+bHeII)/sqD)
           HeIIanal = (sqD*tan((dt+kHeII)*0.5d0*sqD)-bHeII)*0.5d0/aHeII

        ! two roots to quadratic
        elseif (dHeII/bHeII/bHeII > 1.0d-8) then
           rt1 = (-bHeII+sqrt(dHeII))*0.5d0/aHeII
           rt2 = (-bHeII-sqrt(dHeII))*0.5d0/aHeII 
           kHeII = log(abs((HeII0-rt2)/(HeII0-rt1)))/aHeII/(rt2-rt1)
           expArg = aHeII*(dt+kHeII)*(rt2-rt1)
           if (expArg < -7.0d2) then
              HeIIanal = rt2
           elseif (expArg > 7.0d2) then
              HeIIanal = rt1
           else
              expV = exp(expArg)
              if ((HeII0-rt1)*(HeII0-rt2) > 0.d0) then
                 HeIIanal = (rt2-rt1*expV)/(1.d0-expV)
              else
                 HeIIanal = (rt2+rt1*expV)/(1.d0+expV)
              end if
           end if

        ! one double root to quadratic
        else
           rt1 = -bHeII*0.5d0/aHeII
           kHeII = 1.d0/aHeII/(rt1-HeII0)
           HeIIanal = rt1 - 1.d0/aHeII/(dt+kHeII)

        end if  ! roots
     end if  ! quadratic

  endif  ! Nchem

  ! enforce bounds
  HIanal = max(min(HIanal, nH/nUn), min_ni/nUn)
  HeIanal = max(min(HeIanal, nHe/nUn), min_ni/nUn)
  HeIIanal = max(min(HeIIanal, nHe/nUn), min_ni/nUn)
  
  
  ! compute residuals
  HIres = HI - HIanal
  HeIres = HeI - HeIanal
  HeIIres = HeII - HeIIanal
  
  ! check some things
  !    this statement checks if HIanal = NaN
  if (HIanal /= HIanal) then
     print *,'NaN encountered in AnalyticChemResid (HI)!!'
     ier = 0
  end if
  
  !    this statement checks if HeIanal = NaN
  if (HeIanal /= HeIanal) then
     print *,'NaN encountered in AnalyticChemResid (HeI)!!'
     ier = 0
  end if
  
  !    this statement checks if HeIIanal = NaN
  if (HeIIanal /= HeIIanal) then
     print *,'NaN encountered in AnalyticChemResid (HeII)!!'
     ier = 0
  end if
  
  ! exit subroutine
  return

end subroutine MFSplit_AnalyticChemResid
!=======================================================================





subroutine MFSplit_AnalyticRadResid(E1res, E2res, E3res, E1, E2, E3, HI, &
     HeI, HeII, E10, E20, E30, HI0, HeI0, HeII0, dt, dx, src_E1, src_E2, &
     src_E3, a, lUn, nUn, Nchem, ier)
  !=======================================================================
  !  written by: Daniel R. Reynolds
  !  date:       August 2009
  !
  !  PURPOSE: 
  !
  !  INPUTS:
  !     E1         - radiation energy freq 1 array
  !     E2         - radiation energy freq 2 array
  !     E3         - radiation energy freq 3 array
  !     HI         - Hydrogen I number density array
  !     HeI        - Helium I number density array
  !     HeII       - Helium II number density array
  !     E10        - old radiation energy freq 1 array
  !     E20        - old radiation energy freq 2 array
  !     E30        - old radiation energy freq 3 array
  !     HI0        - old Hydrogen I number density array
  !     HeI0       - old Helium I number density array
  !     HeII0      - old Helium II number density array
  !     dt         - time step size
  !     dx         - spatial mesh size
  !     rho        - density array
  !     src_E1     - source function values for radiation equation freq 1
  !     src_E2     - source function values for radiation equation freq 2
  !     src_E3     - source function values for radiation equation freq 3
  !     a          - cosmological expansion parameter (new & old)
  !     *Un        - variable scaling constants
  !     Nchem      - number of chemical species in simulation
  !
  !  OUTPUT ARGUMENTS: 
  !     E1res   - radiation energy density freq 1 residual
  !     E2res   - radiation energy density freq 2 residual
  !     E3res   - radiation energy density freq 3 residual
  !     ier     - success/failure flag (0->failure, 1->success)
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
  integer, intent(in)  :: Nchem
  integer, intent(out) :: ier
  REALSUB,  intent(in) :: a
  real, intent(in) :: dt, dx
  real, intent(in) :: lUn, nUn
  real, intent(in) :: E1, E2, E3, HI, HeI, HeII, E10, E20, E30, &
       HI0, HeI0, HeII0, src_E1, src_E2, src_E3
  real, intent(out) :: E1res, E2res, E3res

  !--------------
  ! locals
  real :: ev2erg, c, hp, nu0_HI, nu0_HeI, nu0_HeII, min_rad
  real :: HIval, HeIval, HeIIval, E1val, E2val, E3val
  real :: E1anal, E2anal, E3anal, kE, P, Q, cond, dx_sc
  real :: wts(7), taus(7), vals(7), fvals(7), ival
  real :: sHI, sHeI, sHeII

  !=======================================================================

  ! initialize outputs to have all zero values, flag to success
  ier = 1
  E1res = 0.d0
  E2res = 0.d0
  E3res = 0.d0

  ! initialize constants
  dx_sc = dx*lUn/a
  ev2erg = 1.60217653e-12      ! conversion constant from eV to ergs
  c  = 2.99792458d10           ! speed of light [cm/s]
  hp = 6.6260693d-27           ! Planck's constant [ergs*s]
  nu0_HI   = 13.6d0*ev2erg/hp  ! ionization frequency of HI   [hz]
  nu0_HeI  = 24.6d0*ev2erg/hp  ! ionization frequency of HeI  [hz]
  nu0_HeII = 54.4d0*ev2erg/hp  ! ionization frequency of HeII [hz]
  min_rad  = 0.d0              ! minimum radiation density [ergs/cm^3]

!!$  ! radiation time integral weights and integration nodes
!!$  wts  = (/ 1.d0, 4.d0, 2.d0, 4.d0, 2.d0, 4.d0, 1.d0   /)*dt/18.d0
!!$  taus = (/ 6.d0, 5.d0, 4.d0, 3.d0, 2.d0, 1.d0, 1.d-11 /)*dt/6.d0

  ! get shortcut values time-centered variables
  E1val = (E10 + E1)*0.5d0
  E2val = (E20 + E2)*0.5d0
  E3val = (E30 + E3)*0.5d0
  HIval = (HI0 + HI)*nUn*0.5d0
  if (Nchem == 1) then
     HeIval = 0.d0
     HeIIval = 0.d0
  else
     HeIval = (HeI0 + HeI)*nUn*0.5d0
     HeIIval = (HeII0 + HeII)*nUn*0.5d0
  endif

  ! compute radiation eq 1 ODE rates
  call MFSplit_HICrossSection(sHI, nu0_HI, 1)
  kE = sHI*HIval
  P = c*kE
  Q = src_E1
  cond = 16.d0*c/kE/3.d0

  ! compute quasi-steady-state solution for E1, place in E1anal
  if (abs(P) < 1.0d-14) then
     E1anal = E10 + Q*dt
  else
     if (P*dt > 7.0d2) then
        E1anal = Q/P
     else
        E1anal = (E10 - Q/P)*exp(-P*dt) + Q/P
     end if
  end if
!!$     ! The following analytical solution is derived using the Green's function
!!$     ! solution for the forced, constant-coefficient heat equation, where we 
!!$     ! only consider the spatial domain of dependence to be this cell (and 
!!$     ! ignore flow from farther distances).  While this is only a crude 
!!$     ! approximation to the solution, it is nonetheless very fast and provides
!!$     ! the exact analytical solution for the cell containing an emissivity 
!!$     ! source; as this is where the local nonlinear problem is hardest to solve, 
!!$     ! it should prove quite beneficial.
!!$     !    perform the time integral only if source is nonzero
!!$     E1anal = 0.d0
!!$     if (abs(Q) > 1.d-13) then
!!$        ! evaluate the Green's function integrand at tau values, 
!!$        ! and combine to form the time-integral
!!$        vals = dx_sc/sqrt(cond*taus)
!!$        call erf_vec(fvals, vals, 7)
!!$        ival = sum(wts * exp(-P*taus) * fvals**3)
!!$
!!$        ! compute contribution due to sources
!!$        E1anal = ival*Q
!!$!!$     else
!!$!!$        vals(1) = dx_sc/sqrt(cond*dt)
!!$!!$        call erf_vec(fvals, vals, 1)
!!$     end if
!!$     !    compute contribution due to previous time solution
!!$!!$     !    (re-use fvals since it was last called with tau=0)
!!$!!$     E1anal = E1anal + E10*exp(-P*dt)*fvals(1)**3
!!$     E1anal = E1anal + E10*exp(-P*dt)
     

  ! compute radiation eq 2 ODE rates
  call MFSplit_HICrossSection(sHI, nu0_HeI, 1)
  call MFSplit_HICrossSection(sHeI, nu0_HeI, 1)
  kE = sHI*HIval + sHeI*HeIval
  P = c*kE
  Q = src_E2
  cond = 16.d0*c/kE/3.d0

  ! compute quasi-steady-state solution for E1, place in E1anal
  if (abs(P) < 1.0d-14) then
     E2anal = E20 + Q*dt
  else
     if (P*dt > 7.0d2) then
        E2anal = Q/P
     else
        E2anal = (E20 - Q/P)*exp(-P*dt) + Q/P
     end if
  end if
!!$     ! The following analytical solution is derived using the Green's function
!!$     ! solution for the forced, constant-coefficient heat equation, where we 
!!$     ! only consider the spatial domain of dependence to be this cell (and 
!!$     ! ignore flow from farther distances).  While this is only a crude 
!!$     ! approximation to the solution, it is nonetheless very fast and provides
!!$     ! the exact analytical solution for the cell containing an emissivity 
!!$     ! source; as this is where the local nonlinear problem is hardest to solve, 
!!$     ! it should prove quite beneficial.
!!$     !    perform the time integral only if source is nonzero
!!$     E2anal = 0.d0
!!$     if (abs(Q) > 1.d-13) then
!!$        ! evaluate the Green's function integrand at tau values, 
!!$        ! and combine to form the time-integral
!!$        vals = dx_sc/sqrt(cond*taus)
!!$        call erf_vec(fvals, vals, 7)
!!$        ival = sum(wts * exp(-P*taus) * fvals**3)
!!$
!!$        ! compute contribution due to sources
!!$        E2anal = ival*Q
!!$!!$     else
!!$!!$        vals(1) = dx_sc/sqrt(cond*dt)
!!$!!$        call erf_vec(fvals, vals, 1)
!!$     end if
!!$     !    compute contribution due to previous time solution
!!$!!$     !    (re-use fvals since it was last called with tau=0)
!!$!!$     E2anal = E2anal + E20*exp(-P*dt)*fvals(1)**3
!!$     E2anal = E2anal + E20*exp(-P*dt)
     

  ! compute radiation eq 3 ODE rates
  call MFSplit_HICrossSection(sHI, nu0_HeII, 1)
  call MFSplit_HICrossSection(sHeI, nu0_HeII, 1)
  call MFSplit_HICrossSection(sHeII, nu0_HeII, 1)
  kE = sHI*HIval + sHeI*HeIval + sHeII*HeIIval
  P = c*kE
  Q = src_E3
  cond = 16.d0*c/kE/3.d0

  ! compute quasi-steady-state solution for E3, place in E3anal
  if (abs(P) < 1.0d-14) then
     E3anal = E30 + Q*dt
  else
     if (P*dt > 7.0d2) then
        E3anal = Q/P
     else
        E3anal = (E30 - Q/P)*exp(-P*dt) + Q/P
     end if
  end if
!!$     ! The following analytical solution is derived using the Green's function
!!$     ! solution for the forced, constant-coefficient heat equation, where we 
!!$     ! only consider the spatial domain of dependence to be this cell (and 
!!$     ! ignore flow from farther distances).  While this is only a crude 
!!$     ! approximation to the solution, it is nonetheless very fast and provides
!!$     ! the exact analytical solution for the cell containing an emissivity 
!!$     ! source; as this is where the local nonlinear problem is hardest to solve, 
!!$     ! it should prove quite beneficial.
!!$     !    perform the time integral only if source is nonzero
!!$     E3anal = 0.d0
!!$     if (abs(Q) > 1.d-13) then
!!$        ! evaluate the Green's function integrand at tau values, 
!!$        ! and combine to form the time-integral
!!$        vals = dx_sc/sqrt(cond*taus)
!!$        call erf_vec(fvals, vals, 7)
!!$        ival = sum(wts * exp(-P*taus) * fvals**3)
!!$
!!$        ! compute contribution due to sources
!!$        E3anal = ival*Q
!!$!!$     else
!!$!!$        vals(1) = dx_sc/sqrt(cond*dt)
!!$!!$        call erf_vec(fvals, vals, 1)
!!$     end if
!!$     !    compute contribution due to previous time solution
!!$!!$     !    (re-use fvals since it was last called with tau=0)
!!$!!$     E3anal = E3anal + E30*exp(-P*dt)*fvals(1)**3
!!$     E3anal = E3anal + E30*exp(-P*dt)
     

  ! enforce bounds
  E1anal = max(E1anal, min_rad)
  E2anal = max(E2anal, min_rad)
  E3anal = max(E3anal, min_rad)
  
  ! compute residuals
  E1res = E1 - E1anal
  E2res = E2 - E2anal
  E3res = E3 - E3anal
  
  ! check some things
  !    this statement checks if E1anal = NaN
  if (E1anal /= E1anal) then
     print *,'NaN encountered in AnalyticRadResid (E1)!!'
     ier = 0
     return
  end if
  !    this statement checks if E2anal = NaN
  if (E2anal /= E2anal) then
     print *,'NaN encountered in AnalyticRadResid (E2)!!'
     ier = 0
     return
  end if
  !    this statement checks if E3anal = NaN
  if (E3anal /= E3anal) then
     print *,'NaN encountered in AnalyticRadResid (E3)!!'
     ier = 0
     return
  end if

  ! exit subroutine
  return

end subroutine MFSplit_AnalyticRadResid
!=======================================================================
