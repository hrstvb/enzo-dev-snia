!=======================================================================
!
! Copyright 2009 Daniel R. Reynolds
!
! This software is released under the terms of the "Enzo Public License"
! in the accompanying LICENSE file.
!
!=======================================================================
subroutine MFProb_RadJac(JE1_E1, JE1_HI, JE2_E2, JE2_HI, JE2_HeI,      &
     JE3_E3, JE3_HI, JE3_HeI, JE3_HeII, E1, E2, E3, HI, HeI, HeII, dt, &
     theta, a, adot, aUn, lUn, rUn, nUn, Nchem, Nx, Ny, Nz, NGxl,      &
     NGxr, NGyl, NGyr, NGzl, NGzr, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       August 2009
!  modified1:  
!
!  PURPOSE: Computes the spatially-local jacobian components from the 
!           nonlinear residual for the multi-frequency radiation problem,
!              E1dot = Div(D(E1)*Grad(E1)) + 4*pi*src1 - c*k1*E1,
!              E2dot = Div(D(E2)*Grad(E2)) + 4*pi*src2 - c*k2*E2,
!              E3dot = Div(D(E3)*Grad(E3)) + 4*pi*src3 - c*k3*E3,
!           where D(E) is a nonlinear flux-limiter depending on E.  
!
!  INPUTS:
!     E*         - Radiation energy density at each frequency (new & old times)
!     HI         - Hydrogen I density (new & old times)
!     HeI        - Helium I density (new & old times)
!     HeII       - Helium II density (new & old times)
!     dt         - time step size
!     theta      - time integration parameter
!     a,a0       - cosmological expansion parameter (new & old times)
!     adot,adot0 - da/dt (new & old times)
!     *Un,*Un0   - variable scaling constants (new & old times)
!     Nx,Ny,Nz   - active mesh size in each direction
!     NG*l/NG*r  - left/right ghost cells in each direction
!
!     Note: the vector inputs are of size (Nx + NGxl + NGxr) in 
!     the x-direction, others are similar.
!
!  OUTPUT ARGUMENTS: 
!     JE1_*      - E1 eqn jacobians wrt inputs
!     JE2_*      - E2 eqn jacobians wrt inputs
!     JE3_*      - E3 eqn jacobians wrt inputs
!     ier        - success/failure output flag (0->failure, 1->success)
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
  integer, intent(in)  :: Nx, NGxl, NGxr, Nchem
  integer, intent(in)  :: Ny, NGyl, NGyr
  integer, intent(in)  :: Nz, NGzl, NGzr
  integer, intent(out) :: ier
  REALSUB, intent(in)  :: a, adot
  real, intent(in) :: aUn, lUn, rUn, nUn
  real, intent(in) :: dt, theta
  real, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr), &
       intent(in) :: E1, E2, E3
  real, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr), &
       intent(in) :: HI, HeI, HeII
  real, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr), &
       intent(out) :: JE1_E1, JE1_HI, JE2_E2, JE2_HI, JE2_HeI,       &
                      JE3_E3, JE3_HI, JE3_HeI, JE3_HeII

  !--------------
  ! local declarations
  real :: c, hp, dtfac, ev2erg, nu0_HI, nu0_HeI, nu0_HeII, kap
  real :: sHI1, sHI2, sHeI2, sHI3, sHeI3, sHeII3

  !=======================================================================

  ! set output flag
  ier = 1

  ! set shortcut values 
  c = 2.99792458d10             ! speed of light [cm/s]
  hp = 6.6260693d-27            ! Planck's constant (ergs*s)
  dtfac  = dt*theta
  ev2erg = 1.60217653d-12       ! conversion constant from eV to ergs
  nu0_HI   = 13.6d0*ev2erg/hp   ! ionization threshold of HI (hz)
  nu0_HeI  = 24.6d0*ev2erg/hp   ! ionization threshold of HeI (hz)
  nu0_HeII = 54.4d0*ev2erg/hp   ! ionization threshold of HeII (hz)

  !   set shortcut variables
  call MFProb_HICrossSection(sHI1, nu0_HI, 1)
  call MFProb_HICrossSection(sHI2, nu0_HeI, 1)
  call MFProb_HeICrossSection(sHeI2, nu0_HeI, 1)
  call MFProb_HICrossSection(sHI3, nu0_HeII, 1)
  call MFProb_HeICrossSection(sHeI3, nu0_HeII, 1)
  call MFProb_HeIICrossSection(sHeII3, nu0_HeII, 1)


  ! compute derivatives depending on chemistry
  if (Nchem == 1) then

     ! compute derivatives over domain
     !    Eloc = (1.d0 + dtfac*4.d0*kap)*E(i,j,k)
     !     kap = sHI*HI(i,j,k)*nUn
     JE1_E1 = 1.d0 + dtfac*4.d0*sHI1*HI*nUn
     JE2_E2 = 1.d0 + dtfac*4.d0*sHI2*HI*nUn
     JE3_E3 = 1.d0 + dtfac*4.d0*sHI3*HI*nUn
     JE1_HI = dtfac*4.d0*sHI1*nUn*E1
     JE2_HI = dtfac*4.d0*sHI2*nUn*E2
     JE3_HI = dtfac*4.d0*sHI3*nUn*E3
     
  else

     ! compute derivatives over domain
     !    Eloc = (1.d0 + dtfac*4.d0*kap)*E(i,j,k)
     !     kap = (sHI*HI(i,j,k) + sHeI*HeI(i,j,k) + sHeII*HeII(i,j,k))*nUn
     JE1_E1   = 1.d0 + dtfac*4.d0*sHI1*HI*nUn
     JE2_E2   = 1.d0 + dtfac*4.d0*(sHI2*HI+sHeI2*HeI)*nUn
     JE3_E3   = 1.d0 + dtfac*4.d0*(sHI3*HI+sHeI3*HeI+sHeII3*HeII)*nUn
     JE1_HI   = dtfac*4.d0*sHI1*nUn*E1
     JE2_HI   = dtfac*4.d0*sHI2*nUn*E2
     JE2_HeI  = dtfac*4.d0*sHeI2*nUn*E2
     JE3_HI   = dtfac*4.d0*sHI3*nUn*E3
     JE3_HeI  = dtfac*4.d0*sHeI3*nUn*E3
     JE3_HeII = dtfac*4.d0*sHeII3*nUn*E3
              
  endif

  return
end subroutine MFProb_RadJac
!=======================================================================
