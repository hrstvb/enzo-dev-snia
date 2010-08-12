!=======================================================================
!
! Copyright 2009 Daniel R. Reynolds
!
! This software is released under the terms of the "Enzo Public License"
! in the accompanying LICENSE file.
!
!=======================================================================
subroutine MFProb_Sources(eta1, eta2, eta3, NGamDot_FS, ecsrc, HIsrc,   &
     HeIsrc, HeIIsrc, time, a, Model, ProbType, ESpectrum, Nchem,       &
     SMEmiss, NGammaDot, EtaRadius, EtaCenter, aUn, dUn, vUn, lUn, tUn, &
     rUn, eUn, nUn, Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr, NGzl, NGzr,     &
     x0L, x0R, x1L, x1R, x2L, x2R, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       August 2009
!
!  PURPOSE: Computes the sources into the radiation energy, gas energy, 
!           and species number density equations.
!
!  INPUTS:
!     time       - simulation time for evaluation
!     a          - cosmological expansion parameter
!     Model      - flag denoting physical model to use
!     ProbType   - flag denoting physical problem to run
!     Nchem      - number of chemical species
!     SMEmiss    - percentage of mass composed of Hydrogen
!     NGammaDot  - ionization source
!     EtaRadius  - ionization source radius in cells
!     EtaCenter  - ionization source center (comoving, 3D coordinates in cm)
!     *Un        - variable scaling constants
!     Nx,Ny,Nz   - active mesh size in each direction
!     NG*l/NG*r  - left/right ghost cells in each direction
!     x*L/x*R    - left/right subdomain boundaries (comoving, no ghosts)
!
!     Note: the vector inputs are of size (Nx + NGxl + NGxr) in 
!           the x-direction, others are similar.
!
!     Note: for now we only supply emissivity sources, and only when 
!           SMEmiss is disabled.
!
!  OUTPUT ARGUMENTS: 
!     eta1       - array of emissivity sources for rad freq. 1
!     eta2       - array of emissivity sources for rad freq. 2
!     eta3       - array of emissivity sources for rad freq. 3
!     NGamDot_FS - ionization source strength for free-streaming radiation
!     ecsrc      - array of gas energy sources
!     HIsrc      - array of HI number density sources
!     HeIsrc     - array of HeI number density sources
!     HeIIsrc    - array of HeII number density sources
!     ier        - success/failure flag (1->success, 0->failure)
!
!  EXTERNALS: 
!
!  LOCALS:
!
!=======================================================================
  implicit none
#include "fortran.def"

  !--------------
  ! argument declarations
  integer, intent(in) :: Nchem, Model, ProbType, SMEmiss, ESpectrum
  integer, intent(in) :: Nx, NGxl, NGxr
  integer, intent(in) :: Ny, NGyl, NGyr
  integer, intent(in) :: Nz, NGzl, NGzr
  integer, intent(out) :: ier
  REALSUB,  intent(in) :: a
  real, intent(in) :: time, NGammaDot, EtaRadius
  real, intent(in) :: EtaCenter(3)
  real, intent(in) :: aUn, dUn, vUn, lUn, tUn, rUn, eUn, nUn
  real, intent(in) :: x0L, x0R, x1L, x1R, x2L, x2R
  real, intent(out) :: NGamDot_FS
  real, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr), &
       intent(out) :: eta1, eta2, eta3, ecsrc, HIsrc, HeIsrc, HeIIsrc
  
  !--------------
  ! locals
  integer :: i, j, k, l, nbins
  integer, parameter :: nnus=99
  real :: pi, h_nu1, h_nu2, h_nu3, etaconst1, etaconst2, etaconst3
  real :: dx, dy, dz, dV, cellXl, cellXr, cellYl, cellYr, cellZl, cellZr
  real :: cellXc, cellYc, cellZc, ev2erg
  real :: hp, nu0_HI, nu0_HeI, nu0_HeII, chi1, chi2, chi3, chiint, chinuint
  real :: Llimit, Ulimit
  real, dimension(nnus) :: nus, etas, nusB, chis, chisB, wts, wtsB

  !=======================================================================

  ! initialize outputs to have all zero values, flag to success
  !    only overwrite emissivity sources if starmaker is not being used
  ier = 1
  ecsrc = 0.d0
  HIsrc = 0.d0
  if (Nchem == 3) then
     HeIsrc  = 0.d0
     HeIIsrc = 0.d0
  endif
#ifdef EMISSIVITY  
  if (SMEmiss > 0)  return
#endif
  eta1 = 0.d0
  eta2 = 0.d0
  eta3 = 0.d0

  
  ! initialize constants
  pi    = 4.D0*datan(1.D0)
  dx    = (x0R-x0L)/Nx                ! mesh spacing (comoving), x0 direction
  dy    = (x1R-x1L)/Ny                ! mesh spacing (comoving), x1 direction
  dz    = (x2R-x2L)/Nz                ! mesh spacing (comoving), x2 direction
  dV    = dx*dy*dz*(lUn)**3           ! cell volume (proper)
  ev2erg = 1.60217653e-12             ! coefficient to convert eV to ergs
  h_nu1 = 13.6d0*ev2erg               ! ionization energy of HI [ergs]
  h_nu2 = 24.6d0*ev2erg               ! ionization energy of HeI [ergs]
  h_nu3 = 54.4d0*ev2erg               ! ionization energy of HeII [ergs]
  hp = 6.6260693d-27                  ! Planck's constant (ergs*s)
  nu0_HI   = h_nu1/hp                 ! ionization threshold of HI (hz)
  nu0_HeI  = h_nu2/hp                 ! ionization threshold of HeI (hz)
  nu0_HeII = h_nu3/hp                 ! ionization threshold of HeII (hz)


  ! get the radiation cross section at each frequency
  call MFProb_RadiationSpectrum(chi1, nu0_HI,   ESpectrum, 1)
  call MFProb_RadiationSpectrum(chi2, nu0_HeI,  ESpectrum, 1)
  call MFProb_RadiationSpectrum(chi3, nu0_HeII, ESpectrum, 1)

  ! compute the scaling factor for the output emissivity
  !    set integration nodes
  do l=1,nnus
     nus(l) = 1 + (l-1)*(nu0_HeII-1.d0)/(nnus-1)
     etas(l) = 0.1d0 + (l-1)*(1.d0 - 1.d-8 - 0.1d0)/(nnus-1)
     nusB(l) = nu0_HeII/etas(l);
  enddo
  call MFProb_RadiationSpectrum(chis,  nus,  ESpectrum, nnus)
  call MFProb_RadiationSpectrum(chisB, nusB, ESpectrum, nnus)
  !    set the quadrature weights (composite Simpson)
  nbins = (nnus-1)/2
  wts = 0.d0
  wtsB = 0.d0
  do l=1,nbins
     i = 2*l-1
     j = 2*l
     k = 2*l+1
     wts(i) = wts(i) +      (nus(k)-nus(i))
     wts(j) = wts(j) + 4.d0*(nus(k)-nus(i))
     wts(k) = wts(k) +      (nus(k)-nus(i))
     wtsB(i) = wtsB(i) +      (etas(k)-etas(i))*nu0_HeII/etas(i)**2
     wtsB(j) = wtsB(j) + 4.d0*(etas(k)-etas(i))*nu0_HeII/etas(j)**2
     wtsB(k) = wtsB(k) +      (etas(k)-etas(i))*nu0_HeII/etas(k)**2
  enddo
  chinuint = (sum(wts*chis/nus/hp) + sum(wtsB*chisB/nusB/hp))/6.d0

  !    set integration nodes for free-streaming emissivity
  do l=1,nnus
     nus(l) = nu0_HI + (l-1)*(nu0_HeII-nu0_HI)/(nnus-1)
  enddo
  call MFProb_RadiationSpectrum(chis, nus, ESpectrum, nnus)
  !    set the quadrature weights (composite Simpson)
  wts = 0.d0
  do l=1,nbins
     i = 2*l-1
     j = 2*l
     k = 2*l+1
     wts(i) = wts(i) +      (nus(k)-nus(i))
     wts(j) = wts(j) + 4.d0*(nus(k)-nus(i))
     wts(k) = wts(k) +      (nus(k)-nus(i))
  enddo
  chiint = sum(wts*chis)/6.d0

  ! compute the effective NGammaDot for the free-streaming emissivity
  NGamDot_FS = NGammaDot*chiint/chinuint/h_nu1

  ! compute point source emissivity for various problems

  !   point-source emissivity at location (EtaCenter(1:3))
  if (ProbType == 460) then

!!$     print *,'MFProb Sources:'
!!$     print '(2(A,es9.2))','  NGammaDot =',NGammaDot,', chinuint =',chinuint*h_nu1
!!$     print '(2(A,es9.2))','  NGamDot_FS =',NGamDot_FS,', chiint =',chiint
!!$     print '(3(A,es9.2))','  chi1 =',chi1,', chi2 =',chi2,', chi3 =',chi3
!!$     print '(A,es9.2)','  EtaRadius =',EtaRadius
!!$     print '(A,3(es9.2))','  EtaCenter =',EtaCenter
!!$     print '(A,es9.2)','  dV =',dV

     ! one-cell source
     if (EtaRadius == 0.d0) then
        
        ! compute eta factor for given ionization source
        etaconst1 = NGammaDot*chi1/chinuint/dV
        etaconst2 = NGammaDot*chi2/chinuint/dV
        etaconst3 = NGammaDot*chi3/chinuint/dV
        
        ! place ionization source in one cell
        do k=1,Nz
           
           ! z-boundaries (comoving) for this cell
           cellZl = x2L + (k-1)*dz
           cellZr = cellZl + dz
           
           do j=1,Ny
              
              ! y-boundaries (comoving) for this cell
              cellYl = x1L + (j-1)*dy
              cellYr = cellYl + dy
              
              do i=1,Nx
                 
                 ! x-boundaries (comoving) for this cell
                 cellXl = x0L + (i-1)*dx
                 cellXr = cellXl + dx
                 
                 ! see if domain center is in cell (or on left edge)
                 if ( (cellXl <= EtaCenter(1)) .and. (cellXr > EtaCenter(1)) .and. &
                      (cellYl <= EtaCenter(2)) .and. (cellYr > EtaCenter(2)) .and. &
                      (cellZl <= EtaCenter(3)) .and. (cellZr > EtaCenter(3)) ) then
                    eta1(i,j,k) = etaconst1
                    eta2(i,j,k) = etaconst2
                    eta3(i,j,k) = etaconst3
                 endif
                 
              enddo
           enddo
        enddo

     else

        ! compute eta factor for given ionization source
        etaconst1 = NGammaDot*chi1/chinuint/dV/8.d0/(EtaRadius**3)
        etaconst2 = NGammaDot*chi2/chinuint/dV/8.d0/(EtaRadius**3)
        etaconst3 = NGammaDot*chi3/chinuint/dV/8.d0/(EtaRadius**3)
        
        ! place ionization source in center of domain
        do k=1,Nz
           
           ! z-center (comoving) for this cell
           cellZc = x2L + (k-0.5d0)*dz
           
           do j=1,Ny
              
              ! y-center (comoving) for this cell
              cellYc = x1L + (j-0.5d0)*dy
              
              do i=1,Nx
                 
                 ! x-center (comoving) for this cell
                 cellXc = x0L + (i-0.5d0)*dx
                 
                 ! see if cell is within source region
                 if ( (abs(cellXc-EtaCenter(1)) < EtaRadius*dx) .and. &
                      (abs(cellYc-EtaCenter(2)) < EtaRadius*dy) .and. &
                      (abs(cellZc-EtaCenter(3)) < EtaRadius*dz) ) then
                    eta1(i,j,k) = etaconst1
                    eta2(i,j,k) = etaconst2
                    eta3(i,j,k) = etaconst3
                 endif
                 
              enddo
           enddo
        enddo

     endif ! EtaRadius == 0

!!$     print *,'  max(eta1) =',maxval(eta1)
!!$     print *,'  max(eta2) =',maxval(eta2)
!!$     print *,'  max(eta3) =',maxval(eta3)

  !   emissivity flux along x=0 wall (NGammaDot photons/s/cm^2)
  else if (ProbType == 461) then

     ! place ionization source along left wall (if on this subdomain)
     if (x1L == 0.d0) then

        ! compute eta factor for given ionization source, and put on wall
        etaconst1 = NGammaDot*chi1/chinuint/dy
        etaconst2 = NGammaDot*chi2/chinuint/dy
        etaconst3 = NGammaDot*chi3/chinuint/dy
        eta1(1,:,:) = etaconst1
        eta2(1,:,:) = etaconst2
        eta3(1,:,:) = etaconst3

     endif
     
  !   point-source emissivity at center of every processor
  elseif (ProbType == 414) then

     ! compute eta factor for given ionization source
     etaconst1 = NGammaDot*chi1/chinuint/dV
     etaconst2 = NGammaDot*chi2/chinuint/dV
     etaconst3 = NGammaDot*chi3/chinuint/dV
        
     ! place ionization source in center of subdomain
     eta1(int(Nx/2),int(Ny/2),int(Nz/2)) = etaconst1
     eta2(int(Nx/2),int(Ny/2),int(Nz/2)) = etaconst2
     eta3(int(Nx/2),int(Ny/2),int(Nz/2)) = etaconst3
     
  endif ! ProbType

  return
end subroutine MFProb_Sources
!=======================================================================
