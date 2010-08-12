!=======================================================================
!
! Copyright 2009 Daniel R. Reynolds
!
! This software is released under the terms of the "Enzo Public License"
! in the accompanying LICENSE file.
!
!=======================================================================
subroutine MFProb_RadInit(Ef, E1, E2, E3, EiScale, ESpectrum, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       October 2009
!
!  PURPOSE: Takes in the free-streaming radiation energy density Ef, 
!  and uses the assumed radiation spectrum defined by ESpectrum to
!  determine values of E1-E3.  In addition, it returns the scaling 
!  factor EiScale that should be used to non-dimensionalize the E1-E3 
!  units.
!
!  INPUTS:
!     Ef         - free-streaming radiation energy density
!     ESpectrum  - flag denoting choice of radiation spectrum
!
!  OUTPUTS: 
!     E1         - monochromatic radiation energy density at nu1
!     E2         - monochromatic radiation energy density at nu2
!     E3         - monochromatic radiation energy density at nu3
!     EiScale    - scaling factor to non-dimensionalize E1-E3 units
!     ier        - success/failure flag (1->success, 0->failure)
!
!=======================================================================
  implicit none
#include "fortran.def"

  !--------------
  ! argument declarations
  integer, intent(in) :: ESpectrum
  integer, intent(out) :: ier
  real, intent(in) :: Ef
  real, intent(out) :: E1, E2, E3, EiScale
  
  !--------------
  ! locals
  integer :: l, nbins
  integer, parameter :: nnus=999
  real :: nu1, nu2, nu3, chi1, chi2, chi3, chiint
  real, dimension(nnus) :: nus, chis, wts

  !=======================================================================

  ! initialize outputs to have all zero values, flag to success
  !    only overwrite emissivity sources if starmaker is not being used
  ier = 1
  E1 = 0.d0
  E2 = 0.d0
  E3 = 0.d0
  EiScale = 1.d0
  
  ! initialize constants
  nu1 = 13.6d0*1.60217653e-12/6.6260693d-27   ! ionization freq. of HI   [hz]
  nu2 = 24.6d0*1.60217653e-12/6.6260693d-27   ! ionization freq. of HeI  [hz]
  nu3 = 54.4d0*1.60217653e-12/6.6260693d-27   ! ionization freq. of HeII [hz]

  ! get the radiation spectrum at each frequency
  call MFProb_RadiationSpectrum(chi1, nu1, ESpectrum, 1)
  call MFProb_RadiationSpectrum(chi2, nu2, ESpectrum, 1)
  call MFProb_RadiationSpectrum(chi3, nu3, ESpectrum, 1)

  ! get the integrated radiation spectrum
  do l=1,nnus
     nus(l) = nu1 + 1.d1**((l-1)*18.d0/(nnus-1))
  enddo
  call MFProb_RadiationSpectrum(chis, nus, ESpectrum, nnus)
  !    set the quadrature weights (composite Simpson)
  wts = 0.d0
  nbins = (nnus+1)/2
  do l=1,nbins
     wts(2*l-1) = wts(2*l-1) + (nus(2*l+1)-nus(2*l-1))/6.d0
     wts(2*l)   = wts(2*l)   + (nus(2*l+1)-nus(2*l-1))*2.d0/3.d0
     wts(2*l+1) = wts(2*l+1) + (nus(2*l+1)-nus(2*l-1))/6.d0
  enddo
  chiint = sum(wts*chis)

  ! compute the monochromatic radiation energies
  E1 = Ef*chi1/chiint
  E2 = Ef*chi2/chiint
  E3 = Ef*chi3/chiint

  ! compute the monochromatic radiation scaling factor
  EiScale = chiint * 3.d0 / (chi1 + chi2 + chi3)
  
  return
end subroutine MFProb_RadInit
!=======================================================================
