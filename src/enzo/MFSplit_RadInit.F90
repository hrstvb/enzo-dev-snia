!=======================================================================
!
! Copyright 2009 Daniel R. Reynolds
!
! This software is released under the terms of the "Enzo Public License"
! in the accompanying LICENSE file.
!
!=======================================================================
subroutine MFSplit_RadInit(Ef, E1, E2, E3, E1Units, E2Units, E3Units, &
                           chiint, chinuint, ESpectrum, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       October 2009
!
!  PURPOSE: Takes in the free-streaming radiation energy density Ef, 
!  and uses the assumed radiation spectrum defined by ESpectrum to
!  determine values of E1-E3.  In addition, it returns the scaling 
!  factors E*Units that should be used to non-dimensionalize the E1-E3 
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
!     E1Units    - scaling factor to non-dimensionalize E1 units
!     E2Units    - scaling factor to non-dimensionalize E2 units
!     E3Units    - scaling factor to non-dimensionalize E3 units
!     chiint     - int_nu1^infty chi(nu) dnu
!     chinuint   - int_0^infty chi(nu) nu1/nu dnu
!     ier        - success/failure flag (1->success, 0->failure)
!
!=======================================================================
  implicit none
#include "fortran.def"

  !--------------
  ! argument declarations
  integer, intent(in) :: ESpectrum
  integer, intent(out) :: ier
  REAL, intent(in) :: Ef
  REAL, intent(out) :: E1, E2, E3, E1Units, E2Units, E3Units, chiint, chinuint
  
  !--------------
  ! locals
  integer :: i, j, k, l, nbins
  integer, parameter :: nnus=999
  REAL :: nu1, nu2, nu3, chi1, chi2, chi3
  REAL, dimension(nnus) :: nus, etas, nusB, chis, chisB, wts, wtsB

  !=======================================================================

  ! flag to success
  ier = 1
  
  ! initialize constants
  nu1 = 13.6d0*1.60217653e-12/6.6260693d-27   ! ionization freq. of HI   [hz]
  nu2 = 24.6d0*1.60217653e-12/6.6260693d-27   ! ionization freq. of HeI  [hz]
  nu3 = 54.4d0*1.60217653e-12/6.6260693d-27   ! ionization freq. of HeII [hz]

  ! get the radiation spectrum at each frequency
  call MFSplit_RadiationSpectrum(chi1, nu1, ESpectrum, 1)
  call MFSplit_RadiationSpectrum(chi2, nu2, ESpectrum, 1)
  call MFSplit_RadiationSpectrum(chi3, nu3, ESpectrum, 1)


  ! compute the emissivity scaling factor 
  ! (computed here once and used repeatedly in emissivity sources)
  !    set integration nodes
  do l=1,nnus
     nus(l)  = 1.d0 + (nu3-1.d0)*(l-1)/(nnus-1)
     etas(l) = 0.1d0 + (1.d0 - 1.d-8 - 0.1d0)*(l-1)/(nnus-1)
     nusB(l) = nu3/etas(l);
  enddo
  call MFSplit_RadiationSpectrum(chis,  nus,  ESpectrum, nnus)
  call MFSplit_RadiationSpectrum(chisB, nusB, ESpectrum, nnus)
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
     wtsB(i) = wtsB(i) +      (etas(k)-etas(i))*nu3/etas(i)**2
     wtsB(j) = wtsB(j) + 4.d0*(etas(k)-etas(i))*nu3/etas(j)**2
     wtsB(k) = wtsB(k) +      (etas(k)-etas(i))*nu3/etas(k)**2
  enddo
  chinuint = (sum(wts*chis/nus) + sum(wtsB*chisB/nusB))*nu1/6.d0


  ! get the integrated radiation spectrum
  !    set integration nodes
  do l=1,nnus
     nus(l) = nu1 + 1.d1**((l-1)*18.d0/(nnus-1))
  enddo
  call MFSplit_RadiationSpectrum(chis, nus, ESpectrum, nnus)
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

  ! compute the monochromatic radiation scaling factors
  E1Units = chi1/chiint
  E2Units = chi2/chiint
  E3Units = chi3/chiint
  
  ! compute the monochromatic radiation energies
  E1 = Ef*E1Units
  E2 = Ef*E2Units
  E3 = Ef*E3Units

  return
end subroutine MFSplit_RadInit
!=======================================================================
