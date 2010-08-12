!=======================================================================
!
! Copyright 2009 Daniel R. Reynolds
!
! This software is released under the terms of the "Enzo Public License"
! in the accompanying LICENSE file.
!
!=======================================================================
subroutine MFSplit_RadiationSpectrum(chi, nu, ESpectrum, nnus)
!=======================================================================
!  PURPOSE: Returns the radiation spectrum chi(nu), that depends on the 
!           user input ESpectrum:
!                1 -> T=1e5 blackbody spectrum
!                default -> power law spectrum with power -1.5
!  INPUTS:
!     nu         - frequencies
!     ESpectrum  - spectrum choice
!     nnus       - length of nu
!  OUTPUT ARGUMENTS: 
!     chi        - spectrum values
!=======================================================================
  implicit none
#include "fortran.def"

  !--------------
  ! argument declarations
  integer, intent(in) :: ESpectrum, nnus
  real, intent(in)    :: nu(nnus)
  real, intent(out)   :: chi(nnus)
  
  !=======================================================================

  ! evaluate the radiation spectrum based on the ESpectrum parameter 
  ! (uses f90 whole array operations)
  if (Espectrum == 1) then
     ! T = 1e5 K blackbody spectrum
     !   8*pi*hp*(nu/c)**3/(exp(nu*hp/kb/1e5)-1)
     chi = 1.66531285080455d-25 * (nu/2.99792458d10)**3 &
          / (exp(nu*4.79923759121064d-16)-1.d0)
  else
     ! simple power law spectrum with power -1.5
     !   (nu/nu0)**(-1.5)
     chi = (nu*3.04093193738880d-16)**(-1.5d0)
  endif
  
  return
end subroutine MFSplit_RadiationSpectrum
!=======================================================================




subroutine MFSplit_HICrossSection(sig, nu, nnus)
!=======================================================================
!  PURPOSE: Returns the HI cross section sig(nu). 
!  INPUTS:    nu - frequencies, nnus - length of nu
!  OUTPUT:   sig - spectrum values
!=======================================================================
  implicit none
#include "fortran.def"

  !--------------
  ! argument declarations
  integer, intent(in) :: nnus
  real, intent(in)    :: nu(nnus)
  real, intent(out)   :: sig(nnus)
  
  !--------------
  ! locals
  integer :: i
  real*8, parameter :: pi=3.141592653589793238d0
  real*8, parameter :: nu0_HI=3.28846557762383d15    ! HI ionization (hz)
  real*8  :: eps

  !=======================================================================

  do i=1,nnus
     if (nu(i) <= nu0_HI) then
        sig(i) = 6.3d-18
     else
        eps = sqrt(nu(i)/nu0_HI - 1.d0)
        sig(i) = 6.3d-18 * (nu0_HI/nu(i))**4    &
                * exp(4.d0-4.d0*atan(eps)/eps)  &
                / (1.d0-exp(-2.d0*pi/eps))
     endif
  enddo
  
  return
end subroutine MFSplit_HICrossSection
!=======================================================================




subroutine MFSplit_HeICrossSection(sig, nu, nnus)
!=======================================================================
!  PURPOSE: Returns the HeI cross section sig(nu). 
!  INPUTS:    nu - frequencies, nnus - length of nu
!  OUTPUT:   sig - spectrum values
!=======================================================================
  implicit none
#include "fortran.def"

  !--------------
  ! argument declarations
  integer, intent(in) :: nnus
  real, intent(in)    :: nu(nnus)
  real, intent(out)   :: sig(nnus)
  
  !--------------
  ! locals
  integer :: i
  real*8, parameter :: pi=3.141592653589793238d0
  real*8, parameter :: nu0_HeI=5.94825391246663d15   ! HeI ionization (hz)
  real*8  :: eps

  !=======================================================================

  do i=1,nnus
     if (nu(i) <= nu0_HeI) then
        sig(i) = 7.42d-18
     else
        sig(i) = 7.42e-18 * (1.66d0*(nu0_HeI/nu(i))**(2.05d0)  &
                           - 0.66d0*(nu0_HeI/nu(i))**(3.05d0))
     endif
  enddo
  
  return
end subroutine MFSplit_HeICrossSection
!=======================================================================




subroutine MFSplit_HeIICrossSection(sig, nu, nnus)
!=======================================================================
!  PURPOSE: Returns the HeII cross section sig(nu). 
!  INPUTS:    nu - frequencies, nnus - length of nu
!  OUTPUT:   sig - spectrum values
!=======================================================================
  implicit none
#include "fortran.def"

  !--------------
  ! argument declarations
  integer, intent(in) :: nnus
  real, intent(in)    :: nu(nnus)
  real, intent(out)   :: sig(nnus)
  
  !--------------
  ! locals
  integer :: i
  real*8, parameter :: pi=3.141592653589793238d0
  real*8, parameter :: nu0_HeII=1.31538623104953d16  ! HeII ionization (hz)
  real*8  :: eps

  !=======================================================================

  do i=1,nnus
     if (nu(i) <= nu0_HeII) then
        sig(i) = 1.575d-18
     else
        eps = sqrt(nu(i)/nu0_HeII - 1.d0)
        sig(i) = 1.575d-18 * (nu0_HeII/nu(i))**4  &
                * exp(4.d0-4.d0*atan(eps)/eps) &
                / (1.d0-exp(-2.d0*pi/eps))
     endif
  enddo
  
  return
end subroutine MFSplit_HeIICrossSection
!=======================================================================
