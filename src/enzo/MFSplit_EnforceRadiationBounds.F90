!=======================================================================
!
! Copyright 2009 Daniel R. Reynolds
!
! This software is released under the terms of the "Enzo Public License"
! in the accompanying LICENSE file.
!
!=======================================================================
subroutine MFSplit_EnforceRadiationBounds(Ef, E1, E2, E3, fsUn, E1Un, &
     E2Un, E3Un, Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr, NGzl, NGzr, ier)
!=======================================================================
!  PURPOSE: Enforces the constraints on the monochromatic radiation 
!           density values: 
!                0 < Ei < Eot(nu_i) = Ef*chi(nu_i)/chiint
!  NOTE: This routine is called *after* the E1, E2 and E3 arrays have 
!        been rescaled to internal solver units.
!  INPUTS:
!     Ef        - free-streaming radiation
!     E1        - monochromatic radiation density at nu_1
!     E2        - monochromatic radiation density at nu_2
!     E3        - monochromatic radiation density at nu_3
!     *Un       - unit conversion factors
!     N*        - grid information
!  OUTPUTS:
!     ier - success/failure flag (1->success, 0->failure)
!=======================================================================
  implicit none
#include "fortran.def"

  !--------------
  ! argument declarations
  integer, intent(in)  :: Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr, NGzl, NGzr
  real, intent(in) :: fsUn, E1Un, E2Un, E3Un
  real, intent(in) :: Ef(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr)
  real, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr) &
       :: E1, E2, E3
  integer, intent(out) :: ier
  
  !--------------
  ! locals
  real :: epsilon, half, zero
  integer :: i

  !=======================================================================

  ! initialize return flag
  ier = 1

  ! set radiation floor factor (since we never use Ei directly, only Ei/Ef)
  epsilon = 1.d0
  half = 0.5d0
  zero = 0.d0
  do i = 1,10000
     if (epsilon*half > zero) then
        epsilon = epsilon*half
     else
        exit
     end if
  enddo
  epsilon = epsilon**(0.33333d0)

  ! enforce bounds on monochromatic radiation fields
  E1 = max(epsilon*Ef,min(E1,Ef))
  E2 = max(epsilon*Ef,min(E2,Ef))
  E3 = max(epsilon*Ef,min(E3,Ef))

  return
end subroutine MFSplit_EnforceRadiationBounds
!=======================================================================
