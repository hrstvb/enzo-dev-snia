!=======================================================================
!
! Copyright 2009 Daniel R. Reynolds
!
! This software is released under the terms of the "Enzo Public License"
! in the accompanying LICENSE file.
!
!=======================================================================
subroutine MFSplit_FSRadiationMask(Ef, E1, E2, E3, fsUn, E1Un, E2Un, &
     E3Un, Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr, NGzl, NGzr, ier)
!=======================================================================
!  PURPOSE: cuts off the free-streaming radiation outside of the 
!           "ionized" region.
!  NOTE: All values are given in normalized solver units.
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
  real, intent(in), dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr) &
       :: E1, E2, E3
  real :: Ef(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr)
  integer, intent(out) :: ier
  
  !--------------
  ! locals
  integer :: i, j, k
  real :: epsilon, half, zero, Ef_floor, E1min, E2min, E3min

  !=======================================================================

  ! initialize return flag to success
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
 
  ! set radiation floor/cutoff values (normalized units)
  Ef_floor = 1.0d-50
  E1min = 1.0d-16
  E2min = 1.0d-16
  E3min = 1.0d-16
!!$  E1min = epsilon
!!$  E2min = epsilon
!!$  E3min = epsilon

  ! loop over domain
  do k = 1,Nz,1
     do j = 1,Ny,1
        do i = 1,Nx,1

           ! How do I define the "ionized" region where Ef is retained???  
           ! If I use HII/H + (HeII+HeIII)/He then we have a chicken/egg 
           !   problem, since radiation must first ionize a cell, but a cell 
           !   must be ionized to get radiation.  
           ! Instead, we should use some measure of the monochromatic 
           !   radiation values, but what?  All approaches are some kind of 
           !   kluge, with little/no theoretical justification whatsoever.  
           ! Moreover, will a cutoff actually save us from the fact that we 
           !   don't have free-streaming BCs to suck radiation out of the 
           !   domain (periodic and Neumann will just fill in the whole 
           !   domain while Dirichlet requires some additional approximation 
           !   to get the *value* to force at the boundary).

!!$           ! Option 1: all monochromatic values must be above cutoff
!!$           if ((E1(i,j,k)<E1min) .or. (E2(i,j,k)<E2min) .or. (E3(i,j,k)<E3min)) then
!!$              Ef(i,j,k) = Ef_floor
!!$           end if

           ! Option 2: any monochromatic value must be above cutoff
           if ((E1(i,j,k)<E1min) .and. (E2(i,j,k)<E2min) .and. (E3(i,j,k)<E3min)) then
              Ef(i,j,k) = Ef_floor
           end if
           
           ! Option 3: do nothing, and just live with the Ef values we get

        end do
     end do
  end do

  return
end subroutine MFSplit_FSRadiationMask
!=======================================================================
