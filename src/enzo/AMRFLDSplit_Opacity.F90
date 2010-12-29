!=======================================================================
!
! Copyright 2010 Daniel R. Reynolds
!
! This software is released under the terms of the "Enzo Public License"
! in the accompanying LICENSE file.
!
!=======================================================================
subroutine AMRFLDSplit_Opacity(kappaE, time, rho, n_HI, n_HeI, n_HeII, &
       a, IsE, IsEsHI, IsEsHInu, IsEsHeI, IsEsHeInu, IsEsHeII,         &
       IsEsHeIInu, aUnits, DenUnits, LenUnits, TimeUnits, NiUnits,     &
       Nchem, Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr, NGzl, NGzr, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       December 2010
!
!  PURPOSE: Computes the local opacities over the domain 
!           (NOTE: these are normalized by NiUnits, and must be 
!            converted to the appropriate redshift-dependent values 
!            in physics routines via kappaE_true = kappaE*NiUnits).
!
!  INPUTS:
!     time       - simulation time for evaluation
!     rho        - gas density
!     n_HI       - proportional density of HI species
!     n_HeI      - proportional density of HeI species
!     n_HeII     - proportional density of HeII species
!     a          - cosmological expansion parameter
!     IsE        - int_{nu0_HI}^{inf} sigE dnu
!     IsEsHI     - int_{nu0_HI}^{inf} sigE*sigHI dnu
!     IsEsHInu   - int_{nu0_HI}^{inf} sigE*sigHI/nu dnu
!     IsEsHeI    - int_{nu0_HeI}^{inf} sigE*sigHeI dnu
!     IsEsHeInu  - int_{nu0_HeI}^{inf} sigE*sigHeI/nu dnu
!     IsEsHeII   - int_{nu0_HeII}^{inf} sigE*sigHeII dnu
!     IsEsHeIInu - int_{nu0_HeII}^{inf} sigE*sigHeII/nu dnu
!     *Units     - variable scaling constants
!     Nchem      - number of chemistry species in problem {0, 1, 3}
!     Nx,Ny,Nz   - active mesh size in each direction
!     NG*l/NG*r  - left/right ghost cells in each direction
!
!     Note: the vector inputs are of size (Nx + NGxl + NGxr) in 
!     the x-direction, others are similar.
!
!  OUTPUT ARGUMENTS: 
!     kappaE     - Opacity coefficient (Energy mean)
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
  integer, intent(in) :: Nchem
  integer, intent(in) :: Nx, NGxl, NGxr
  integer, intent(in) :: Ny, NGyl, NGyr
  integer, intent(in) :: Nz, NGzl, NGzr
  integer, intent(out) :: ier
  REALSUB, intent(in) :: a
  real,    intent(in) :: time, IsE, IsEsHI, IsEsHInu, IsEsHeI
  real,    intent(in) :: IsEsHeInu, IsEsHeII, IsEsHeIInu
  real,    intent(in) :: aUnits, DenUnits, LenUnits, TimeUnits, NiUnits
  real, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr), intent(in) &
       :: rho, n_HI, n_HeI, n_HeII
  real, intent(out) :: kappaE(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr)

  !--------------
  ! locals
  integer :: i, j, k
  real    :: mp, eV2ergs, HIconst, HeIconst, HeIIconst
  real    :: rhoval, Tval

  !=======================================================================

!!$  write(*,*) 'Entering AMRFLDSplit::Opacity routine'

  ! initialize outputs to zero, flag to success
  kappaE = 0.d0
  ier = 1

  ! set shortcut values
  eV2ergs = 1.60217733d-12       ! coversion factor from eV to ergs
  mp = 1.67262171d-24            ! mass of a proton [g]

  ! compute opacity shortcuts
  HIconst   = IsEsHI/IsE
  HeIconst  = IsEsHeI/IsE/4.d0
  HeIIconst = IsEsHeII/IsE/4.d0

  ! compute opacity over domain depending on number of chemical species 

  ! Hydrogen only
  if (Nchem == 1) then
     do k=1-NGzl,Nz+NGzr,1
        do j=1-NGyl,Ny+NGyr,1
           do i=1-NGxl,Nx+NGxr,1
              kappaE(i,j,k) = n_HI(i,j,k)*HIconst
           enddo
        enddo
     enddo
     
  ! Hydrogen plus Helium
  else if (Nchem == 3) then
     do k=1-NGzl,Nz+NGzr,1
        do j=1-NGyl,Ny+NGyr,1
           do i=1-NGxl,Nx+NGxr,1
              kappaE(i,j,k) = n_HI(i,j,k)  *HIconst   &
                            + n_HeI(i,j,k) *HeIconst  &
                            + n_HeII(i,j,k)*HeIIconst
           enddo
        enddo
     enddo

  else
     write(0,*) 'AMRFLDSplit_Opacity ERROR: illegal Nchem =',Nchem, &
          ', requires Nchem = {1, 3}'
     ier = 0
  endif

  return

end subroutine AMRFLDSplit_Opacity
!=======================================================================
