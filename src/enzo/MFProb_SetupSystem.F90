!=======================================================================
!
! Copyright 2009 Daniel R. Reynolds
!
! This software is released under the terms of the "Enzo Public License"
! in the accompanying LICENSE file.
!
!=======================================================================
subroutine MFProb_SetupSystem(mat11, mat12, mat13, mat21, mat22, mat23, &
     mat31, mat32, mat33, rhs1, rhs2, rhs3, E1, E2, E3, E10, E20, E30,  &
     HI, HI0, HeI, HeI0, HeII, HeII0, adj11, adj12, adj13, adj21,       &
     adj22, adj23, adj31, adj32, adj33, LType, LImp, dt, theta, a, lUn, &
     rUn, nUn, nUn0, dx, dy, dz, BCValsXl, BCValsXr, BCValsYl,          &
     BCValsYr, BCValsZl, BCValsZr, BCXl, BCXr, BCYl, BCYr, BCZl, BCZr,  &
     x0s, x0e, x1s, x1e, x2s, x2e, Nchem, Nx, Ny, Nz, NGxl, NGxr, NGyl, &
     NGyr, NGzl, NGzr, xlface, xrface, ylface, yrface, zlface, zrface, ier)
  !=======================================================================
  !  written by: Daniel R. Reynolds
  !  date:       August 2009
  !  modified1:  
  !
  !  PURPOSE: Computes the array of matrix stencil elements for the 
  !           multi-frequency problem,
  !              E1dot = Div(D(E1)*Grad(E1)) + 4*pi*src1 - c*k1*E1,
  !              E2dot = Div(D(E2)*Grad(E2)) + 4*pi*src2 - c*k2*E2,
  !              E3dot = Div(D(E3)*Grad(E3)) + 4*pi*src3 - c*k3*E3,
  !           where D(E) is a nonlinear flux-limiter depending on E.  
  !           We define the values
  !              R_i = |Grad(E)_i|/k/E,
  !              k = absorption coefficient
  !           The '_i' subscript implies the gradient in the ith 
  !           direction; these quantities are all required at cell faces, 
  !           as that is the location of the divergence calculations.
  !           With these components, we allow any of the following three 
  !           forms of the limiter, 
  !             [Levermore-Pomraning, 1981],
  !                 D_i(E) = c/k/R_i*[coth(R_i)-1/R_i],
  !             [rational approx. to above, Levermore-Pomraning, 1981],
  !                 D_i(E) = c/k*(2+R_i)/(6+3*R_i+R_i**2),
  !             [Reynolds approximation to LP],
  !                 D_i(E) = 2/pi*c/k/R_i*atan(R_i*pi/6),
  !             [Zeus form of rational approx. to LP],
  !                 D_i(E) = 2/pi*c/k/R_i*atan(R_i*pi/6),
  !
  !           As the stencil has 7 non-zero elements per matrix row, we 
  !           set these entries over the computational domain, with the 
  !           proper adjustments due to the choice of limiter.
  !
  !  INPUTS:
  !     rhs*   - arrays of rhs values, same size as variables
  !     E1,E10     - Radiation energy density, frequency 1 (new & old)
  !     E2,E20     - Radiation energy density, frequency 2 (new & old)
  !     E3,E30     - Radiation energy density, frequency 3 (new & old)
  !     HI,HI0     - Hydrogen I density (new & old)
  !     HeI,HeI0   - Helium I density (new & old)
  !     HeII,HeII0 - Helium II density (new & old)
  !     adj1*      - Schur complement adjustment vectors (for E1)
  !     adj2*      - Schur complement adjustment vectors (for E2)
  !     adj3*      - Schur complement adjustment vectors (for E3)
  !     LType      - integer flag denoting type of flux limiter:
  !                       0 -> standard Levermore-Pomraning lim. (LP, 1981)
  !                       1 -> rational approx. to LP lim. (LP, 1981)
  !                       2 -> Reynolds approx to LP lim.
  !                       3 -> turns off the limiter (constant of 1/3)
  !                       4 -> Zeus limiter
  !     LImp       - integer flag denoting implicitness of flux limiter:
  !                       0 -> fully lagged to previous time step
  !                       1 -> fully lagged to previous newton iterate
  !                       2 -> lag only opacity dependence (unimplemented)
  !     a          - cosmological expansion parameter
  !     dt         - time step size
  !     theta      - overall implicitness parameter
  !     *Un        - variable scaling constants
  !     dx,dy,dz   - mesh spacing in each direction
  !     BCVals*    - boundary condition value arrays for each face
  !     BC*        - boundary condition type in each direction, face
  !                     0->periodic
  !                     1->Dirichlet
  !                     2->Neumann
  !     x*{s,e}    - start/end indices of linear solver domain; 
  !                  typically 1:Nx for standard dims, but Dirichlet 
  !                  BCs may move these to 0:Nx, 1:Nx+1, etc.
  !     Nx,Ny,Nz   - active mesh size in each direction
  !     NG*l/NG*r  - left/right ghost cells in each direction
  !     *{l,r}face - integer flag denoting whether direction/face 
  !                  is external to the domain (0->int, 1->ext)
  !
  !     Note: the vector inputs are of size (Nx + NGxl + NGxr) in 
  !     the x-direction, others are similar.
  !
  !  OUTPUT ARGUMENTS: 
  !     mat11  - array of stencil values over the active domain 
  !                  for the frequency 1 wrt 1. Since each spatial stencil has
  !                  7 nonzero entries, and as this array should not include 
  !                  ghost cells, it has dimensions (7,Nx,Ny,Nz).
  !     mat12  - array of stencil values over the active domain 
  !                  for the frequency 1 wrt 2, has dimensions (Nx,Ny,Nz).
  !     mat13  - array of stencil values over the active domain 
  !                  for the frequency 1 wrt 3, has dimensions (Nx,Ny,Nz).
  !     mat2*  - array of stencil values for frequency 2
  !     mat3*  - array of stencil values for frequency 3
  !     rhs*   - arrays of rhs values adjusted for BCs, same size as variables
  !     ier    - success/failure flag (0->failure, 1->success)
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
  integer, intent(in)  :: LType, LImp, Nchem
  integer,  intent(in) :: BCXl, BCXr, x0s, x0e, Nx, NGxl, NGxr, xlface, xrface
  integer,  intent(in) :: BCYl, BCYr, x1s, x1e, Ny, NGyl, NGyr, ylface, yrface
  integer,  intent(in) :: BCZl, BCZr, x2s, x2e, Nz, NGzl, NGzr, zlface, zrface
  integer, intent(out) :: ier
  REALSUB, intent(in)  :: a
  real, intent(in) :: dx, dy, dz, dt, theta
  real, intent(in) :: lUn, rUn, nUn, nUn0
  real, intent(in),  dimension(*) :: E1, E2, E3, E10, E20, E30
  real, intent(in),  dimension(*) :: HI, HI0, HeI, HeI0, HeII, HeII0
  real, intent(in) :: BCValsXl(*), BCValsYl(*), BCValsZl(*)
  real, intent(in) :: BCValsXr(*), BCValsYr(*), BCValsZr(*)
  real*8,   intent(out), dimension(*) :: mat11, mat12, mat13, mat21, mat22, &
       mat23, mat31, mat32, mat33
  real*8,   dimension(*) :: rhs1, rhs2, rhs3
  real, intent(in),  dimension(*) :: adj11, adj12, adj13, adj21, adj22, &
       adj23, adj31, adj32, adj33
 
  !--------------
  ! local declarations
  real :: hp, ev2erg, nu0_HI, nu0_HeI, nu0_HeII, sHI, sHeI, sHeII

  !=======================================================================

  ! set shortcut values 
  hp = 6.6260693d-27            ! Planck's constant (ergs*s)
  ev2erg = 1.60217653d-12       ! conversion constant from eV to ergs
  nu0_HI   = 13.6d0*ev2erg/hp   ! ionization threshold of HI (hz)
  nu0_HeI  = 24.6d0*ev2erg/hp   ! ionization threshold of HeI (hz)
  nu0_HeII = 54.4d0*ev2erg/hp   ! ionization threshold of HeII (hz)

  ! pass arguments based on Nchem
  if (Nchem == 1) then
  
     !!!!!!! E1 matrix !!!!!!!
     !   set shortcut variables
     call MFProb_HICrossSection(sHI, nu0_HI, 1)
     sHeI  = 0.d0
     sHeII = 0.d0

     ! call the appropriate dimension-specific routine
!!$     print *, 'Setting up E1 matrix, sHI =',sHI
     if ((Nz == 1) .and. (Ny == 1)) then
        call MFProb_SetupSystem_1D(mat11, mat12, mat13, rhs1, E1, E10, HI, HI0,    &
             HI, HI0, HI, HI0, adj11, adj12, adj13, LType, LImp, dt, theta, sHI,   &
             sHeI, sHeII, a, lUn, rUn, nUn, nUn0, dx, BCValsXl, BCValsXr, BCXl,    &
             BCXr, x0s, x0e, Nx, NGxl, NGxr, xlface, xrface, ier)
     elseif (Nz == 1) then
        call MFProb_SetupSystem_2D(mat11, mat12, mat13, rhs1, E1, E10, HI, HI0,    &
             HI, HI0, HI, HI0, adj11, adj12, adj13, LType, LImp, dt, theta, sHI,   &
             sHeI, sHeII, a, lUn, rUn, nUn, nUn0, dx, dy, BCValsXl, BCValsXr,      &
             BCValsYl, BCValsYr, BCXl, BCXr, BCYl, BCYr, x0s, x0e, x1s, x1e, Nx,   &
             Ny, NGxl, NGxr, NGyl, NGyr, xlface, xrface, ylface, yrface, ier)
     else
        call MFProb_SetupSystem_3D(mat11, mat12, mat13, rhs1, E1, E10, HI, HI0,    &
             HI, HI0, HI, HI0, adj11, adj12, adj13, LType, LImp, dt, theta, sHI,   &
             sHeI, sHeII, a, lUn, rUn, nUn, nUn0, dx, dy, dz, BCValsXl, BCValsXr,  &
             BCValsYl, BCValsYr, BCValsZl, BCValsZr, BCXl, BCXr, BCYl, BCYr, BCZl, &
             BCZr, x0s, x0e, x1s, x1e, x2s, x2e, Nx, Ny, Nz, NGxl, NGxr, NGyl,     &
             NGyr, NGzl, NGzr, xlface, xrface, ylface, yrface, zlface, zrface, ier)
     end if
     
     !!!!!!! E2 evaluation !!!!!!!
     !   set pointers and shortcut variables
     call MFProb_HICrossSection(sHI, nu0_HeI, 1)
     
     ! call the appropriate dimension-specific routine
!!$     print *, 'Setting up E2 matrix, sHI =',sHI
     if ((Nz == 1) .and. (Ny == 1)) then
        call MFProb_SetupSystem_1D(mat22, mat21, mat23, rhs2, E2, E20, HI, HI0,    &
             HI, HI0, HI, HI0, adj22, adj21, adj23, LType, LImp, dt, theta, sHI,   &
             sHeI, sHeII, a, lUn, rUn, nUn, nUn0, dx, BCValsXl, BCValsXr, BCXl,    &
             BCXr, x0s, x0e, Nx, NGxl, NGxr, xlface, xrface, ier)
     elseif (Nz == 1) then
        call MFProb_SetupSystem_2D(mat22, mat21, mat23, rhs2, E2, E20, HI, HI0,    &
             HI, HI0, HI, HI0, adj22, adj21, adj23, LType, LImp, dt, theta, sHI,   &
             sHeI, sHeII, a, lUn, rUn, nUn, nUn0, dx, dy, BCValsXl, BCValsXr,      &
             BCValsYl, BCValsYr, BCXl, BCXr, BCYl, BCYr, x0s, x0e, x1s, x1e, Nx,   &
             Ny, NGxl, NGxr, NGyl, NGyr, xlface, xrface, ylface, yrface, ier)
     else
        call MFProb_SetupSystem_3D(mat22, mat21, mat23, rhs2, E2, E20, HI, HI0,    &
             HI, HI0, HI, HI0, adj22, adj21, adj23, LType, LImp, dt, theta, sHI,   &
             sHeI, sHeII, a, lUn, rUn, nUn, nUn0, dx, dy, dz, BCValsXl, BCValsXr,  &
             BCValsYl, BCValsYr, BCValsZl, BCValsZr, BCXl, BCXr, BCYl, BCYr, BCZl, &
             BCZr, x0s, x0e, x1s, x1e, x2s, x2e, Nx, Ny, Nz, NGxl, NGxr, NGyl,     &
             NGyr, NGzl, NGzr, xlface, xrface, ylface, yrface, zlface, zrface, ier)
     end if
     

     !!!!!!! E3 evaluation !!!!!!!
     !   set pointers and shortcut variables
     call MFProb_HICrossSection(sHI, nu0_HeII, 1)
     
     ! call the appropriate dimension-specific routine
!!$     print *, 'Setting up E3 matrix, sHI =',sHI
     if ((Nz == 1) .and. (Ny == 1)) then
        call MFProb_SetupSystem_1D(mat33, mat31, mat32, rhs3, E3, E30, HI, HI0,    &
             HI, HI0, HI, HI0, adj33, adj31, adj32, LType, LImp, dt, theta, sHI,   &
             sHeI, sHeII, a, lUn, rUn, nUn, nUn0, dx, BCValsXl, BCValsXr, BCXl,    &
             BCXr, x0s, x0e, Nx, NGxl, NGxr, xlface, xrface, ier)
     elseif (Nz == 1) then
        call MFProb_SetupSystem_2D(mat33, mat31, mat32, rhs3, E3, E30, HI, HI0,    &
             HI, HI0, HI, HI0, adj33, adj31, adj32, LType, LImp, dt, theta, sHI,   &
             sHeI, sHeII, a, lUn, rUn, nUn, nUn0, dx, dy, BCValsXl, BCValsXr,      &
             BCValsYl, BCValsYr, BCXl, BCXr, BCYl, BCYr, x0s, x0e, x1s, x1e, Nx,   &
             Ny, NGxl, NGxr, NGyl, NGyr, xlface, xrface, ylface, yrface, ier)
     else
        call MFProb_SetupSystem_3D(mat33, mat31, mat32, rhs3, E3, E30, HI, HI0,    &
             HI, HI0, HI, HI0, adj33, adj31, adj32, LType, LImp, dt, theta, sHI,   &
             sHeI, sHeII, a, lUn, rUn, nUn, nUn0, dx, dy, dz, BCValsXl, BCValsXr,  &
             BCValsYl, BCValsYr, BCValsZl, BCValsZr, BCXl, BCXr, BCYl, BCYr, BCZl, &
             BCZr, x0s, x0e, x1s, x1e, x2s, x2e, Nx, Ny, Nz, NGxl, NGxr, NGyl,     &
             NGyr, NGzl, NGzr, xlface, xrface, ylface, yrface, zlface, zrface, ier)
     end if

  else  ! Nchem == 3

     !!!!!!! E1 matrix !!!!!!!
     !   set shortcut variables
     call MFProb_HICrossSection(sHI, nu0_HI, 1)
     sHeI  = 0.d0
     sHeII = 0.d0

     ! call the appropriate dimension-specific routine
!!$     print *, 'Setting up E1 matrix, sHI =',sHI
     if ((Nz == 1) .and. (Ny == 1)) then
        call MFProb_SetupSystem_1D(mat11, mat12, mat13, rhs1, E1, E10, HI, HI0, HeI, &
             HeI0, HeII, HeII0, adj11, adj12, adj13, LType, LImp, dt, theta, sHI,    &
             sHeI, sHeII, a, lUn, rUn, nUn, nUn0, dx, BCValsXl, BCValsXr, BCXl,      &
             BCXr, x0s, x0e, Nx, NGxl, NGxr, xlface, xrface, ier)
     elseif (Nz == 1) then
        call MFProb_SetupSystem_2D(mat11, mat12, mat13, rhs1, E1, E10, HI, HI0, HeI, &
             HeI0, HeII, HeII0, adj11, adj12, adj13, LType, LImp, dt, theta, sHI,    &
             sHeI, sHeII, a, lUn, rUn, nUn, nUn0, dx, dy, BCValsXl, BCValsXr,        &
             BCValsYl, BCValsYr, BCXl, BCXr, BCYl, BCYr, x0s, x0e, x1s, x1e, Nx, Ny, &
             NGxl, NGxr, NGyl, NGyr, xlface, xrface, ylface, yrface, ier)
     else
        call MFProb_SetupSystem_3D(mat11, mat12, mat13, rhs1, E1, E10, HI, HI0, HeI, &
             HeI0, HeII, HeII0, adj11, adj12, adj13, LType, LImp, dt, theta, sHI,    &
             sHeI, sHeII, a, lUn, rUn, nUn, nUn0, dx, dy, dz, BCValsXl, BCValsXr,    &
             BCValsYl, BCValsYr, BCValsZl, BCValsZr, BCXl, BCXr, BCYl, BCYr, BCZl,   &
             BCZr, x0s, x0e, x1s, x1e, x2s, x2e, Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr, &
             NGzl, NGzr, xlface, xrface, ylface, yrface, zlface, zrface, ier)
     end if
     
     !!!!!!! E2 evaluation !!!!!!!
     !   set pointers and shortcut variables
     call MFProb_HICrossSection(sHI, nu0_HeI, 1)
     call MFProb_HeICrossSection(sHeI, nu0_HeI, 1)
     
     ! call the appropriate dimension-specific routine
!!$     print *, 'Setting up E2 matrix, sHI =',sHI
!!$     print *, 'Setting up E2 matrix, sHeI =',sHeI
     if ((Nz == 1) .and. (Ny == 1)) then
        call MFProb_SetupSystem_1D(mat22, mat21, mat23, rhs2, E2, E20, HI, HI0, HeI, &
             HeI0, HeII, HeII0, adj22, adj21, adj23, LType, LImp, dt, theta, sHI,    &
             sHeI, sHeII, a, lUn, rUn, nUn, nUn0, dx, BCValsXl, BCValsXr, BCXl,      &
             BCXr, x0s, x0e, Nx, NGxl, NGxr, xlface, xrface, ier)
     elseif (Nz == 1) then
        call MFProb_SetupSystem_2D(mat22, mat21, mat23, rhs2, E2, E20, HI, HI0, HeI, &
             HeI0, HeII, HeII0, adj22, adj21, adj23, LType, LImp, dt, theta, sHI,    &
             sHeI, sHeII, a, lUn, rUn, nUn, nUn0, dx, dy, BCValsXl, BCValsXr,        &
             BCValsYl, BCValsYr, BCXl, BCXr, BCYl, BCYr, x0s, x0e, x1s, x1e, Nx, Ny, &
             NGxl, NGxr, NGyl, NGyr, xlface, xrface, ylface, yrface, ier)
     else
        call MFProb_SetupSystem_3D(mat22, mat21, mat23, rhs2, E2, E20, HI, HI0, HeI, &
             HeI0, HeII, HeII0, adj22, adj21, adj23, LType, LImp, dt, theta, sHI,    &
             sHeI, sHeII, a, lUn, rUn, nUn, nUn0, dx, dy, dz, BCValsXl, BCValsXr,    &
             BCValsYl, BCValsYr, BCValsZl, BCValsZr, BCXl, BCXr, BCYl, BCYr, BCZl,   &
             BCZr, x0s, x0e, x1s, x1e, x2s, x2e, Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr, &
             NGzl, NGzr, xlface, xrface, ylface, yrface, zlface, zrface, ier)
     end if
     

     !!!!!!! E3 evaluation !!!!!!!
     !   set pointers and shortcut variables
     call MFProb_HICrossSection(sHI, nu0_HeII, 1)
     call MFProb_HeICrossSection(sHeI, nu0_HeII, 1)
     call MFProb_HeIICrossSection(sHeII, nu0_HeII, 1)
     
     ! call the appropriate dimension-specific routine
!!$     print *, 'Setting up E2 matrix, sHI =',sHI
!!$     print *, 'Setting up E2 matrix, sHeI =',sHeI
!!$     print *, 'Setting up E2 matrix, sHeII =',sHeII
     if ((Nz == 1) .and. (Ny == 1)) then
        call MFProb_SetupSystem_1D(mat33, mat31, mat32, rhs3, E3, E30, HI, HI0, HeI, &
             HeI0, HeII, HeII0, adj33, adj31, adj32, LType, LImp, dt, theta, sHI,    &
             sHeI, sHeII, a, lUn, rUn, nUn, nUn0, dx, BCValsXl, BCValsXr, BCXl,      &
             BCXr, x0s, x0e, Nx, NGxl, NGxr, xlface, xrface, ier)
     elseif (Nz == 1) then
        call MFProb_SetupSystem_2D(mat33, mat31, mat32, rhs3, E3, E30, HI, HI0, HeI, &
             HeI0, HeII, HeII0, adj33, adj31, adj32, LType, LImp, dt, theta, sHI,    &
             sHeI, sHeII, a, lUn, rUn, nUn, nUn0, dx, dy, BCValsXl, BCValsXr,        &
             BCValsYl, BCValsYr, BCXl, BCXr, BCYl, BCYr, x0s, x0e, x1s, x1e, Nx, Ny, &
             NGxl, NGxr, NGyl, NGyr, xlface, xrface, ylface, yrface, ier)
     else
        call MFProb_SetupSystem_3D(mat33, mat31, mat32, rhs3, E3, E30, HI, HI0, HeI, &
             HeI0, HeII, HeII0, adj33, adj31, adj32, LType, LImp, dt, theta, sHI,    &
             sHeI, sHeII, a, lUn, rUn, nUn, nUn0, dx, dy, dz, BCValsXl, BCValsXr,    &
             BCValsYl, BCValsYr, BCValsZl, BCValsZr, BCXl, BCXr, BCYl, BCYr, BCZl,   &
             BCZr, x0s, x0e, x1s, x1e, x2s, x2e, Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr, &
             NGzl, NGzr, xlface, xrface, ylface, yrface, zlface, zrface, ier)
     end if

  endif  ! Nchem

  
  return
end subroutine MFProb_SetupSystem
!=======================================================================






subroutine MFProb_SetupSystem_3D(mat1, mat2, mat3, rhs, E, E0, HI, HI0, HeI, &
     HeI0, HeII, HeII0, adj1, adj2, adj3, LType, LImp, dt, theta, sHI,    &
     sHeI, sHeII, a, lUn, rUn, nUn, nUn0, dx, dy, dz, BCValsXl, BCValsXr,    &
     BCValsYl, BCValsYr, BCValsZl, BCValsZr, BCXl, BCXr, BCYl, BCYr, BCZl,   &
     BCZr, x0s, x0e, x1s, x1e, x2s, x2e, Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr, &
     NGzl, NGzr, xlface, xrface, ylface, yrface, zlface, zrface, ier)
  !=======================================================================
  !  PURPOSE: 3D version of the routine
  !=======================================================================
#include "fortran.def"
  implicit none
  
  !--------------
  ! argument declarations
  integer, intent(in)  :: LType, LImp
  integer,  intent(in) :: BCXl, BCXr, x0s, x0e, Nx, NGxl, NGxr, xlface, xrface
  integer,  intent(in) :: BCYl, BCYr, x1s, x1e, Ny, NGyl, NGyr, ylface, yrface
  integer,  intent(in) :: BCZl, BCZr, x2s, x2e, Nz, NGzl, NGzr, zlface, zrface
  integer, intent(out) :: ier
  REALSUB, intent(in)  :: a
  real, intent(in) :: dx, dy, dz, dt, theta, sHI, sHeI, sHeII
  real, intent(in) :: lUn, rUn, nUn, nUn0
  real, intent(in) :: BCValsXl(Ny,Nz), BCValsXr(Ny,Nz)
  real, intent(in) :: BCValsYl(Nx,Nz), BCValsYr(Nx,Nz)
  real, intent(in) :: BCValsZl(Nx,Ny), BCValsZr(Nx,Ny)
  real, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr), &
       intent(in), target :: E, E0, HI, HI0, HeI, HeI0, HeII, HeII0
  real, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr), &
       intent(in) :: adj1, adj2, adj3
  real*8, intent(out), dimension(7,x0s:x0e,x1s:x1e,x2s:x2e) :: mat1
  real*8, intent(out), dimension(x0s:x0e,x1s:x1e,x2s:x2e) :: mat2
  real*8, intent(out), dimension(x0s:x0e,x1s:x1e,x2s:x2e) :: mat3
  real*8 :: rhs(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr)

  !--------------
  ! locals
  integer :: i, j, k
  real :: dtfac, dtfac0, c, pi, dxi, dyi, dzi, dxfac, dyfac, dzfac
  real :: AGradEr, Erf, kap, Eavg, R, Rmin, nUnit, Dlim
  real, pointer :: Er(:,:,:), nHI(:,:,:), nHeI(:,:,:), nHeII(:,:,:)


  !=======================================================================
  
  ! initialize output flag, and set mat to have all zero values
  ier = 1
  mat1 = 0.d0
  mat2 = 0.d0
  mat3 = 0.d0

  ! set shortcut values
  dxi   = a/dx/lUn
  dyi   = a/dy/lUn
  dzi   = a/dz/lUn
  dxfac = dt*theta*dxi*dxi
  dyfac = dt*theta*dyi*dyi
  dzfac = dt*theta*dzi*dzi
  c  = 2.99792458d10     ! speed of light [cm/s]
  pi = 4.d0*datan(1.d0)
  Rmin = 1.0e-20

  ! set pointers for Ef1 evaluation (and appropriate limiter construction)
  select case (LImp)
  case(1)       ! lagged to previous newton iterate
     Er => E
     nHI => HI
     nHeI => HeI
     nHeII => HeII
     nUnit = nUn
  case(2)       ! lag only opacity to previous Newton iterate
     Er => E
     nHI => HI
     nHeI => HeI
     nHeII => HeII
     nUnit = nUn
  case default  ! fully lagged to previous time step
     Er => E0
     nHI => HI0
     nHeI => HeI0
     nHeII => HeII0
     nUnit = nUn0
  end select



!!$  print *, 'entering MFProb_SetupSystem_3D:'
!!$  print *, '    LType =',LType
!!$  print *, '    LImp =',LImp
!!$  print *, '    BCx* =',BCXl,BCXr
!!$  print *, '    BCy* =',BCYl,BCYr
!!$  print *, '    BCz* =',BCZl,BCZr
!!$  print *, '    x0* =',x0s,x0e
!!$  print *, '    x1* =',x1s,x1e
!!$  print *, '    x2* =',x2s,x2e
!!$  print *, '    Nx* =',Nx,NGxl,NGxr
!!$  print *, '    Ny* =',Ny,NGyl,NGyr
!!$  print *, '    Nz* =',Nz,NGzl,NGzr
!!$  print *, '    x*face =',xlface,xrface
!!$  print *, '    y*face =',ylface,yrface
!!$  print *, '    z*face =',zlface,zrface
!!$  print *, '    a =',a
!!$  print *, '    dx* =',dx,dy,dz
!!$  print *, '    dt =',dt
!!$  print *, '    theta =',theta
!!$  print *, '    sH* =',sHI,sHeI,sHeII
!!$  print *, '    units =',lUn,rUn,nUn,nUn0
!!$  print *, '    BCValsX* =',sum(BCValsXl),sum(BCValsXr)
!!$  print *, '    BCValsY* =',sum(BCValsYl),sum(BCValsYr)
!!$  print *, '    BCValsZ* =',sum(BCValsZl),sum(BCValsZr)
!!$  print *, '    E,E0 =',sum(E),sum(E0)
!!$  print *, '    HI,HI0 =',sum(HI),sum(HI0)
!!$  print *, '    HeI,HeI0 =',sum(HeI),sum(HeI0)
!!$  print *, '    HeII,HeII0 =',sum(HeII),sum(HeII0)
!!$  print *, '    adj* =',sum(adj1),sum(adj2),sum(adj3)


  ! iterate over the active domain
  do k=1,Nz,1
     do j=1,Ny,1
        do i=1,Nx,1

           ! initialize matrix entries
           mat1(:,i,j,k) = 0.d0
           mat1(4,i,j,k) = adj1(i,j,k)
           mat2(i,j,k)   = adj2(i,j,k)
           mat3(i,j,k)   = adj3(i,j,k)

           !--------------
           ! z-directional limiter, lower face
           !    face-centered gradient, value, opacity
           AGradEr = abs(Er(i,j,k) - Er(i,j,k-1))*dzi
           Erf = (Er(i,j,k) + Er(i,j,k-1))/2.d0

           !    compute R for limiter 
           R = max(AGradEr/Erf, Rmin)

           !    compute opacity
           kap = (sHI*(  nHI(i,j,k)   + nHI(i,j,k-1))  &
                + sHeI*( nHeI(i,j,k)  + nHeI(i,j,k-1)) &
                + sHeII*(nHeII(i,j,k) + nHeII(i,j,k-1)))*0.5d0*nUnit

           !    compute limiter
           if (LType == 1) then       ! rational approx. to LP lim. (LP, 1981)
              Dlim = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
           else if (LType == 2) then  ! Reynolds approx to LP lim.
              Dlim = 2.d0*c/pi*datan(R*pi/6.d0/kap)/R
           else if (LType == 3) then  ! no limiter
              Dlim = c/kap/3.d0
           else if (LType == 4) then  ! Zeus limiter
              Dlim = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
           else                       ! standard Levermore-Pomraning (LP, 1981)
              Dlim = c*(cosh(R/kap)/sinh(R/kap)-kap/R)/R
           endif

           !    set the relevant matrix entries. Note: the diffusive component 
           !    need not be rescaled, since scaling and chain rule cancel 
           !       dep. on z-left Er
           mat1(1,i,j,k) = mat1(1,i,j,k) - dzfac*Dlim
           !       dep. on self Er
           mat1(4,i,j,k) = mat1(4,i,j,k) + dzfac*Dlim


           !--------------
           ! y-directional limiter, lower face
           !    face-centered gradient, value, opacity
           AGradEr = abs(Er(i,j,k) - Er(i,j-1,k))*dyi
           Erf = (Er(i,j,k) + Er(i,j-1,k))/2.d0

           !    compute R for limiter 
           R = max(AGradEr/Erf, Rmin)

           !    compute opacity
           kap = (sHI*(  nHI(i,j,k)   + nHI(i,j-1,k))  &
                + sHeI*( nHeI(i,j,k)  + nHeI(i,j-1,k)) &
                + sHeII*(nHeII(i,j,k) + nHeII(i,j-1,k)))*0.5d0*nUnit

           !    compute limiter
           if (LType == 1) then       ! rational approx. to LP lim. (LP, 1981)
              Dlim = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
           else if (LType == 2) then  ! Reynolds approx to LP lim.
              Dlim = 2.d0*c/pi*datan(R*pi/6.d0/kap)/R
           else if (LType == 3) then  ! no limiter
              Dlim = c/kap/3.d0
           else if (LType == 4) then  ! Zeus limiter
              Dlim = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
           else                       ! standard Levermore-Pomraning (LP, 1981)
              Dlim = c*(cosh(R/kap)/sinh(R/kap)-kap/R)/R
           endif

           !    set the relevant matrix entries. Note: the diffusive component 
           !    need not be rescaled, since scaling and chain rule cancel 
           !       dep. on y-left Er
           mat1(2,i,j,k) = mat1(2,i,j,k) - dyfac*Dlim
           !       dep. on self Er
           mat1(4,i,j,k) = mat1(4,i,j,k) + dyfac*Dlim


           !--------------
           ! x-directional limiter, lower face
           !    face-centered gradient, value, opacity
           AGradEr = abs(Er(i,j,k) - Er(i-1,j,k))*dxi
           Erf = (Er(i,j,k) + Er(i-1,j,k))/2.d0

           !    compute R for limiter 
           R = max(AGradEr/Erf, Rmin)

           !    compute opacity
           kap = (sHI*(  nHI(i,j,k)   + nHI(i-1,j,k))  &
                + sHeI*( nHeI(i,j,k)  + nHeI(i-1,j,k)) &
                + sHeII*(nHeII(i,j,k) + nHeII(i-1,j,k)))*0.5d0*nUnit

           !    compute limiter
           if (LType == 1) then       ! rational approx. to LP lim. (LP, 1981)
              Dlim = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
           else if (LType == 2) then  ! Reynolds approx to LP lim.
              Dlim = 2.d0*c/pi*datan(R*pi/6.d0/kap)/R
           else if (LType == 3) then  ! no limiter
              Dlim = c/kap/3.d0
           else if (LType == 4) then  ! Zeus limiter
              Dlim = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
           else                       ! standard Levermore-Pomraning (LP, 1981)
              Dlim = c*(cosh(R/kap)/sinh(R/kap)-kap/R)/R
           endif

           !    set the relevant matrix entries. Note: the diffusive component 
           !    need not be rescaled, since scaling and chain rule cancel 
           !       dep. on x-left Er
           mat1(3,i,j,k) = mat1(3,i,j,k) - dxfac*Dlim
           !       dep. on self Er
           mat1(4,i,j,k) = mat1(4,i,j,k) + dxfac*Dlim


           !--------------
           ! x-directional limiter, upper face
           !    face-centered gradient, value, opacity
           AGradEr = abs(Er(i+1,j,k) - Er(i,j,k))*dxi
           Erf = (Er(i,j,k) + Er(i+1,j,k))/2.d0

           !    compute R for limiter 
           R = max(AGradEr/Erf, Rmin)

           !    compute opacity
           kap = (sHI*(  nHI(i,j,k)   + nHI(i+1,j,k))  &
                + sHeI*( nHeI(i,j,k)  + nHeI(i+1,j,k)) &
                + sHeII*(nHeII(i,j,k) + nHeII(i+1,j,k)))*0.5d0*nUnit
           
           !    compute limiter
           if (LType == 1) then       ! rational approx. to LP lim. (LP, 1981)
              Dlim = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
           else if (LType == 2) then  ! Reynolds approx to LP lim.
              Dlim = 2.d0*c/pi*datan(R*pi/6.d0/kap)/R
           else if (LType == 3) then  ! no limiter
              Dlim = c/kap/3.d0
           else if (LType == 4) then  ! Zeus limiter
              Dlim = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
           else                       ! standard Levermore-Pomraning (LP, 1981)
              Dlim = c*(cosh(R/kap)/sinh(R/kap)-kap/R)/R
           endif

           !    set the relevant matrix entries. Note: the diffusive component 
           !    need not be rescaled, since scaling and chain rule cancel 
           !       dep. on x-right Er
           mat1(5,i,j,k) = mat1(5,i,j,k) - dxfac*Dlim
           !       dep. on self Er
           mat1(4,i,j,k) = mat1(4,i,j,k) + dxfac*Dlim


           !--------------
           ! y-directional limiter, upper face
           !    face-centered gradient, value, opacity
           AGradEr = abs(Er(i,j+1,k) - Er(i,j,k))*dyi
           Erf = (Er(i,j,k) + Er(i,j+1,k))/2.d0

           !    compute R for limiter 
           R = max(AGradEr/Erf, Rmin)

           kap = (sHI*(  nHI(i,j,k)   + nHI(i,j+1,k))  &
                + sHeI*( nHeI(i,j,k)  + nHeI(i,j+1,k)) &
                + sHeII*(nHeII(i,j,k) + nHeII(i,j+1,k)))*0.5d0*nUnit

           !    compute limiter
           if (LType == 1) then       ! rational approx. to LP lim. (LP, 1981)
              Dlim = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
           else if (LType == 2) then  ! Reynolds approx to LP lim.
              Dlim = 2.d0*c/pi*datan(R*pi/6.d0/kap)/R
           else if (LType == 3) then  ! no limiter
              Dlim = c/kap/3.d0
           else if (LType == 4) then  ! Zeus limiter
              Dlim = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
           else                       ! standard Levermore-Pomraning (LP, 1981)
              Dlim = c*(cosh(R/kap)/sinh(R/kap)-kap/R)/R
           endif

           !    set the relevant matrix entries. Note: the diffusive component 
           !    need not be rescaled, since scaling and chain rule cancel 
           !       dep. on y-right Er
           mat1(6,i,j,k) = mat1(6,i,j,k) - dyfac*Dlim
           !       dep. on self Er
           mat1(4,i,j,k) = mat1(4,i,j,k) + dyfac*Dlim


           !--------------
           ! z-directional limiter, upper face
           !    face-centered gradient, value, opacity
           AGradEr = abs(Er(i,j,k+1) - Er(i,j,k))*dzi
           Erf = (Er(i,j,k) + Er(i,j,k+1))/2.d0

           !    compute R for limiter 
           R = max(AGradEr/Erf, Rmin)

           !    compute opacity
           kap = (sHI*(  nHI(i,j,k)   + nHI(i,j+1,k))  &
                + sHeI*( nHeI(i,j,k)  + nHeI(i,j+1,k)) &
                + sHeII*(nHeII(i,j,k) + nHeII(i,j+1,k)))*0.5d0*nUnit

           !    compute limiter
           if (LType == 1) then       ! rational approx. to LP lim. (LP, 1981)
              Dlim = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R) 
           else if (LType == 2) then  ! Reynolds approx to LP lim.
              Dlim = 2.d0*c/pi*datan(R*pi/6.d0/kap)/R
           else if (LType == 3) then  ! no limiter
              Dlim = c/kap/3.d0
           else if (LType == 4) then  ! Zeus limiter
              Dlim = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R) 
           else                       ! standard Levermore-Pomraning (LP, 1981)
              Dlim = c*(cosh(R/kap)/sinh(R/kap)-kap/R)/R
           endif

           !    set the relevant matrix entries. Note: the diffusive component 
           !    need not be rescaled, since scaling and chain rule cancel 
           !       dep. on z-right Er
           mat1(7,i,j,k) = mat1(7,i,j,k) - dzfac*Dlim
           !       dep. on self Er
           mat1(4,i,j,k) = mat1(4,i,j,k) + dzfac*Dlim

        enddo
     enddo
  enddo

!!$  print *,'  sum(mat1) =',sum(mat1)
!!$  print *,'  sum(mat2) =',sum(mat2)
!!$  print *,'  sum(mat3) =',sum(mat3)
!!$  print *,'  sum(rhs) =',sum(rhs)



  ! update matrix/rhs based on boundary conditions/location
  !    z-left face
  if (zlface == 1) then
     ! Dirichlet
     if (BCZl==1) then
        k = 0
        do j=1,Ny
           do i=1,Nx
              mat1(:,i,j,k) = 0.d0
              mat1(4,i,j,k) = 1.d0
              rhs(i,j,k)    = BCvalsZl(i,j)/rUn
           enddo
        enddo
     ! Neumann
     else if (BCZl==2) then
        k = 1
        do j=1,Ny
           do i=1,Nx
              R = mat1(1,i,j,k)
              mat1(1,i,j,k) = 0.d0
              mat1(4,i,j,k) = mat1(4,i,j,k) + R
              rhs(i,j,k)    = rhs(i,j,k) - R*BCvalsZl(i,j)/rUn/dzi
           enddo
        enddo
     endif
  end if

  !    y-left face
  if (ylface == 1) then
     ! Dirichlet
     if (BCYl==1) then
        j = 0
        do k=1,Nz
           do i=1,Nx
              mat1(:,i,j,k) = 0.d0
              mat1(4,i,j,k) = 1.d0
              rhs(i,j,k)    = BCvalsYl(i,k)/rUn
           enddo
        enddo
     ! Neumann
     else if (BCYl==2) then
        j = 1
        do k=1,Nz
           do i=1,Nx
              R = mat1(2,i,j,k)
              mat1(2,i,j,k) = 0.d0
              mat1(4,i,j,k) = mat1(4,i,j,k) + R
              rhs(i,j,k)    = rhs(i,j,k) - R*BCvalsYl(i,k)/rUn/dyi
           enddo
        enddo
     endif
  end if

  !    x-left face
  if (xlface == 1) then
     ! Dirichlet
     if (BCXl==1) then
        i = 0
        do k=1,Nz
           do j=1,Ny
              mat1(:,i,j,k) = 0.d0
              mat1(4,i,j,k) = 1.d0
              rhs(i,j,k)    = BCvalsXl(j,k)/rUn
           enddo
        enddo
     ! Neumann
     else if (BCXl==2) then
        i = 1
        do k=1,Nz
           do j=1,Ny
              R = mat1(3,i,j,k)
              mat1(3,i,j,k) = 0.d0
              mat1(4,i,j,k) = mat1(4,i,j,k) + R
              rhs(i,j,k)    = rhs(i,j,k) - R*BCvalsXl(j,k)/rUn/dxi
           enddo
        enddo
     endif
  end if

  !    x-right face
  if (xrface==1) then
     ! Dirichlet
     if (BCXr==1) then
        i = Nx+1
        do k=1,Nz
           do j=1,Ny
              mat1(:,i,j,k) = 0.d0
              mat1(4,i,j,k) = 1.d0
              rhs(i,j,k)    = BCvalsXr(j,k)/rUn
           enddo
        enddo
     ! Neumann
     else if (BCXr==2) then
        i = Nx
        do k=1,Nz
           do j=1,Ny
              R = mat1(5,i,j,k)
              mat1(4,i,j,k) = mat1(4,i,j,k) + R
              mat1(5,i,j,k) = 0.d0
              rhs(i,j,k)    = rhs(i,j,k) - R*BCvalsXr(j,k)/rUn/dxi
           enddo
        enddo
     endif
  endif

  !    y-right face
  if (yrface==1) then
     ! Dirichlet
     if (BCYr==1) then
        j = Ny+1
        do k=1,Nz
           do i=1,Nx
              mat1(:,i,j,k) = 0.d0
              mat1(4,i,j,k) = 1.d0
              rhs(i,j,k)    = BCvalsYr(i,k)/rUn
           enddo
        enddo
     ! Neumann
     else if (BCYr==2) then
        j = Ny
        do k=1,Nz
           do i=1,Nx
              R = mat1(6,i,j,k)
              mat1(4,i,j,k) = mat1(4,i,j,k) + R
              mat1(6,i,j,k) = 0.d0
              rhs(i,j,k)    = rhs(i,j,k) - R*BCvalsYr(i,k)/rUn/dyi
           enddo
        enddo
     endif
  endif

  !    z-right face
  if (zrface==1) then
     ! Dirichlet
     if (BCZr==1) then
        k = Nz+1
        do j=1,Ny
           do i=1,Nx
              mat1(:,i,j,k) = 0.d0
              mat1(4,i,j,k) = 1.d0
              rhs(i,j,k)    = BCvalsZr(i,j)/rUn
           enddo
        enddo
     ! Neumann
     else if (BCZr==2) then
        k = Nz
        do j=1,Ny
           do i=1,Nx
              R = mat1(7,i,j,k)
              mat1(4,i,j,k) = mat1(4,i,j,k) + R
              mat1(7,i,j,k) = 0.d0
              rhs(i,j,k)    = rhs(i,j,k) - R*BCvalsZr(i,j)/rUn/dzi
           enddo
        enddo
     endif
  endif


!!$  print *,'  sum(mat1) =',sum(mat1)
!!$  print *,'  sum(rhs) =',sum(rhs)


  nullify(Er)
  
  return
end subroutine MFProb_SetupSystem_3D
!=======================================================================






subroutine MFProb_SetupSystem_2D(mat1, mat2, mat3, rhs, E, E0, HI, HI0,     &
     HeI, HeI0, HeII, HeII0, adj1, adj2, adj3, LType, LImp, dt, theta, sHI, &
     sHeI, sHeII, a, lUn, rUn, nUn, nUn0, dx, dy, BCValsXl, BCValsXr,       &
     BCValsYl, BCValsYr, BCXl, BCXr, BCYl, BCYr, x0s, x0e, x1s, x1e, Nx,    &
     Ny, NGxl, NGxr, NGyl, NGyr, xlface, xrface, ylface, yrface, ier)
  !=======================================================================
  !  PURPOSE: 2D version of the routine
  !=======================================================================
#include "fortran.def"
  implicit none
  
  !--------------
  ! argument declarations
  integer, intent(in)  :: LType, LImp
  integer,  intent(in) :: BCXl, BCXr, x0s, x0e, Nx, NGxl, NGxr, xlface, xrface
  integer,  intent(in) :: BCYl, BCYr, x1s, x1e, Ny, NGyl, NGyr, ylface, yrface
  integer, intent(out) :: ier
  REALSUB, intent(in)  :: a
  real, intent(in) :: dx, dy, dt, theta, sHI, sHeI, sHeII
  real, intent(in) :: lUn, rUn, nUn, nUn0
  real, intent(in) :: BCValsXl(Ny), BCValsXr(Ny)
  real, intent(in) :: BCValsYl(Nx), BCValsYr(Nx)
  real, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr), &
       intent(in), target :: E, E0, HI, HI0, HeI, HeI0, HeII, HeII0
  real, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr), &
       intent(in) :: adj1, adj2, adj3
  real*8, intent(out), dimension(5,x0s:x0e,x1s:x1e) :: mat1
  real*8, intent(out), dimension(x0s:x0e,x1s:x1e) :: mat2
  real*8, intent(out), dimension(x0s:x0e,x1s:x1e) :: mat3
  real*8 :: rhs(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr)

  !--------------
  ! locals
  integer :: i, j
  real :: dtfac, dtfac0, c, pi, dxi, dyi, dxfac, dyfac
  real :: AGradEr, Erf, kap, Eavg, R, Rmin, nUnit, Dlim
  real, pointer :: Er(:,:), nHI(:,:), nHeI(:,:), nHeII(:,:)


  !=======================================================================
  
  ! initialize output flag, and set mat to have all zero values
  ier = 1
  mat1 = 0.d0
  mat2 = 0.d0
  mat3 = 0.d0

  ! set shortcut values
  dxi   = a/dx/lUn
  dyi   = a/dy/lUn
  dxfac = dt*theta*dxi*dxi
  dyfac = dt*theta*dyi*dyi
  c  = 2.99792458d10     ! speed of light [cm/s]
  pi = 4.d0*datan(1.d0)
  Rmin = 1.0e-20

  ! set pointers for Ef1 evaluation (and appropriate limiter construction)
  select case (LImp)
  case(1)       ! lagged to previous newton iterate
     Er => E
     nHI => HI
     nHeI => HeI
     nHeII => HeII
     nUnit = nUn
  case(2)       ! lag only opacity to previous Newton iterate
     Er => E
     nHI => HI
     nHeI => HeI
     nHeII => HeII
     nUnit = nUn
  case default  ! fully lagged to previous time step
     Er => E0
     nHI => HI0
     nHeI => HeI0
     nHeII => HeII0
     nUnit = nUn0
  end select

  
  ! iterate over the active domain
  do j=1,Ny,1
     do i=1,Nx,1

        ! initialize matrix entries
        mat1(:,i,j) = 0.d0
        mat1(4,i,j) = adj1(i,j)
        mat2(i,j)   = adj2(i,j)
        mat3(i,j)   = adj3(i,j)

        !--------------
        ! y-directional limiter, lower face
        !    face-centered gradient, value, opacity
        AGradEr = abs(Er(i,j) - Er(i,j-1))*dyi
        Erf = (Er(i,j) + Er(i,j-1))/2.d0

        !    compute R for limiter 
        R = max(AGradEr/Erf, Rmin)

        !    compute opacity
        kap = (sHI*(  nHI(i,j)   + nHI(i,j-1))  &
             + sHeI*( nHeI(i,j)  + nHeI(i,j-1)) &
             + sHeII*(nHeII(i,j) + nHeII(i,j-1)))*0.5d0*nUnit

        !    compute limiter
        if (LType == 1) then       ! rational approx. to LP lim. (LP, 1981)
           Dlim = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
        else if (LType == 2) then  ! Reynolds approx to LP lim.
           Dlim = 2.d0*c/pi*datan(R*pi/6.d0/kap)/R
        else if (LType == 3) then  ! no limiter
           Dlim = c/kap/3.d0
        else if (LType == 4) then  ! Zeus limiter
           Dlim = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
        else                       ! standard Levermore-Pomraning (LP, 1981)
           Dlim = c*(cosh(R/kap)/sinh(R/kap)-kap/R)/R
        endif

        !    set the relevant matrix entries. Note: the diffusive component 
        !    need not be rescaled, since scaling and chain rule cancel 
        !       dep. on y-left Er
        mat1(1,i,j) = mat1(1,i,j) - dyfac*Dlim
        !       dep. on self Er
        mat1(3,i,j) = mat1(3,i,j) + dyfac*Dlim


        !--------------
        ! x-directional limiter, lower face
        !    face-centered gradient, value, opacity
        AGradEr = abs(Er(i,j) - Er(i-1,j))*dxi
        Erf = (Er(i,j) + Er(i-1,j))/2.d0

        !    compute R for limiter 
        R = max(AGradEr/Erf, Rmin)

        !    compute opacity
        kap = (sHI*(  nHI(i,j)   + nHI(i-1,j))  &
             + sHeI*( nHeI(i,j)  + nHeI(i-1,j)) &
             + sHeII*(nHeII(i,j) + nHeII(i-1,j)))*0.5d0*nUnit

        !    compute limiter
        if (LType == 1) then       ! rational approx. to LP lim. (LP, 1981)
           Dlim = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
        else if (LType == 2) then  ! Reynolds approx to LP lim.
           Dlim = 2.d0*c/pi*datan(R*pi/6.d0/kap)/R
        else if (LType == 3) then  ! no limiter
           Dlim = c/kap/3.d0
        else if (LType == 4) then  ! Zeus limiter
           Dlim = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
        else                       ! standard Levermore-Pomraning (LP, 1981)
           Dlim = c*(cosh(R/kap)/sinh(R/kap)-kap/R)/R
        endif

        !    set the relevant matrix entries. Note: the diffusive component 
        !    need not be rescaled, since scaling and chain rule cancel 
        !       dep. on x-left Er
        mat1(2,i,j) = mat1(2,i,j) - dxfac*Dlim
        !       dep. on self Er
        mat1(3,i,j) = mat1(3,i,j) + dxfac*Dlim


        !--------------
        ! x-directional limiter, upper face
        !    face-centered gradient, value, opacity
        AGradEr = abs(Er(i+1,j) - Er(i,j))*dxi
        Erf = (Er(i,j) + Er(i+1,j))/2.d0

        !    compute R for limiter 
        R = max(AGradEr/Erf, Rmin)

        !    compute opacity
        kap = (sHI*(  nHI(i,j)   + nHI(i+1,j))  &
             + sHeI*( nHeI(i,j)  + nHeI(i+1,j)) &
             + sHeII*(nHeII(i,j) + nHeII(i+1,j)))*0.5d0*nUnit

        !    compute limiter
        if (LType == 1) then       ! rational approx. to LP lim. (LP, 1981)
           Dlim = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
        else if (LType == 2) then  ! Reynolds approx to LP lim.
           Dlim = 2.d0*c/pi*datan(R*pi/6.d0/kap)/R
        else if (LType == 3) then  ! no limiter
           Dlim = c/kap/3.d0
        else if (LType == 4) then  ! Zeus limiter
           Dlim = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
        else                       ! standard Levermore-Pomraning (LP, 1981)
           Dlim = c*(cosh(R/kap)/sinh(R/kap)-kap/R)/R
        endif

        !    set the relevant matrix entries. Note: the diffusive component 
        !    need not be rescaled, since scaling and chain rule cancel 
        !       dep. on x-right Er
        mat1(4,i,j) = mat1(4,i,j) - dxfac*Dlim
        !       dep. on self Er
        mat1(3,i,j) = mat1(3,i,j) + dxfac*Dlim


        !--------------
        ! y-directional limiter, upper face
        !    face-centered gradient, value, opacity
        AGradEr = abs(Er(i,j+1) - Er(i,j))*dyi
        Erf = (Er(i,j) + Er(i,j+1))/2.d0

        !    compute R for limiter 
        R = max(AGradEr/Erf, Rmin)

        kap = (sHI*(  nHI(i,j)   + nHI(i,j+1))  &
             + sHeI*( nHeI(i,j)  + nHeI(i,j+1)) &
             + sHeII*(nHeII(i,j) + nHeII(i,j+1)))*0.5d0*nUnit

        !    compute limiter
        if (LType == 1) then       ! rational approx. to LP lim. (LP, 1981)
           Dlim = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
        else if (LType == 2) then  ! Reynolds approx to LP lim.
           Dlim = 2.d0*c/pi*datan(R*pi/6.d0/kap)/R
        else if (LType == 3) then  ! no limiter
           Dlim = c/kap/3.d0
        else if (LType == 4) then  ! Zeus limiter
           Dlim = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
        else                       ! standard Levermore-Pomraning (LP, 1981)
           Dlim = c*(cosh(R/kap)/sinh(R/kap)-kap/R)/R
        endif

        !    set the relevant matrix entries. Note: the diffusive component 
        !    need not be rescaled, since scaling and chain rule cancel 
        !       dep. on y-right Er
        mat1(5,i,j) = mat1(5,i,j) - dyfac*Dlim
        !       dep. on self Er
        mat1(3,i,j) = mat1(3,i,j) + dyfac*Dlim

     enddo
  enddo


!!$  print *,'  sum(mat1) =',sum(mat1)
!!$  print *,'  sum(mat2) =',sum(mat2)
!!$  print *,'  sum(mat3) =',sum(mat3)
!!$  print *,'  sum(rhs) =',sum(rhs)


  ! update matrix/rhs based on boundary conditions/location
  !    y-left face
  if (ylface == 1) then
     ! Dirichlet
     if (BCYl==1) then
        j = 0
        do i=1,Nx
           mat1(:,i,j) = 0.d0
           mat1(3,i,j) = 1.d0
           rhs(i,j)    = BCvalsYl(i)/rUn
        enddo
     ! Neumann
     else if (BCYl==2) then
        j = 1
        do i=1,Nx
           R = mat1(1,i,j)
           mat1(1,i,j) = 0.d0
           mat1(3,i,j) = mat1(3,i,j) + R
           rhs(i,j)    = rhs(i,j) - R*BCvalsYl(i)/rUn/dyi
        enddo
     endif
  end if

  !    x-left face
  if (xlface == 1) then
     ! Dirichlet
     if (BCXl==1) then
        i = 0
        do j=1,Ny
           mat1(:,i,j) = 0.d0
           mat1(3,i,j) = 1.d0
           rhs(i,j)    = BCvalsXl(j)/rUn
        enddo
     ! Neumann
     else if (BCXl==2) then
        i = 1
        do j=1,Ny
           R = mat1(2,i,j)
           mat1(2,i,j) = 0.d0
           mat1(3,i,j) = mat1(3,i,j) + R
           rhs(i,j)    = rhs(i,j) - R*BCvalsXl(j)/rUn/dxi
        enddo
     endif
  end if

  !    x-right face
  if (xrface==1) then
     ! Dirichlet
     if (BCXr==1) then
        i = Nx+1
        do j=1,Ny
           mat1(:,i,j) = 0.d0
           mat1(3,i,j) = 1.d0
           rhs(i,j)    = BCvalsXr(j)/rUn
        enddo
     ! Neumann
     else if (BCXr==2) then
        i = Nx
        do j=1,Ny
           R = mat1(4,i,j)
           mat1(3,i,j) = mat1(3,i,j) + R
           mat1(4,i,j) = 0.d0
           rhs(i,j)    = rhs(i,j) - R*BCvalsXr(j)/rUn/dxi
        enddo
     endif
  endif

  !    y-right face
  if (yrface==1) then
     ! Dirichlet
     if (BCYr==1) then
        j = Ny+1
        do i=1,Nx
           mat1(:,i,j) = 0.d0
           mat1(3,i,j) = 1.d0
           rhs(i,j)    = BCvalsYr(i)/rUn
        enddo
     ! Neumann
     else if (BCYr==2) then
        j = Ny
        do i=1,Nx
           R = mat1(5,i,j)
           mat1(3,i,j) = mat1(3,i,j) + R
           mat1(5,i,j) = 0.d0
           rhs(i,j)    = rhs(i,j) - R*BCvalsYr(i)/rUn/dyi
        enddo
     endif
  endif


!!$  print *,'  sum(mat1) =',sum(mat1)
!!$  print *,'  sum(rhs) =',sum(rhs)


  nullify(Er)
  
  return
end subroutine MFProb_SetupSystem_2D
!=======================================================================






subroutine MFProb_SetupSystem_1D(mat1, mat2, mat3, rhs, E, E0, HI, HI0,   &
     HeI, HeI0, HeII, HeII0, adj1, adj2, adj3, LType, LImp, dt, theta,    &
     sHI, sHeI, sHeII, a, lUn, rUn, nUn, nUn0, dx, BCValsXl, BCValsXr,    &
     BCXl, BCXr, x0s, x0e, Nx, NGxl, NGxr, xlface, xrface, ier)
  !=======================================================================
  !  PURPOSE: 1D version of the routine
  !=======================================================================
#include "fortran.def"
  implicit none
  
  !--------------
  ! argument declarations
  integer, intent(in)  :: LType, LImp
  integer,  intent(in) :: BCXl, BCXr, x0s, x0e, Nx, NGxl, NGxr, xlface, xrface
  integer, intent(out) :: ier
  REALSUB, intent(in)  :: a
  real, intent(in) :: dx, dt, theta, sHI, sHeI, sHeII
  real, intent(in) :: lUn, rUn, nUn, nUn0
  real, intent(in) :: BCValsXl(1), BCValsXr(1)
  real, dimension(1-NGxl:Nx+NGxr), intent(in), target :: E, E0, HI, HI0, &
       HeI, HeI0, HeII, HeII0
  real, dimension(1-NGxl:Nx+NGxr), intent(in) :: adj1, adj2, adj3
  real*8, intent(out), dimension(3,x0s:x0e) :: mat1
  real*8, intent(out), dimension(x0s:x0e) :: mat2
  real*8, intent(out), dimension(x0s:x0e) :: mat3
  real*8 :: rhs(1-NGxl:Nx+NGxr)

  !--------------
  ! locals
  integer :: i
  real :: dtfac, dtfac0, c, pi, dxi, dxfac
  real :: AGradEr, Erf, kap, Eavg, R, Rmin, nUnit, Dlim
  real, pointer :: Er(:), nHI(:), nHeI(:), nHeII(:)


  !=======================================================================
  
  ! initialize output flag, and set mat to have all zero values
  ier = 1
  mat1 = 0.d0
  mat2 = 0.d0
  mat3 = 0.d0

  ! set shortcut values
  dxi   = a/dx/lUn
  dxfac = dt*theta*dxi*dxi
  c  = 2.99792458d10     ! speed of light [cm/s]
  pi = 4.d0*datan(1.d0)
  Rmin = 1.0e-20

  ! set pointers for Ef1 evaluation (and appropriate limiter construction)
  select case (LImp)
  case(1)       ! lagged to previous newton iterate
     Er => E
     nHI => HI
     nHeI => HeI
     nHeII => HeII
     nUnit = nUn
  case(2)       ! lag only opacity to previous Newton iterate
     Er => E
     nHI => HI
     nHeI => HeI
     nHeII => HeII
     nUnit = nUn
  case default  ! fully lagged to previous time step
     Er => E0
     nHI => HI0
     nHeI => HeI0
     nHeII => HeII0
     nUnit = nUn0
  end select


  ! iterate over the active domain
  do i=1,Nx,1

     ! initialize matrix entries
     mat1(:,i) = 0.d0
     mat1(2,i) = adj1(i)
     mat2(i) = adj2(i)
     mat3(i) = adj3(i)

     !--------------
     ! x-directional limiter, lower face
     !    face-centered gradient, value, opacity
     AGradEr = abs(Er(i) - Er(i-1))*dxi
     Erf = (Er(i) + Er(i-1))/2.d0

     !    compute R for limiter 
     R = max(AGradEr/Erf, Rmin)

     !    compute opacity
     kap = (sHI*(  nHI(i)   + nHI(i-1))  &
          + sHeI*( nHeI(i)  + nHeI(i-1)) &
          + sHeII*(nHeII(i) + nHeII(i-1)))*0.5d0*nUnit

     !    compute limiter
     if (LType == 1) then       ! rational approx. to LP lim. (LP, 1981)
        Dlim = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
     else if (LType == 2) then  ! Reynolds approx to LP lim.
        Dlim = 2.d0*c/pi*datan(R*pi/6.d0/kap)/R
     else if (LType == 3) then  ! no limiter
        Dlim = c/kap/3.d0
     else if (LType == 4) then  ! Zeus limiter
        Dlim = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
     else                       ! standard Levermore-Pomraning (LP, 1981)
        Dlim = c*(cosh(R/kap)/sinh(R/kap)-kap/R)/R
     endif

     !    set the relevant matrix entries. Note: the diffusive component 
     !    need not be rescaled, since scaling and chain rule cancel 
     !       dep. on x-left Er
     mat1(1,i) = mat1(1,i) - dxfac*Dlim
     !       dep. on self Er
     mat1(2,i) = mat1(2,i) + dxfac*Dlim


     !--------------
     ! x-directional limiter, upper face
     !    face-centered gradient, value, opacity
     AGradEr = abs(Er(i+1) - Er(i))*dxi
     Erf = (Er(i) + Er(i+1))/2.d0

     !    compute R for limiter 
     R = max(AGradEr/Erf, Rmin)

     !    compute opacity
     kap = (sHI*(  nHI(i)   + nHI(i+1))  &
          + sHeI*( nHeI(i)  + nHeI(i+1)) &
          + sHeII*(nHeII(i) + nHeII(i+1)))*0.5d0*nUnit

     !    compute limiter
     if (LType == 1) then       ! rational approx. to LP lim. (LP, 1981)
        Dlim = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
     else if (LType == 2) then  ! Reynolds approx to LP lim.
        Dlim = 2.d0*c/pi*datan(R*pi/6.d0/kap)/R
     else if (LType == 3) then  ! no limiter
        Dlim = c/kap/3.d0
     else if (LType == 4) then  ! Zeus limiter
        Dlim = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
     else                       ! standard Levermore-Pomraning (LP, 1981)
        Dlim = c*(cosh(R/kap)/sinh(R/kap)-kap/R)/R
     endif

     !    set the relevant matrix entries. Note: the diffusive component 
     !    need not be rescaled, since scaling and chain rule cancel 
     !       dep. on x-right Er
     mat1(3,i) = mat1(3,i) - dxfac*Dlim
     !       dep. on self Er
     mat1(2,i) = mat1(2,i) + dxfac*Dlim

  enddo


!!$  print *,'  sum(mat1) =',sum(mat1)
!!$  print *,'  sum(mat2) =',sum(mat2)
!!$  print *,'  sum(mat3) =',sum(mat3)
!!$  print *,'  sum(rhs) =',sum(rhs)


  ! update matrix/rhs based on boundary conditions/location
  !    x-left face
  if (xlface == 1) then
     ! Dirichlet
     if (BCXl==1) then
        i = 0
        mat1(:,i) = 0.d0
        mat1(2,i) = 1.d0
        rhs(i)    = BCvalsXl(1)/rUn
     ! Neumann
     else if (BCXl==2) then
        i = 1
        R = mat1(1,i)
        mat1(1,i) = 0.d0
        mat1(2,i) = mat1(2,i) + R
        rhs(i)    = rhs(i) - R*BCvalsXl(1)/rUn/dxi
     endif
  end if

  !    x-right face
  if (xrface==1) then
     ! Dirichlet
     if (BCXr==1) then
        i = Nx+1
        mat1(:,i) = 0.d0
        mat1(2,i) = 1.d0
        rhs(i)    = BCvalsXr(1)/rUn
     ! Neumann
     else if (BCXr==2) then
        i = Nx
        R = mat1(3,i)
        mat1(2,i) = mat1(2,i) + R
        mat1(3,i) = 0.d0
        rhs(i)    = rhs(i) - R*BCvalsXr(1)/rUn/dxi
     endif
  endif


!!$  print *,'  sum(mat1) =',sum(mat1)
!!$  print *,'  sum(rhs) =',sum(rhs)


  nullify(Er)
  
  return
end subroutine MFProb_SetupSystem_1D
!=======================================================================
