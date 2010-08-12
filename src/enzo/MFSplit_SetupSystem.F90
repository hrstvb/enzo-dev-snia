!=======================================================================
!
! Copyright 2009 Daniel R. Reynolds
!
! This software is released under the terms of the "Enzo Public License"
! in the accompanying LICENSE file.
!
!=======================================================================
subroutine MFSplit_SetupSystem(mat, rhs, rhsnorm, freq, E, E0, HI, HI0, &
     HeI, HeI0, HeII, HeII0, src, LType, dt, a, a0, theta, aUn, lUn,    &
     lUn0, nUn, nUn0, rank, dx, dy, dz, BCXl, BCXr, BCYl, BCYr, BCZl,   &
     BCZr, x0s, x0e, x1s, x1e, x2s, x2e, Nchem, Nx, Ny, Nz, NGxl, NGxr, &
     NGyl, NGyr, NGzl, NGzr, xlface, xrface, ylface, yrface, zlface, zrface, ier)
  !=======================================================================
  !  written by: Daniel R. Reynolds
  !  date:       December 2009
  !  modified1:  
  !
  !  PURPOSE: Computes the array of matrix stencil elements and vector 
  !           of rhs entries for the multi-frequency problems,
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
  !           As the stencil has {7,5,3} non-zero elements per matrix row
  !           (depending on whether the problem is 3D, 2D or 1D), we 
  !           set these entries over the computational domain, with the 
  !           proper adjustments due to the choice of limiter.
  !
  !           We in fact solve a scaled version of the equation. Moreover, we 
  !           do not solve the equations directly, and instead solve for a 
  !           correction to the current state such that the corrected solutions
  !           satisfy the above equations.  This helps with enforcement of 
  !           boundary conditions, since they may be directly placed onto the 
  !           current state, and the corrections need only refrain from interfering.
  !
  !  INPUTS:
  !     freq       - {1,2,3} monochromatic equation we wish to set up
  !     E,E0       - Radiation energy density (new & old)
  !     HI,HI0     - Hydrogen I density (new & old)
  !     HeI,HeI0   - Helium I density (new & old)
  !     HeII,HeII0 - Helium II density (new & old)
  !     src        - spatially-dependent radiation source
  !     LType      - integer flag denoting type of flux limiter:
  !                       0 -> standard Levermore-Pomraning lim. (LP, 1981)
  !                       1 -> rational approx. to LP lim. (LP, 1981)
  !                       2 -> Reynolds approx to LP lim.
  !                       3 -> turns off the limiter (constant of 1/3)
  !                       4 -> Zeus limiter
  !     dt         - time step size
  !     a,a0       - cosmological expansion parameter (new and old time steps)
  !     theta      - overall implicitness parameter
  !     *Un,*Un0   - variable scaling constants (new and old time steps)
  !     rank       - 1, 2 or 3; the dimensionality of the problem
  !     dx,dy,dz   - mesh spacing in each direction
  !     BC*        - boundary condition type in each direction, face
  !                     0->periodic
  !                     1->Dirichlet
  !                     2->Neumann
  !     x*{s,e}    - start/end indices of linear solver domain; 
  !                  typically 1:Nx for standard dims, but Dirichlet 
  !                  BCs may move these to 0:Nx, 1:Nx+1, etc.
  !     Nchem      - number of chemical species in the problem
  !     Nx,Ny,Nz   - active mesh size in each direction
  !     NG*l/NG*r  - left/right ghost cells in each direction
  !     *{l,r}face - integer flag denoting whether direction/face 
  !                  is external to the domain (0->int, 1->ext)
  !
  !     Note: the vector inputs are of size (Nx + NGxl + NGxr) in 
  !     the x-direction, others are similar.
  !
  !  OUTPUT ARGUMENTS: 
  !     mat   - array of stencil values over the active domain. Since each 
  !                  spatial stencil has 7 nonzero entries, and as this 
  !                  array should not include ghost cells, it has 
  !                  dimensions (7,Nx,Ny,Nz).
  !     rhs    - array of rhs values, same size as variables
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
  integer,  intent(in) :: freq, LType, Nchem, rank
  integer,  intent(in) :: BCXl, BCXr, x0s, x0e, Nx, NGxl, NGxr, xlface, xrface
  integer,  intent(in) :: BCYl, BCYr, x1s, x1e, Ny, NGyl, NGyr, ylface, yrface
  integer,  intent(in) :: BCZl, BCZr, x2s, x2e, Nz, NGzl, NGzr, zlface, zrface
  REALSUB,  intent(in) :: a, a0
  REAL, intent(in) :: dx, dy, dz, dt, theta
  REAL, intent(in) :: aUn, lUn, lUn0, nUn, nUn0
  REAL, intent(in), dimension(*) :: E, E0, HI, HI0, HeI, HeI0, HeII, HeII0
  REAL, intent(in), dimension(*) :: src
  real*8,   intent(out) :: mat(*), rhs(*)
  REAL, intent(out) :: rhsnorm
  integer,  intent(out) :: ier
 
  !--------------
  ! local declarations
  REAL :: hp, ev2erg, nu0_HI, nu0_HeI, nu0_HeII, sHI, sHeI, sHeII

  !=======================================================================

  ! set shortcut values 
  hp = 6.6260693d-27            ! Planck's constant (ergs*s)
  ev2erg = 1.60217653d-12       ! conversion constant from eV to ergs
  nu0_HI   = 13.6d0*ev2erg/hp   ! ionization threshold of HI (hz)
  nu0_HeI  = 24.6d0*ev2erg/hp   ! ionization threshold of HeI (hz)
  nu0_HeII = 54.4d0*ev2erg/hp   ! ionization threshold of HeII (hz)

  ! pass arguments based on Nchem
  if (Nchem == 1) then
  
     !   set shortcut variables depending on which matrix we are setting up
     if (freq == 1) then
        call MFSplit_HICrossSection(sHI, nu0_HI, 1)
        sHeI  = 0.d0
        sHeII = 0.d0
     elseif (freq == 2) then
        call MFSplit_HICrossSection(sHI, nu0_HeI, 1)
        sHeI  = 0.d0
        sHeII = 0.d0
     elseif (freq == 3) then
        call MFSplit_HICrossSection(sHI, nu0_HeII, 1)
        sHeI  = 0.d0
        sHeII = 0.d0
     endif

     ! call the appropriate dimension-specific routine
     if (rank == 1) then
        call MFSplit_SetupSystem_1D(mat, rhs, rhsnorm, E, E0, HI, HI0, HI,  &
             HI0, HI, HI0, src, LType, dt, theta, sHI, sHeI, sHeII, a, a0,  &
             aUn, lUn, lUn0, nUn, nUn0, dx, BCXl, BCXr, x0s, x0e, Nx, NGxl, &
             NGxr, xlface, xrface, ier)
     elseif (rank == 2) then
        call MFSplit_SetupSystem_2D(mat, rhs, rhsnorm, E, E0, HI, HI0, HI,  &
             HI0, HI, HI0, src, LType, dt, theta, sHI, sHeI, sHeII, a, a0,  &
             aUn, lUn, lUn0, nUn, nUn0, dx, dy, BCXl, BCXr, BCYl, BCYr,     &
             x0s, x0e, x1s, x1e, Nx, Ny, NGxl, NGxr, NGyl, NGyr, xlface,    &
             xrface, ylface, yrface, ier)
     else
        call MFSplit_SetupSystem_3D(mat, rhs, rhsnorm, E, E0, HI, HI0, HI,  &
             HI0, HI, HI0, src, LType, dt, theta, sHI, sHeI, sHeII, a, a0,  &
             aUn, lUn, lUn0, nUn, nUn0, dx, dy, dz, BCXl, BCXr, BCYl, BCYr, &
             BCZl, BCZr, x0s, x0e, x1s, x1e, x2s, x2e, Nx, Ny, Nz, NGxl,    &
             NGxr, NGyl, NGyr, NGzl, NGzr, xlface, xrface, ylface, yrface,  &
             zlface, zrface, ier)
     end if
    
  else  ! Nchem == 3

     !   set shortcut variables depending on which matrix we are setting up
     if (freq == 1) then
        call MFSplit_HICrossSection(sHI, nu0_HI, 1)
        sHeI  = 0.d0
        sHeII = 0.d0
     elseif (freq == 2) then
        call MFSplit_HICrossSection(sHI, nu0_HeI, 1)
        call MFSplit_HeICrossSection(sHeI, nu0_HeI, 1)
        sHeII = 0.d0
     elseif (freq == 3) then
        call MFSplit_HICrossSection(sHI, nu0_HeII, 1)
        call MFSplit_HeICrossSection(sHeI, nu0_HeII, 1)
        call MFSplit_HeIICrossSection(sHeII, nu0_HeII, 1)
     else
        write(0,*) 'MFSplit_SetupSystem error, freq =',freq,'is undefined!'
        return
     endif

     ! call the appropriate dimension-specific routine
     if (rank == 1) then
        call MFSplit_SetupSystem_1D(mat, rhs, rhsnorm, E, E0, HI, HI0, HeI, &
             HeI0, HeII, HeII0, src, LType, dt, theta, sHI, sHeI, sHeII, a, &
             a0, aUn, lUn, lUn0, nUn, nUn0, dx, BCXl, BCXr, x0s, x0e, Nx,   &
             NGxl, NGxr, xlface, xrface, ier)
     elseif (rank == 2) then
        call MFSplit_SetupSystem_2D(mat, rhs, rhsnorm, E, E0, HI, HI0, HeI, &
             HeI0, HeII, HeII0, src, LType, dt, theta, sHI, sHeI, sHeII, a, &
             a0, aUn, lUn, lUn0, nUn, nUn0, dx, dy, BCXl, BCXr, BCYl, BCYr, &
             x0s, x0e, x1s, x1e, Nx, Ny, NGxl, NGxr, NGyl, NGyr, xlface,    &
             xrface, ylface, yrface, ier)
     else
        call MFSplit_SetupSystem_3D(mat, rhs, rhsnorm, E, E0, HI, HI0, HeI, &
             HeI0, HeII, HeII0, src, LType, dt, theta, sHI, sHeI, sHeII, a, &
             a0, aUn, lUn, lUn0, nUn, nUn0, dx, dy, dz, BCXl, BCXr, BCYl,   &
             BCYr, BCZl, BCZr, x0s, x0e, x1s, x1e, x2s, x2e, Nx, Ny, Nz,    &
             NGxl, NGxr, NGyl, NGyr, NGzl, NGzr, xlface, xrface, ylface,    &
             yrface, zlface, zrface, ier)
     end if
     
  endif  ! Nchem
  
  return
end subroutine MFSplit_SetupSystem
!=======================================================================






subroutine MFSplit_SetupSystem_3D(mat, rhs, rhsnorm, E, E0, HI, HI0, HeI, &
     HeI0, HeII, HeII0, src, LType, dt, theta, sHI, sHeI, sHeII, a, a0,   &
     aUn, lUn, lUn0, nUn, nUn0, dx, dy, dz, BCXl, BCXr, BCYl, BCYr, BCZl, &
     BCZr, x0s, x0e, x1s, x1e, x2s, x2e, Nx, Ny, Nz, NGxl, NGxr, NGyl,    &
     NGyr, NGzl, NGzr, xlface, xrface, ylface, yrface, zlface, zrface, ier)
  !=======================================================================
  !  PURPOSE: 3D version of the routine
  !=======================================================================
#include "fortran.def"
  implicit none
  
  !--------------
  ! argument declarations
  integer, intent(in)  :: LType
  integer,  intent(in) :: BCXl, BCXr, x0s, x0e, Nx, NGxl, NGxr, xlface, xrface
  integer,  intent(in) :: BCYl, BCYr, x1s, x1e, Ny, NGyl, NGyr, ylface, yrface
  integer,  intent(in) :: BCZl, BCZr, x2s, x2e, Nz, NGzl, NGzr, zlface, zrface
  integer, intent(out) :: ier
  REALSUB, intent(in)  :: a, a0
  REAL, intent(in) :: dx, dy, dz, dt, theta, sHI, sHeI, sHeII
  REAL, intent(in) :: aUn, lUn, lUn0, nUn, nUn0
  REAL, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr), &
       intent(in), target :: E, E0, HI, HI0, HeI, HeI0, HeII, HeII0, src
  real*8, intent(out) :: mat(7,x0s:x0e,x1s:x1e,x2s:x2e)
  real*8, intent(out) :: rhs(x0s:x0e,x1s:x1e,x2s:x2e)
  REAL, intent(out) :: rhsnorm

  !--------------
  ! locals
  integer :: i, j, k
  real*8 :: dtfac, dtfac0, c, pi
  real*8 :: dxi, dxi0, dyi, dyi0, dzi, dzi0
  real*8 :: kap, kap0, E0avg, R, R0, Rmin
  real*8 :: D_xl, D0_xl, D_xr, D0_xr, E0d_xl, E0d_xr, Ed_xl, Ed_xr
  real*8 :: D_yl, D0_yl, D_yr, D0_yr, E0d_yl, E0d_yr, Ed_yl, Ed_yr
  real*8 :: D_zl, D0_zl, D_zr, D0_zr, E0d_zl, E0d_zr, Ed_zl, Ed_zr

  real*8 :: kap_max, kap_avg, nHI_max, nHI_avg

  !=======================================================================
  
  ! initialize outputs to zero, flag to success
  mat = 0.d0
  rhs = 0.d0
  ier = 1

  ! set shortcut values
  dtfac  = dt*theta
  dtfac0 = dt*(1.d0-theta)
  dxi    = a/dx/lUn
  dyi    = a/dy/lUn
  dzi    = a/dz/lUn
  dxi0   = a0/dx/lUn0
  dyi0   = a0/dy/lUn0
  dzi0   = a0/dz/lUn0
  c      = 2.99792458d10     ! speed of light [cm/s]
  pi     = 4.d0*datan(1.d0)
  Rmin   = 1.0d-20

!!$  print *, 'entering MFSplit_SetupSystem_3D:'
!!$  print *, '    LType =',LType
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
!!$  print *, '    units =',lUn,nUn,nUn0
!!$  print *, '    E,E0 =',sum(E),sum(E0)
!!$  print *, '    HI,HI0 =',sum(HI),sum(HI0)
!!$  print *, '    HeI,HeI0 =',sum(HeI),sum(HeI0)
!!$  print *, '    HeII,HeII0 =',sum(HeII),sum(HeII0)
!!$  pause


  kap_max = 0.d0
  kap_avg = 0.d0
  nHI_max = 0.d0
  nHI_avg = 0.d0

  ! iterate over the active domain
  do k=1,Nz,1
     do j=1,Ny,1
        do i=1,Nx,1

           !--------------
           ! z-directional limiter, lower face
           ! compute gradients of E0, E
           E0avg  = (E0(i,j,k) + E0(i,j,k-1))*0.5d0
           E0d_zl = E0(i,j,k) - E0(i,j,k-1)
           Ed_zl  = E(i,j,k)  - E(i,j,k-1)

           !    compute R for limiter 
           R  = max(dzi *abs(E0d_zl)/E0avg, Rmin)
           R0 = max(dzi0*abs(E0d_zl)/E0avg, Rmin)

           !    compute average opacity over face
           kap = (sHI*(  HI(i,j,k)   + HI(i,j,k-1))  &
                + sHeI*( HeI(i,j,k)  + HeI(i,j,k-1)) &
                + sHeII*(HeII(i,j,k) + HeII(i,j,k-1)))*0.5d0*nUn
           kap0 = (sHI*( HI0(i,j,k)    + HI0(i,j,k-1))  &
                 + sHeI*( HeI0(i,j,k)  + HeI0(i,j,k-1)) &
                 + sHeII*(HeII0(i,j,k) + HeII0(i,j,k-1)))*0.5d0*nUn0

           !    compute limiter
           if (LType == 1) then       ! rational approx. to LP lim. (LP, 1981)
              D_zl  = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
              D0_zl = c*(2.d0*kap0+R0)/(6.d0*kap0*kap0+3.d0*kap0*R0+R0*R0)
           else if (LType == 2) then  ! Larsen n=2 limiter
              D_zl  = c/sqrt(kap*kap*9.d0 + R*R)
              D0_zl = c/sqrt(kap0*kap0*9.d0 + R0*R0)
           else if (LType == 3) then  ! no limiter
              D_zl  = c/kap/3.d0
              D0_zl = c/kap0/3.d0
           else if (LType == 4) then  ! Zeus limiter
              D_zl  = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
              D0_zl = c*(2.d0*kap0+R0)/(6.d0*kap0*kap0+3.d0*kap0*R0+R0*R0)
           else                       ! standard Levermore-Pomraning (LP, 1981)
              D_zl  = c*(cosh(R/kap)/sinh(R/kap)-kap/R)/R
              D0_zl = c*(cosh(R0/kap0)/sinh(R0/kap0)-kap0/R0)/R0
           endif


           !--------------
           ! y-directional limiter, lower face
           ! compute gradients of E0, E
           E0avg  = (E0(i,j,k) + E0(i,j-1,k))*0.5d0
           E0d_yl = E0(i,j,k) - E0(i,j-1,k)
           Ed_yl  = E(i,j,k)  - E(i,j-1,k)

           !    compute R for limiter 
           R  = max(dyi *abs(E0d_yl)/E0avg, Rmin)
           R0 = max(dyi0*abs(E0d_yl)/E0avg, Rmin)

           !    compute average opacity over face
           kap = (sHI*(  HI(i,j,k)   + HI(i,j-1,k))  &
                + sHeI*( HeI(i,j,k)  + HeI(i,j-1,k)) &
                + sHeII*(HeII(i,j,k) + HeII(i,j-1,k)))*0.5d0*nUn
           kap0 = (sHI*( HI0(i,j,k)    + HI0(i,j-1,k))  &
                 + sHeI*( HeI0(i,j,k)  + HeI0(i,j-1,k)) &
                 + sHeII*(HeII0(i,j,k) + HeII0(i,j-1,k)))*0.5d0*nUn0

           !    compute limiter
           if (LType == 1) then       ! rational approx. to LP lim. (LP, 1981)
              D_yl  = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
              D0_yl = c*(2.d0*kap0+R0)/(6.d0*kap0*kap0+3.d0*kap0*R0+R0*R0)
           else if (LType == 2) then  ! Larsen n=2 limiter
              D_yl  = c/sqrt(kap*kap*9.d0 + R*R)
              D0_yl = c/sqrt(kap0*kap0*9.d0 + R0*R0)
           else if (LType == 3) then  ! no limiter
              D_yl  = c/kap/3.d0
              D0_yl = c/kap0/3.d0
           else if (LType == 4) then  ! Zeus limiter
              D_yl  = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
              D0_yl = c*(2.d0*kap0+R0)/(6.d0*kap0*kap0+3.d0*kap0*R0+R0*R0)
           else                       ! standard Levermore-Pomraning (LP, 1981)
              D_yl  = c*(cosh(R/kap)/sinh(R/kap)-kap/R)/R
              D0_yl = c*(cosh(R0/kap0)/sinh(R0/kap0)-kap0/R0)/R0
           endif

           !--------------
           ! x-directional limiter, lower face
           ! compute gradients of E0, E
           E0avg  = (E0(i,j,k) + E0(i-1,j,k))*0.5d0
           E0d_xl = E0(i,j,k) - E0(i-1,j,k)
           Ed_xl  = E(i,j,k)  - E(i-1,j,k)

           !    compute R for limiter 
           R  = max(dxi *abs(E0d_xl)/E0avg, Rmin)
           R0 = max(dxi0*abs(E0d_xl)/E0avg, Rmin)

           !    compute opacity
           kap = (sHI*(  HI(i,j,k)   + HI(i-1,j,k))  &
                + sHeI*( HeI(i,j,k)  + HeI(i-1,j,k)) &
                + sHeII*(HeII(i,j,k) + HeII(i-1,j,k)))*0.5d0*nUn
           kap0 = (sHI*( HI0(i,j,k)    + HI0(i-1,j,k))  &
                 + sHeI*( HeI0(i,j,k)  + HeI0(i-1,j,k)) &
                 + sHeII*(HeII0(i,j,k) + HeII0(i-1,j,k)))*0.5d0*nUn0

           !    compute limiter
           if (LType == 1) then       ! rational approx. to LP lim. (LP, 1981)
              D_xl  = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
              D0_xl = c*(2.d0*kap0+R0)/(6.d0*kap0*kap0+3.d0*kap0*R0+R0*R0)
           else if (LType == 2) then  ! Larsen n=2 limiter
              D_xl  = c/sqrt(kap*kap*9.d0 + R*R)
              D0_xl = c/sqrt(kap0*kap0*9.d0 + R0*R0)
           else if (LType == 3) then  ! no limiter
              D_xl  = c/kap/3.d0
              D0_xl = c/kap0/3.d0
           else if (LType == 4) then  ! Zeus limiter
              D_xl  = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
              D0_xl = c*(2.d0*kap0+R0)/(6.d0*kap0*kap0+3.d0*kap0*R0+R0*R0)
           else                       ! standard Levermore-Pomraning (LP, 1981)
              D_xl  = c*(cosh(R/kap)/sinh(R/kap)-kap/R)/R
              D0_xl = c*(cosh(R0/kap0)/sinh(R0/kap0)-kap0/R0)/R0
           endif

           !--------------
           ! x-directional limiter, upper face
           ! compute gradients of E0, E
           E0avg  = (E0(i+1,j,k) + E0(i,j,k))*0.5d0
           E0d_xr = E0(i+1,j,k) - E0(i,j,k)
           Ed_xr  = E(i+1,j,k)  - E(i,j,k)

           !    compute R for limiter 
           R  = max(dxi *abs(E0d_xr)/E0avg, Rmin)
           R0 = max(dxi0*abs(E0d_xr)/E0avg, Rmin)

           !    compute opacity
           kap = (sHI*(  HI(i+1,j,k)   + HI(i,j,k))  &
                + sHeI*( HeI(i+1,j,k)  + HeI(i,j,k)) &
                + sHeII*(HeII(i+1,j,k) + HeII(i,j,k)))*0.5d0*nUn
           kap0 = (sHI*( HI0(i+1,j,k)    + HI0(i,j,k))  &
                 + sHeI*( HeI0(i+1,j,k)  + HeI0(i,j,k)) &
                 + sHeII*(HeII0(i+1,j,k) + HeII0(i,j,k)))*0.5d0*nUn0

           !    compute limiter
           if (LType == 1) then       ! rational approx. to LP lim. (LP, 1981)
              D_xr  = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
              D0_xr = c*(2.d0*kap0+R0)/(6.d0*kap0*kap0+3.d0*kap0*R0+R0*R0)
           else if (LType == 2) then  ! Larsen n=2 limiter
              D_xr  = c/sqrt(kap*kap*9.d0 + R*R)
              D0_xr = c/sqrt(kap0*kap0*9.d0 + R0*R0)
           else if (LType == 3) then  ! no limiter
              D_xr  = c/kap/3.d0
              D0_xr = c/kap0/3.d0
           else if (LType == 4) then  ! Zeus limiter
              D_xr  = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
              D0_xr = c*(2.d0*kap0+R0)/(6.d0*kap0*kap0+3.d0*kap0*R0+R0*R0)
           else                       ! standard Levermore-Pomraning (LP, 1981)
              D_xr  = c*(cosh(R/kap)/sinh(R/kap)-kap/R)/R
              D0_xr = c*(cosh(R0/kap0)/sinh(R0/kap0)-kap0/R0)/R0
           endif

           !--------------
           ! y-directional limiter, upper face
           ! compute gradients of E0, E
           E0avg  = (E0(i,j+1,k) + E0(i,j,k))*0.5d0
           E0d_yr = E0(i,j+1,k) - E0(i,j,k)
           Ed_yr  = E(i,j+1,k)  - E(i,j,k)

           !    compute R for limiter 
           R  = max(dyi *abs(E0d_yr)/E0avg, Rmin)
           R0 = max(dyi0*abs(E0d_yr)/E0avg, Rmin)

           !    compute opacity
           kap = (sHI*(  HI(i,j+1,k)   + HI(i,j,k))  &
                + sHeI*( HeI(i,j+1,k)  + HeI(i,j,k)) &
                + sHeII*(HeII(i,j+1,k) + HeII(i,j,k)))*0.5d0*nUn
           kap0 = (sHI*( HI0(i,j+1,k)    + HI0(i,j,k))  &
                 + sHeI*( HeI0(i,j+1,k)  + HeI0(i,j,k)) &
                 + sHeII*(HeII0(i,j+1,k) + HeII0(i,j,k)))*0.5d0*nUn0

           !    compute limiter
           if (LType == 1) then       ! rational approx. to LP lim. (LP, 1981)
              D_yr  = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
              D0_yr = c*(2.d0*kap0+R0)/(6.d0*kap0*kap0+3.d0*kap0*R0+R0*R0)
           else if (LType == 2) then  ! Larsen n=2 limiter
              D_yr  = c/sqrt(kap*kap*9.d0 + R*R)
              D0_yr = c/sqrt(kap0*kap0*9.d0 + R0*R0)
           else if (LType == 3) then  ! no limiter
              D_yr  = c/kap/3.d0
              D0_yr = c/kap0/3.d0
           else if (LType == 4) then  ! Zeus limiter
              D_yr  = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
              D0_yr = c*(2.d0*kap0+R0)/(6.d0*kap0*kap0+3.d0*kap0*R0+R0*R0)
           else                       ! standard Levermore-Pomraning (LP, 1981)
              D_yr  = c*(cosh(R/kap)/sinh(R/kap)-kap/R)/R
              D0_yr = c*(cosh(R0/kap0)/sinh(R0/kap0)-kap0/R0)/R0
           endif

           !--------------
           ! z-directional limiter, upper face
           ! compute gradients of E0, E
           E0avg  = (E0(i,j,k+1) + E0(i,j,k))*0.5d0
           E0d_zr = E0(i,j,k+1) - E0(i,j,k)
           Ed_zr  = E(i,j,k+1)  - E(i,j,k)

           !    compute R for limiter 
           R  = max(dzi *abs(E0d_zr)/E0avg, Rmin)
           R0 = max(dzi0*abs(E0d_zr)/E0avg, Rmin)

           !    compute opacity
           kap = (sHI*(  HI(i,j,k+1)   + HI(i,j,k))  &
                + sHeI*( HeI(i,j,k+1)  + HeI(i,j,k)) &
                + sHeII*(HeII(i,j,k+1) + HeII(i,j,k)))*0.5d0*nUn
           kap0 = (sHI*( HI0(i,j,k+1)    + HI0(i,j,k))  &
                 + sHeI*( HeI0(i,j,k+1)  + HeI0(i,j,k)) &
                 + sHeII*(HeII0(i,j,k+1) + HeII0(i,j,k)))*0.5d0*nUn0

           !    compute limiter
           if (LType == 1) then       ! rational approx. to LP lim. (LP, 1981)
              D_zr  = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
              D0_zr = c*(2.d0*kap0+R0)/(6.d0*kap0*kap0+3.d0*kap0*R0+R0*R0)
           else if (LType == 2) then  ! Larsen n=2 limiter
              D_zr  = c/sqrt(kap*kap*9.d0 + R*R)
              D0_zr = c/sqrt(kap0*kap0*9.d0 + R0*R0)
           else if (LType == 3) then  ! no limiter
              D_zr  = c/kap/3.d0
              D0_zr = c/kap0/3.d0
           else if (LType == 4) then  ! Zeus limiter
              D_zr  = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
              D0_zr = c*(2.d0*kap0+R0)/(6.d0*kap0*kap0+3.d0*kap0*R0+R0*R0)
           else                       ! standard Levermore-Pomraning (LP, 1981)
              D_zr  = c*(cosh(R/kap)/sinh(R/kap)-kap/R)/R
              D0_zr = c*(cosh(R0/kap0)/sinh(R0/kap0)-kap0/R0)/R0
           endif

           ! opacity values in this cell
           kap = (sHI*HI(i,j,k) + sHeI*HeI(i,j,k) + sHeII*HeII(i,j,k))*nUn
           kap0 = (sHI*HI0(i,j,k) + sHeI*HeI0(i,j,k) + sHeII*HeII0(i,j,k))*nUn0

           kap_max = max(kap_max, kap)
           kap_avg = kap_avg + kap
           nHI_max = max(nHI_max, HI(i,j,k))
           nHI_avg = nHI_avg + HI(i,j,k)

           ! set the matrix entries
           mat(:,i,j,k) = (/  &
                -dtfac*dzi*dzi*D_zl, &       ! z-left
                -dtfac*dyi*dyi*D_yl, &       ! y-left
                -dtfac*dxi*dxi*D_xl, &       ! x-left
                1.d0 + dtfac*(c*kap + dxi*dxi*(D_xl+D_xr)          &  ! self
                     + dyi*dyi*(D_yl+D_yr) + dzi*dzi*(D_zl+D_zr)), &
                -dtfac*dxi*dxi*D_xr, &       ! x-right
                -dtfac*dyi*dyi*D_yr, &       ! y-right
                -dtfac*dzi*dzi*D_zr /)       ! z-right
                
           ! set the rhs entries
           rhs(i,j,k) = ( (dtfac + dtfac0)*src(i,j,k)                   &
                        + (1.d0 - dtfac0*c*kap0)*E0(i,j,k)              &
                        + dtfac0*dzi0*dzi0*(D0_zr*E0d_zr-D0_zl*E0d_zl)  &
                        + dtfac0*dyi0*dyi0*(D0_yr*E0d_yr-D0_yl*E0d_yl)  &
                        + dtfac0*dxi0*dxi0*(D0_xr*E0d_xr-D0_xl*E0d_xl)  &
                        - (1.d0 + dtfac*c*kap)*E(i,j,k)                 &
                        + dtfac*dzi*dzi*(D_zr*Ed_zr-D_zl*Ed_zl)         &
                        + dtfac*dyi*dyi*(D_yr*Ed_yr-D_yl*Ed_yl)         &
                        + dtfac*dxi*dxi*(D_xr*Ed_xr-D_xl*Ed_xl) )
           
        enddo
     enddo
  enddo

  print *, 'MFSplit_SetupSystem:'
  print *, '  sHI =',sHI,',  nUn =',nUn
  print *, '  max(kap) =',kap_max,',  avg(kap) =',kap_avg/Nx/Ny/Nz
  print *, '  max(HI)  =',nHI_max,',  avg(HI)  =',nHI_avg/Nx/Ny/Nz

!!$  print *,'  sum(mat) =',sum(mat)
!!$  print *,'  sum(rhs) =',sum(rhs)



  ! update matrix/rhs based on boundary conditions/location
  !    z-left face
  if (zlface == 1) then
     ! Dirichlet
     if (BCZl==1) then
        k = x2s
        do j=1,Ny
           do i=1,Nx
              mat(1,i,j,k) = 0.d0
           enddo
        enddo
     ! Neumann
     else if (BCZl==2) then
        k = x2s
        do j=1,Ny
           do i=1,Nx
              mat(4,i,j,k) = mat(4,i,j,k) + mat(1,i,j,k)
              mat(1,i,j,k) = 0.d0
           enddo
        enddo
     endif
  end if

  !    y-left face
  if (ylface == 1) then
     ! Dirichlet
     if (BCYl==1) then
        j = x1s
        do k=1,Nz
           do i=1,Nx
              mat(2,i,j,k) = 0.d0
           enddo
        enddo
     ! Neumann
     else if (BCYl==2) then
        j = x1s
        do k=1,Nz
           do i=1,Nx
              mat(4,i,j,k) = mat(4,i,j,k) + mat(2,i,j,k)
              mat(2,i,j,k) = 0.d0
           enddo
        enddo
     endif
  end if

  !    x-left face
  if (xlface == 1) then
     ! Dirichlet
     if (BCXl==1) then
        i = x0s
        do k=1,Nz
           do j=1,Ny
              mat(3,i,j,k) = 0.d0
           enddo
        enddo
     ! Neumann
     else if (BCXl==2) then
        i = x0s
        do k=1,Nz
           do j=1,Ny
              mat(4,i,j,k) = mat(4,i,j,k) + mat(3,i,j,k)
              mat(3,i,j,k) = 0.d0
           enddo
        enddo
     endif
  end if

  !    x-right face
  if (xrface==1) then
     ! Dirichlet
     if (BCXr==1) then
        i = x0e
        do k=1,Nz
           do j=1,Ny
              mat(5,i,j,k) = 0.d0
           enddo
        enddo
     ! Neumann
     else if (BCXr==2) then
        i = x0e
        do k=1,Nz
           do j=1,Ny
              mat(4,i,j,k) = mat(4,i,j,k) + mat(5,i,j,k)
              mat(5,i,j,k) = 0.d0
           enddo
        enddo
     endif
  endif

  !    y-right face
  if (yrface==1) then
     ! Dirichlet
     if (BCYr==1) then
        j = x1e
        do k=1,Nz
           do i=1,Nx
              mat(6,i,j,k) = 0.d0
           enddo
        enddo
     ! Neumann
     else if (BCYr==2) then
        j = x1e
        do k=1,Nz
           do i=1,Nx
              mat(4,i,j,k) = mat(4,i,j,k) + mat(6,i,j,k)
              mat(6,i,j,k) = 0.d0
           enddo
        enddo
     endif
  endif

  !    z-right face
  if (zrface==1) then
     ! Dirichlet
     if (BCZr==1) then
        k = x2e
        do j=1,Ny
           do i=1,Nx
              mat(7,i,j,k) = 0.d0
           enddo
        enddo
     ! Neumann
     else if (BCZr==2) then
        k = x2e
        do j=1,Ny
           do i=1,Nx
              mat(4,i,j,k) = mat(4,i,j,k) + mat(7,i,j,k)
              mat(7,i,j,k) = 0.d0
           enddo
        enddo
     endif
  endif

  rhsnorm = sum(rhs*rhs)

!!$  print *,'MFSplit_SetupSystem_3D:'
!!$  print *,'  sum(mat) =',sum(mat)
!!$  print *,'  sum(rhs) =',sum(rhs)
!!$  print *,'   rhsnorm =',rhsnorm

  return
end subroutine MFSplit_SetupSystem_3D
!=======================================================================






subroutine MFSplit_SetupSystem_2D(mat, rhs, rhsnorm, E, E0, HI, HI0, HeI, &
     HeI0, HeII, HeII0, src, LType, dt, theta, sHI, sHeI, sHeII, a, a0,   &
     aUn, lUn, lUn0, nUn, nUn0, dx, dy, BCXl, BCXr, BCYl, BCYr, x0s, x0e, &
     x1s, x1e, Nx, Ny, NGxl, NGxr, NGyl, NGyr, xlface, xrface, ylface,    &
     yrface, ier)
  !=======================================================================
  !  PURPOSE: 2D version of the routine
  !=======================================================================
#include "fortran.def"
  implicit none
  
  !--------------
  ! argument declarations
  integer, intent(in)  :: LType
  integer,  intent(in) :: BCXl, BCXr, x0s, x0e, Nx, NGxl, NGxr, xlface, xrface
  integer,  intent(in) :: BCYl, BCYr, x1s, x1e, Ny, NGyl, NGyr, ylface, yrface
  integer, intent(out) :: ier
  REALSUB, intent(in)  :: a, a0
  REAL, intent(in) :: dx, dy, dt, theta, sHI, sHeI, sHeII
  REAL, intent(in) :: aUn, lUn, lUn0, nUn, nUn0
  REAL, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr), &
       intent(in), target :: E, E0, HI, HI0, HeI, HeI0, HeII, HeII0, src
  real*8, intent(out) :: mat(5,x0s:x0e,x1s:x1e)
  real*8, intent(out) :: rhs(x0s:x0e,x1s:x1e)
  REAL, intent(out) :: rhsnorm

  !--------------
  ! locals
  integer :: i, j
  real*8 :: dtfac, dtfac0, c, pi
  real*8 :: dxi, dxi0, dyi, dyi0
  real*8 :: kap, kap0, E0avg, R, R0, Rmin
  real*8 :: D_xl, D0_xl, D_xr, D0_xr, E0d_xl, E0d_xr, Ed_xl, Ed_xr
  real*8 :: D_yl, D0_yl, D_yr, D0_yr, E0d_yl, E0d_yr, Ed_yl, Ed_yr


  !=======================================================================
  
  ! initialize outputs to zero, flag to success
  mat = 0.d0
  rhs = 0.d0
  ier = 1

  ! set shortcut values
  dtfac  = dt*theta
  dtfac0 = dt*(1.d0-theta)
  dxi    = a/dx/lUn
  dyi    = a/dy/lUn
  dxi0   = a0/dx/lUn0
  dyi0   = a0/dy/lUn0
  c      = 2.99792458d10     ! speed of light [cm/s]
  pi     = 4.d0*datan(1.d0)
  Rmin   = 1.0d-20

!!$  print *, 'entering MFSplit_SetupSystem_3D:'
!!$  print *, '    LType =',LType
!!$  print *, '    BCx* =',BCXl,BCXr
!!$  print *, '    BCy* =',BCYl,BCYr
!!$  print *, '    x0* =',x0s,x0e
!!$  print *, '    x1* =',x1s,x1e
!!$  print *, '    x2* =',x2s,x2e
!!$  print *, '    Nx* =',Nx,NGxl,NGxr
!!$  print *, '    Ny* =',Ny,NGyl,NGyr
!!$  print *, '    x*face =',xlface,xrface
!!$  print *, '    y*face =',ylface,yrface
!!$  print *, '    a =',a
!!$  print *, '    dx* =',dx,dy
!!$  print *, '    dt =',dt
!!$  print *, '    theta =',theta
!!$  print *, '    sH* =',sHI,sHeI,sHeII
!!$  print *, '    units =',lUn,nUn,nUn0
!!$  print *, '    E,E0 =',sum(E),sum(E0)
!!$  print *, '    HI,HI0 =',sum(HI),sum(HI0)
!!$  print *, '    HeI,HeI0 =',sum(HeI),sum(HeI0)
!!$  print *, '    HeII,HeII0 =',sum(HeII),sum(HeII0)


  ! iterate over the active domain
  do j=1,Ny,1
     do i=1,Nx,1
        
        !--------------
        ! y-directional limiter, lower face
        ! compute gradients of E0, E
        E0avg  = (E0(i,j) + E0(i,j-1))*0.5d0
        E0d_yl = E0(i,j) - E0(i,j-1)
        Ed_yl  = E(i,j)  - E(i,j-1)
        
        !    compute R for limiter 
        R  = max(dyi *abs(E0d_yl)/E0avg, Rmin)
        R0 = max(dyi0*abs(E0d_yl)/E0avg, Rmin)
        
        !    compute average opacity over face
        kap = (sHI*(  HI(i,j)   + HI(i,j-1))  &
             + sHeI*( HeI(i,j)  + HeI(i,j-1)) &
             + sHeII*(HeII(i,j) + HeII(i,j-1)))*0.5d0*nUn
        kap0 = (sHI*( HI0(i,j)    + HI0(i,j-1))  &
              + sHeI*( HeI0(i,j)  + HeI0(i,j-1)) &
              + sHeII*(HeII0(i,j) + HeII0(i,j-1)))*0.5d0*nUn0

        !    compute limiter
        if (LType == 1) then       ! rational approx. to LP lim. (LP, 1981)
           D_yl  = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
           D0_yl = c*(2.d0*kap0+R0)/(6.d0*kap0*kap0+3.d0*kap0*R0+R0*R0)
        else if (LType == 2) then  ! Larsen n=2 limiter
           D_yl  = c/sqrt(kap*kap*9.d0 + R*R)
           D0_yl = c/sqrt(kap0*kap0*9.d0 + R0*R0)
        else if (LType == 3) then  ! no limiter
           D_yl  = c/kap/3.d0
           D0_yl = c/kap0/3.d0
        else if (LType == 4) then  ! Zeus limiter
           D_yl  = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
           D0_yl = c*(2.d0*kap0+R0)/(6.d0*kap0*kap0+3.d0*kap0*R0+R0*R0)
        else                       ! standard Levermore-Pomraning (LP, 1981)
           D_yl  = c*(cosh(R/kap)/sinh(R/kap)-kap/R)/R
           D0_yl = c*(cosh(R0/kap0)/sinh(R0/kap0)-kap0/R0)/R0
        endif
        
        !--------------
        ! x-directional limiter, lower face
        ! compute gradients of E0, E
        E0avg  = (E0(i,j) + E0(i-1,j))*0.5d0
        E0d_xl = E0(i,j) - E0(i-1,j)
        Ed_xl  = E(i,j)  - E(i-1,j)
        
        !    compute R for limiter 
        R  = max(dxi *abs(E0d_xl)/E0avg, Rmin)
        R0 = max(dxi0*abs(E0d_xl)/E0avg, Rmin)
        
        !    compute opacity
        kap = (sHI*(  HI(i,j)   + HI(i-1,j))  &
             + sHeI*( HeI(i,j)  + HeI(i-1,j)) &
             + sHeII*(HeII(i,j) + HeII(i-1,j)))*0.5d0*nUn
        kap0 = (sHI*( HI0(i,j)    + HI0(i-1,j))  &
              + sHeI*( HeI0(i,j)  + HeI0(i-1,j)) &
              + sHeII*(HeII0(i,j) + HeII0(i-1,j)))*0.5d0*nUn0
        
        !    compute limiter
        if (LType == 1) then       ! rational approx. to LP lim. (LP, 1981)
           D_xl  = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
           D0_xl = c*(2.d0*kap0+R0)/(6.d0*kap0*kap0+3.d0*kap0*R0+R0*R0)
        else if (LType == 2) then  ! Larsen n=2 limiter
           D_xl  = c/sqrt(kap*kap*9.d0 + R*R)
           D0_xl = c/sqrt(kap0*kap0*9.d0 + R0*R0)
        else if (LType == 3) then  ! no limiter
           D_xl  = c/kap/3.d0
           D0_xl = c/kap0/3.d0
        else if (LType == 4) then  ! Zeus limiter
           D_xl  = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
           D0_xl = c*(2.d0*kap0+R0)/(6.d0*kap0*kap0+3.d0*kap0*R0+R0*R0)
        else                       ! standard Levermore-Pomraning (LP, 1981)
           D_xl  = c*(cosh(R/kap)/sinh(R/kap)-kap/R)/R
           D0_xl = c*(cosh(R0/kap0)/sinh(R0/kap0)-kap0/R0)/R0
        endif
        
        !--------------
        ! x-directional limiter, upper face
        ! compute gradients of E0, E
        E0avg  = (E0(i+1,j) + E0(i,j))*0.5d0
        E0d_xr = E0(i+1,j) - E0(i,j)
        Ed_xr  = E(i+1,j)  - E(i,j)
        
        !    compute R for limiter 
        R  = max(dxi *abs(E0d_xr)/E0avg, Rmin)
        R0 = max(dxi0*abs(E0d_xr)/E0avg, Rmin)
        
        !    compute opacity
        kap = (sHI*(  HI(i+1,j)   + HI(i,j))  &
             + sHeI*( HeI(i+1,j)  + HeI(i,j)) &
             + sHeII*(HeII(i+1,j) + HeII(i,j)))*0.5d0*nUn
        kap0 = (sHI*( HI0(i+1,j)    + HI0(i,j))  &
              + sHeI*( HeI0(i+1,j)  + HeI0(i,j)) &
              + sHeII*(HeII0(i+1,j) + HeII0(i,j)))*0.5d0*nUn0

        !    compute limiter
        if (LType == 1) then       ! rational approx. to LP lim. (LP, 1981)
           D_xr  = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
           D0_xr = c*(2.d0*kap0+R0)/(6.d0*kap0*kap0+3.d0*kap0*R0+R0*R0)
        else if (LType == 2) then  ! Larsen n=2 limiter
           D_xr  = c/sqrt(kap*kap*9.d0 + R*R)
           D0_xr = c/sqrt(kap0*kap0*9.d0 + R0*R0)
        else if (LType == 3) then  ! no limiter
           D_xr  = c/kap/3.d0
           D0_xr = c/kap0/3.d0
        else if (LType == 4) then  ! Zeus limiter
           D_xr  = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
           D0_xr = c*(2.d0*kap0+R0)/(6.d0*kap0*kap0+3.d0*kap0*R0+R0*R0)
        else                       ! standard Levermore-Pomraning (LP, 1981)
           D_xr  = c*(cosh(R/kap)/sinh(R/kap)-kap/R)/R
           D0_xr = c*(cosh(R0/kap0)/sinh(R0/kap0)-kap0/R0)/R0
        endif
        
        !--------------
        ! y-directional limiter, upper face
        ! compute gradients of E0, E
        E0avg  = (E0(i,j+1) + E0(i,j))*0.5d0
        E0d_yr = E0(i,j+1) - E0(i,j)
        Ed_yr  = E(i,j+1)  - E(i,j)
        
        !    compute R for limiter 
        R  = max(dyi *abs(E0d_yr)/E0avg, Rmin)
        R0 = max(dyi0*abs(E0d_yr)/E0avg, Rmin)
        
        !    compute opacity
        kap = (sHI*(  HI(i,j+1)   + HI(i,j))  &
             + sHeI*( HeI(i,j+1)  + HeI(i,j)) &
             + sHeII*(HeII(i,j+1) + HeII(i,j)))*0.5d0*nUn
        kap0 = (sHI*( HI0(i,j+1)    + HI0(i,j))  &
              + sHeI*( HeI0(i,j+1)  + HeI0(i,j)) &
              + sHeII*(HeII0(i,j+1) + HeII0(i,j)))*0.5d0*nUn0

        !    compute limiter
        if (LType == 1) then       ! rational approx. to LP lim. (LP, 1981)
           D_yr  = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
           D0_yr = c*(2.d0*kap0+R0)/(6.d0*kap0*kap0+3.d0*kap0*R0+R0*R0)
        else if (LType == 2) then  ! Larsen n=2 limiter
           D_yr  = c/sqrt(kap*kap*9.d0 + R*R)
           D0_yr = c/sqrt(kap0*kap0*9.d0 + R0*R0)
        else if (LType == 3) then  ! no limiter
           D_yr  = c/kap/3.d0
           D0_yr = c/kap0/3.d0
        else if (LType == 4) then  ! Zeus limiter
           D_yr  = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
           D0_yr = c*(2.d0*kap0+R0)/(6.d0*kap0*kap0+3.d0*kap0*R0+R0*R0)
        else                       ! standard Levermore-Pomraning (LP, 1981)
           D_yr  = c*(cosh(R/kap)/sinh(R/kap)-kap/R)/R
           D0_yr = c*(cosh(R0/kap0)/sinh(R0/kap0)-kap0/R0)/R0
        endif
        
        ! opacity values in this cell
        kap = (sHI*HI(i,j) + sHeI*HeI(i,j) + sHeII*HeII(i,j))*nUn
        kap0 = (sHI*HI0(i,j) + sHeI*HeI0(i,j) + sHeII*HeII0(i,j))*nUn0
        
        ! set the matrix entries
        mat(:,i,j) = (/  &
             -dtfac*dyi*dyi*D_yl, &       ! y-left
             -dtfac*dxi*dxi*D_xl, &       ! x-left
             1.d0 + dtfac*(c*kap + dxi*dxi*(D_xl+D_xr)  &  ! self
                  + dyi*dyi*(D_yl+D_yr)), &
             -dtfac*dxi*dxi*D_xr, &       ! x-right
             -dtfac*dyi*dyi*D_yr /)       ! y-right
                
        ! set the rhs entries
        rhs(i,j) = ( (dtfac + dtfac0)*src(i,j)                       &
                     + (1.d0 - dtfac0*c*kap0)*E0(i,j)                &
                     + dtfac0*dyi0*dyi0*(D0_yr*E0d_yr-D0_yl*E0d_yl)  &
                     + dtfac0*dxi0*dxi0*(D0_xr*E0d_xr-D0_xl*E0d_xl)  &
                     - (1.d0 + dtfac*c*kap)*E(i,j)                   &
                     + dtfac*dyi*dyi*(D_yr*Ed_yr-D_yl*Ed_yl)         &
                     + dtfac*dxi*dxi*(D_xr*Ed_xr-D_xl*Ed_xl) )
           
     enddo
  enddo

!!$  print *,'  sum(mat) =',sum(mat)
!!$  print *,'  sum(rhs) =',sum(rhs)



  ! update matrix/rhs based on boundary conditions/location
  !    y-left face
  if (ylface == 1) then
     ! Dirichlet
     if (BCYl==1) then
        j = x1s
        do i=1,Nx
           mat(1,i,j) = 0.d0
        enddo
     ! Neumann
     else if (BCYl==2) then
        j = x1s
        do i=1,Nx
           mat(3,i,j) = mat(3,i,j) + mat(1,i,j)
           mat(1,i,j) = 0.d0
        enddo
     endif
  end if
  
  !    x-left face
  if (xlface == 1) then
     ! Dirichlet
     if (BCXl==1) then
        i = x0s
        do j=1,Ny
           mat(2,i,j) = 0.d0
        enddo
     ! Neumann
     else if (BCXl==2) then
        i = x0s
        do j=1,Ny
           mat(3,i,j) = mat(3,i,j) + mat(2,i,j)
           mat(2,i,j) = 0.d0
        enddo
     endif
  end if
  
  !    x-right face
  if (xrface==1) then
     ! Dirichlet
     if (BCXr==1) then
        i = x0e
        do j=1,Ny
           mat(4,i,j) = 0.d0
        enddo
     ! Neumann
     else if (BCXr==2) then
        i = x0e
        do j=1,Ny
           mat(3,i,j) = mat(3,i,j) + mat(4,i,j)
           mat(4,i,j) = 0.d0
        enddo
     endif
  endif

  !    y-right face
  if (yrface==1) then
     ! Dirichlet
     if (BCYr==1) then
        j = x1e
        do i=1,Nx
           mat(5,i,j) = 0.d0
        enddo
     ! Neumann
     else if (BCYr==2) then
        j = x1e
        do i=1,Nx
           mat(3,i,j) = mat(3,i,j) + mat(5,i,j)
           mat(5,i,j) = 0.d0
        enddo
     endif
  endif

  rhsnorm = sum(rhs*rhs)

!!$  print *,'  sum(mat) =',sum(mat)
!!$  print *,'  sum(rhs) =',sum(rhs)
  
  return
end subroutine MFSplit_SetupSystem_2D
!=======================================================================






subroutine MFSplit_SetupSystem_1D(mat, rhs, rhsnorm, E, E0, HI, HI0, HeI, &
     HeI0, HeII, HeII0, src, LType, dt, theta, sHI, sHeI, sHeII, a, a0,   &
     aUn, lUn, lUn0, nUn, nUn0, dx, BCXl, BCXr, x0s, x0e, Nx, NGxl, NGxr, &
     xlface, xrface, ier)
  !=======================================================================
  !  PURPOSE: 1D version of the routine
  !=======================================================================
#include "fortran.def"
  implicit none
  
  !--------------
  ! argument declarations
  integer, intent(in)  :: LType
  integer,  intent(in) :: BCXl, BCXr, x0s, x0e, Nx, NGxl, NGxr, xlface, xrface
  integer, intent(out) :: ier
  REALSUB, intent(in)  :: a, a0
  REAL, intent(in) :: dx, dt, theta, sHI, sHeI, sHeII
  REAL, intent(in) :: aUn, lUn, lUn0, nUn, nUn0
  REAL, dimension(1-NGxl:Nx+NGxr), intent(in), target :: E, E0, HI, HI0, &
       HeI, HeI0, HeII, HeII0, src
  real*8, intent(out) :: mat(3,x0s:x0e)
  real*8, intent(out) :: rhs(x0s:x0e)
  REAL, intent(out) :: rhsnorm

  !--------------
  ! locals
  integer :: i
  real*8 :: dtfac, dtfac0, c, pi
  real*8 :: dxi, dxi0
  real*8 :: kap, kap0, E0avg, R, R0, Rmin
  real*8 :: D_xl, D0_xl, D_xr, D0_xr, E0d_xl, E0d_xr, Ed_xl, Ed_xr


  !=======================================================================
  
  ! initialize outputs to zero, flag to success
  mat = 0.d0
  rhs = 0.d0
  ier = 1

  ! set shortcut values
  dtfac  = dt*theta
  dtfac0 = dt*(1.d0-theta)
  dxi    = a/dx/lUn
  dxi0   = a0/dx/lUn0
  c      = 2.99792458d10     ! speed of light [cm/s]
  pi     = 4.d0*datan(1.d0)
  Rmin   = 1.0d-20

!!$  print *, 'entering MFSplit_SetupSystem_3D:'
!!$  print *, '    LType =',LType
!!$  print *, '    BCx* =',BCXl,BCXr
!!$  print *, '    x0* =',x0s,x0e
!!$  print *, '    x1* =',x1s,x1e
!!$  print *, '    x2* =',x2s,x2e
!!$  print *, '    Nx* =',Nx,NGxl,NGxr
!!$  print *, '    x*face =',xlface,xrface
!!$  print *, '    a =',a
!!$  print *, '    dx =',dx
!!$  print *, '    dt =',dt
!!$  print *, '    theta =',theta
!!$  print *, '    sH* =',sHI,sHeI,sHeII
!!$  print *, '    units =',lUn,nUn,nUn0
!!$  print *, '    E,E0 =',sum(E),sum(E0)
!!$  print *, '    HI,HI0 =',sum(HI),sum(HI0)
!!$  print *, '    HeI,HeI0 =',sum(HeI),sum(HeI0)
!!$  print *, '    HeII,HeII0 =',sum(HeII),sum(HeII0)


  ! iterate over the active domain
  do i=1,Nx,1
        
     !--------------
     ! x-directional limiter, lower face
     ! compute gradients of E0, E
     E0avg  = (E0(i) + E0(i-1))*0.5d0
     E0d_xl = E0(i) - E0(i-1)
     Ed_xl  = E(i)  - E(i-1)
     
     !    compute R for limiter 
     R  = max(dxi *abs(E0d_xl)/E0avg, Rmin)
     R0 = max(dxi0*abs(E0d_xl)/E0avg, Rmin)
     
     !    compute opacity
     kap = (sHI*(  HI(i)   + HI(i-1))  &
          + sHeI*( HeI(i)  + HeI(i-1)) &
          + sHeII*(HeII(i) + HeII(i-1)))*0.5d0*nUn
     kap0 = (sHI*( HI0(i)    + HI0(i-1))  &
           + sHeI*( HeI0(i)  + HeI0(i-1)) &
           + sHeII*(HeII0(i) + HeII0(i-1)))*0.5d0*nUn0
        
     !    compute limiter
     if (LType == 1) then       ! rational approx. to LP lim. (LP, 1981)
        D_xl  = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
        D0_xl = c*(2.d0*kap0+R0)/(6.d0*kap0*kap0+3.d0*kap0*R0+R0*R0)
     else if (LType == 2) then  ! Larsen n=2 limiter
        D_xl  = c/sqrt(kap*kap*9.d0 + R*R)
        D0_xl = c/sqrt(kap0*kap0*9.d0 + R0*R0)
     else if (LType == 3) then  ! no limiter
        D_xl  = c/kap/3.d0
        D0_xl = c/kap0/3.d0
     else if (LType == 4) then  ! Zeus limiter
        D_xl  = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
        D0_xl = c*(2.d0*kap0+R0)/(6.d0*kap0*kap0+3.d0*kap0*R0+R0*R0)
     else                       ! standard Levermore-Pomraning (LP, 1981)
        D_xl  = c*(cosh(R/kap)/sinh(R/kap)-kap/R)/R
        D0_xl = c*(cosh(R0/kap0)/sinh(R0/kap0)-kap0/R0)/R0
     endif
     
     !--------------
     ! x-directional limiter, upper face
     ! compute gradients of E0, E
     E0avg  = (E0(i+1) + E0(i))*0.5d0
     E0d_xr = E0(i+1) - E0(i)
     Ed_xr  = E(i+1)  - E(i)
     
     !    compute R for limiter 
     R  = max(dxi *abs(E0d_xr)/E0avg, Rmin)
     R0 = max(dxi0*abs(E0d_xr)/E0avg, Rmin)
     
     !    compute opacity
     kap = (sHI*(  HI(i+1)   + HI(i))  &
          + sHeI*( HeI(i+1)  + HeI(i)) &
          + sHeII*(HeII(i+1) + HeII(i)))*0.5d0*nUn
     kap0 = (sHI*( HI0(i+1)    + HI0(i))  &
           + sHeI*( HeI0(i+1)  + HeI0(i)) &
           + sHeII*(HeII0(i+1) + HeII0(i)))*0.5d0*nUn0
     
     !    compute limiter
     if (LType == 1) then       ! rational approx. to LP lim. (LP, 1981)
        D_xr  = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
        D0_xr = c*(2.d0*kap0+R0)/(6.d0*kap0*kap0+3.d0*kap0*R0+R0*R0)
     else if (LType == 2) then  ! Larsen n=2 limiter
        D_xr  = c/sqrt(kap*kap*9.d0 + R*R)
        D0_xr = c/sqrt(kap0*kap0*9.d0 + R0*R0)
     else if (LType == 3) then  ! no limiter
        D_xr  = c/kap/3.d0
        D0_xr = c/kap0/3.d0
     else if (LType == 4) then  ! Zeus limiter
        D_xr  = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
        D0_xr = c*(2.d0*kap0+R0)/(6.d0*kap0*kap0+3.d0*kap0*R0+R0*R0)
     else                       ! standard Levermore-Pomraning (LP, 1981)
        D_xr  = c*(cosh(R/kap)/sinh(R/kap)-kap/R)/R
        D0_xr = c*(cosh(R0/kap0)/sinh(R0/kap0)-kap0/R0)/R0
     endif
     
     ! opacity values in this cell
     kap = (sHI*HI(i) + sHeI*HeI(i) + sHeII*HeII(i))*nUn
     kap0 = (sHI*HI0(i) + sHeI*HeI0(i) + sHeII*HeII0(i))*nUn0
     
     ! set the matrix entries
     mat(:,i) = (/  &
          -dtfac*dxi*dxi*D_xl, &       ! x-left
          1.d0 + dtfac*(c*kap + dxi*dxi*(D_xl+D_xr)),  &  ! self
          -dtfac*dxi*dxi*D_xr /)       ! x-right
                
     ! set the rhs entries
     rhs(i) = ( (dtfac + dtfac0)*src(i)                         &
                + (1.d0 - dtfac0*c*kap0)*E0(i)                  &
                + dtfac0*dxi0*dxi0*(D0_xr*E0d_xr-D0_xl*E0d_xl)  &
                - (1.d0 + dtfac*c*kap)*E(i)                     &
                + dtfac*dxi*dxi*(D_xr*Ed_xr-D_xl*Ed_xl) )

  enddo

!!$  print *,'  sum(mat) =',sum(mat)
!!$  print *,'  sum(rhs) =',sum(rhs)



  ! update matrix/rhs based on boundary conditions/location
  !    x-left face
  if (xlface == 1) then
     ! Dirichlet
     if (BCXl==1) then
        i = x0s
        mat(1,i) = 0.d0
     ! Neumann
     else if (BCXl==2) then
        i = x0s
        mat(2,i) = mat(2,i) + mat(1,i)
        mat(1,i) = 0.d0
     endif
  end if
  
  !    x-right face
  if (xrface==1) then
     ! Dirichlet
     if (BCXr==1) then
        i = x0e
        mat(3,i) = 0.d0
     ! Neumann
     else if (BCXr==2) then
        i = x0e
        mat(2,i) = mat(2,i) + mat(3,i)
        mat(3,i) = 0.d0
     endif
  endif

  rhsnorm = sum(rhs*rhs)

!!$  print *,'  sum(mat) =',sum(mat)
!!$  print *,'  sum(rhs) =',sum(rhs)
  
  return
end subroutine MFSplit_SetupSystem_1D
!=======================================================================
