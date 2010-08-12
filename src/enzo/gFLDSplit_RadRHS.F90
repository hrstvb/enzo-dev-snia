!=======================================================================
!
! Copyright 2010 Daniel R. Reynolds
!
! This software is released under the terms of the "Enzo Public License"
! in the accompanying LICENSE file.
!
!=======================================================================
subroutine gFLDSplit_RadRHS(rhs, E0, E, Temp, kappa, src, a, adot, ESpec, &
     aUn, lUn, rUn, nUn, rank, dx, dy, dz, Nx, Ny, Nz, NGxl, NGxr, NGyl,  &
     NGyr, NGzl, NGzr, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       June 2010
!  modified:   
!
!  PURPOSE: Computes the ODE RHS for the Grey FLD radiation problem,
!              RHS = Div(D(E)*Grad(E)) - adot/a*E - c*kappa*E + eta + src
!           where D(E) is a nonlinear flux-limiter 
!           depending on E0 (time lagged).  We define the values
!              R_i = |Grad(E)_i|/E,
!           The '_i' subscript implies the gradient in the ith 
!           direction; these quantities are all required at cell faces, 
!           as that is the location of the divergence calculations.
!           With these components, we compute the limiter, 
!                 D_i(E) = c/((3*kappa)**2 + R_i**2).
!
!  INPUTS:
!     E0         - Grey radiation energy density (prev time step)
!     E          - Grey radiation energy density (current guess)
!     Temp       - gas temperature for black-body radiation
!     kappa      - opacity array
!     src        - spatially-dependent radiation source
!     a          - cosmological expansion factor
!     adot       - da/dt
!     ESpec      - flag denoting what type of radiation field we have:
!                  (-1=>monochromatic)
!     *Un        - variable scaling constants
!     rank       - 1, 2 or 3; the dimensionality of the problem
!     dx,dy,dz   - mesh spacing in each direction
!     Nx,Ny,Nz   - active mesh size in each direction
!     NG*l/NG*r  - left/right ghost cells in each direction
!
!     Note: the vector inputs are of size (Nx + NGxl + NGxr) in 
!     the x-direction, others are similar.
!
!  OUTPUT ARGUMENTS: 
!     rhs        - rhs values
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
  integer, intent(in)  :: rank, ESpec
  integer, intent(in)  :: Nx, NGxl, NGxr
  integer, intent(in)  :: Ny, NGyl, NGyr
  integer, intent(in)  :: Nz, NGzl, NGzr
  REALSUB, intent(in)  :: a, adot
  real,    intent(in)  :: dx, dy, dz
  real,    intent(in)  :: aUn, lUn, rUn, nUn
  real,    intent(in)  :: E0(*), E(*), Temp(*), kappa(*), src(*), rhs(*)
  integer, intent(out) :: ier

  !=======================================================================
  

  ! call the apprpriate routine based on rank
  if (rank == 3) then

     call gFLDSplit_RadRHS3D(rhs, E0, E, Temp, kappa, src, a, adot, ESpec, &
          aUn, lUn, rUn, nUn, dx, dy, dz, Nx, Ny, Nz, NGxl, NGxr, NGyl,    &
          NGyr, NGzl, NGzr, ier)

  elseif (rank == 2) then

     call gFLDSplit_RadRHS2D(rhs, E0, E, Temp, kappa, src, a, adot, ESpec, &
          aUn, lUn, rUn, nUn, dx, dy, Nx, Ny, NGxl, NGxr, NGyl, NGyr, ier)

  elseif (rank == 1) then

     call gFLDSplit_RadRHS1D(rhs, E0, E, Temp, kappa, src, a, adot, ESpec, &
          aUn, lUn, rUn, nUn, dx, Nx, NGxl, NGxr, ier)

  else
     write(0,*) 'gFLDSplit_RadRHS error: illegal rank =',rank
  end if

end subroutine gFLDSplit_RadRHS
!=======================================================================






subroutine gFLDSplit_RadRHS3D(rhs, E0, E, Temp, kappa, src, a, adot, ESpec, &
     aUn, lUn, rUn, nUn, dx, dy, dz, Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr,    &
     NGzl, NGzr, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       June 2010
!  modified:   
!
!  PURPOSE: 3D version of the routine
!=======================================================================
#include "fortran.def"
  implicit none
  
  !--------------
  ! argument declarations
  integer, intent(in) :: ESpec
  integer, intent(in) :: Nx, NGxl, NGxr
  integer, intent(in) :: Ny, NGyl, NGyr
  integer, intent(in) :: Nz, NGzl, NGzr
  REALSUB, intent(in) :: a, adot
  real,    intent(in) :: dx, dy, dz
  real,    intent(in) :: aUn, lUn, rUn, nUn
  real, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr), intent(in) &
                      :: E0, E, kappa, src, Temp
  real,    intent(out) :: rhs(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr) 
  integer, intent(out) :: ier

  !--------------
  ! locals
  integer :: i, j, k
  real*8  :: kap, eta, c, pi, StBz, dxi, dyi, dzi, afac
  real*8  :: D, E0d, Ed, E0avg, R, Rmin

!=======================================================================
  
  ! initialize outputs to zero, flag to success
  rhs = 0.d0
  ier = 1

  ! set shortcut values
  if (ESpec == -1) then
     afac = 0.d0
  else
     afac = adot/a
  endif
  dxi   = a/dx/lUn
  dyi   = a/dy/lUn
  dzi   = a/dz/lUn
  c     = 2.99792458d10     ! speed of light [cm/s]
  pi    = 4.d0*atan(1.d0)
  Rmin  = 1.d-20
  StBz  = 5.6704d-5         ! Stefan-Boltzmann constant [ergs/(s cm^2 K^4)]


  ! iterate over the active domain
  do k=1,Nz,1
     do j=1,Ny,1
        do i=1,Nx,1

           !--------------
           ! this cell
           !    opacity values
           kap = kappa(i,j,k)*nUn
           !    black-body radiation (if applicable; otherwise Temp=0)
           eta = 4.d0*kap*StBz/rUn*Temp(i,j,k)**4
           !    initialize the rhs value
           rhs(i,j,k) = src(i,j,k)/rUn + eta - (afac + c*kap)*E(i,j,k) 

           !--------------
           ! z-direction, lower face
           !    compute gradients of E0, E
           E0d   = E0(i,j,k) - E0(i,j,k-1)
           Ed    = E(i,j,k) - E(i,j,k-1)
           E0avg = (E0(i,j,k) + E0(i,j,k-1))*0.5d0
           !    compute R for limiters
           R = max(dzi*abs(E0d)/E0avg, Rmin)
           !    compute average opacity over face
           kap = (kappa(i,j,k) + kappa(i,j,k-1))*0.5d0*nUn
           !    compute limiter
!!$           D = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
           D = c/sqrt(9.d0*kap*kap + R*R)
           !    update rhs
           rhs(i,j,k) = rhs(i,j,k) - dzi*dzi*D*Ed

           !--------------
           ! y-direction, lower face
           !    compute gradients of E0, E
           E0d   = E0(i,j,k) - E0(i,j-1,k)
           Ed    = E(i,j,k) - E(i,j-1,k)
           E0avg = (E0(i,j,k) + E0(i,j-1,k))*0.5d0
           !    compute R for limiters
           R = max(dyi*abs(E0d)/E0avg, Rmin)
           !    compute average opacity over face
           kap = (kappa(i,j,k) + kappa(i,j-1,k))*0.5d0*nUn
           !    compute limiter
!!$           D = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
           D = c/sqrt(9.d0*kap*kap + R*R)
           !    update rhs
           rhs(i,j,k) = rhs(i,j,k) - dyi*dyi*D*Ed

           !--------------
           ! x-direction, lower face
           !    compute gradients of E0, E
           E0d   = E0(i,j,k) - E0(i-1,j,k)
           Ed    = E(i,j,k) - E(i-1,j,k)
           E0avg = (E0(i,j,k) + E0(i-1,j,k))*0.5d0
           !    compute R for limiters
           R = max(dxi*abs(E0d)/E0avg, Rmin)
           !    compute average opacity over face
           kap = (kappa(i,j,k) + kappa(i-1,j,k))*0.5d0*nUn
           !    compute limiter
!!$           D = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
           D = c/sqrt(9.d0*kap*kap + R*R)
           !    update rhs
           rhs(i,j,k) = rhs(i,j,k) - dxi*dxi*D*Ed

           !--------------
           ! x-direction, upper face
           !    compute gradients of E0, E
           E0d   = E0(i+1,j,k) - E0(i,j,k)
           Ed    = E(i+1,j,k) - E(i,j,k)
           E0avg = (E0(i+1,j,k) + E0(i,j,k))*0.5d0
           !    compute R for limiters
           R = max(dxi*abs(E0d)/E0avg, Rmin)
           !    compute average opacity over face
           kap = (kappa(i,j,k) + kappa(i+1,j,k))*0.5d0*nUn
           !    compute limiter
!!$           D = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
           D = c/sqrt(9.d0*kap*kap + R*R)
           !    update rhs
           rhs(i,j,k) = rhs(i,j,k) + dxi*dxi*D*Ed

           !--------------
           ! y-direction, upper face
           !    compute gradients of E0, E
           E0d   = E0(i,j+1,k) - E0(i,j,k)
           Ed    = E(i,j+1,k) - E(i,j,k)
           E0avg = (E0(i,j+1,k) + E0(i,j,k))*0.5d0
           !    compute R for limiters
           R = max(dyi*abs(E0d)/E0avg, Rmin)
           !    compute average opacity over face
           kap = (kappa(i,j,k) + kappa(i,j+1,k))*0.5d0*nUn
           !    compute limiter
!!$           D = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
           D = c/sqrt(9.d0*kap*kap + R*R)
           !    update rhs
           rhs(i,j,k) = rhs(i,j,k) + dyi*dyi*D*Ed

           !--------------
           ! z-direction, upper face
           !    compute gradients of E0, E
           E0d   = E0(i,j,k+1) - E0(i,j,k)
           Ed    = E(i,j,k+1) - E(i,j,k)
           E0avg = (E0(i,j,k+1) + E0(i,j,k))*0.5d0
           !    compute R for limiters
           R = max(dzi*abs(E0d)/E0avg, Rmin)
           !    compute average opacity over face
           kap = (kappa(i,j,k) + kappa(i,j,k+1))*0.5d0*nUn
           !    compute limiter
!!$           D = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
           D = c/sqrt(9.d0*kap*kap + R*R)
           !    update rhs
           rhs(i,j,k) = rhs(i,j,k) + dzi*dzi*D*Ed

        enddo
     enddo
  enddo

  return
end subroutine gFLDSplit_RadRHS3D
!=======================================================================






subroutine gFLDSplit_RadRHS2D(rhs, E0, E, Temp, kappa, src, a, adot, ESpec, &
     aUn, lUn, rUn, nUn, dx, dy, Nx, Ny, NGxl, NGxr, NGyl, NGyr, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       June 2010
!  modified:   
!
!  PURPOSE: 2D version of the routine
!=======================================================================
#include "fortran.def"
  implicit none
  
  !--------------
  ! argument declarations
  integer, intent(in) :: ESpec
  integer, intent(in) :: Nx, NGxl, NGxr
  integer, intent(in) :: Ny, NGyl, NGyr
  REALSUB, intent(in) :: a, adot
  real,    intent(in) :: dx, dy
  real,    intent(in) :: aUn, lUn, rUn, nUn
  real, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr), intent(in) &
                      :: E0, E, kappa, src, Temp
  real,    intent(out) :: rhs(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr) 
  integer, intent(out) :: ier

  !--------------
  ! locals
  integer :: i, j
  real*8  :: kap, eta, c, pi, StBz, dxi, dyi, afac
  real*8  :: D, E0d, Ed, E0avg, R, Rmin

!=======================================================================
  
  ! initialize outputs to zero, flag to success
  rhs = 0.d0
  ier = 1

  ! set shortcut values
  if (ESpec == -1) then
     afac = 0.d0
  else
     afac = adot/a
  endif
  dxi   = a/dx/lUn
  dyi   = a/dy/lUn
  c     = 2.99792458d10     ! speed of light [cm/s]
  pi    = 4.d0*atan(1.d0)
  Rmin  = 1.d-20
  StBz  = 5.6704d-5         ! Stefan-Boltzmann constant [ergs/(s cm^2 K^4)]


  ! iterate over the active domain
  do j=1,Ny,1
     do i=1,Nx,1

        !--------------
        ! this cell
        !    opacity values
        kap = kappa(i,j)*nUn
        !    black-body radiation (if applicable; otherwise Temp=0)
        eta = 4.d0*kap*StBz/rUn*Temp(i,j)**4
        !    initialize the rhs value
        rhs(i,j) = src(i,j)/rUn + eta - (afac + c*kap)*E(i,j) 
        
        !--------------
        ! y-direction, lower face
        !    compute gradients of E0, E
        E0d   = E0(i,j) - E0(i,j-1)
        Ed    = E(i,j) - E(i,j-1)
        E0avg = (E0(i,j) + E0(i,j-1))*0.5d0
        !    compute R for limiters
        R = max(dyi*abs(E0d)/E0avg, Rmin)
        !    compute average opacity over face
        kap = (kappa(i,j) + kappa(i,j-1))*0.5d0*nUn
        !    compute limiter
!!$           D = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
        D = c/sqrt(9.d0*kap*kap + R*R)
        !    update rhs
        rhs(i,j) = rhs(i,j) - dyi*dyi*D*Ed

        !--------------
        ! x-direction, lower face
        !    compute gradients of E0, E
        E0d   = E0(i,j) - E0(i-1,j)
        Ed    = E(i,j) - E(i-1,j)
        E0avg = (E0(i,j) + E0(i-1,j))*0.5d0
        !    compute R for limiters
        R = max(dxi*abs(E0d)/E0avg, Rmin)
        !    compute average opacity over face
        kap = (kappa(i,j) + kappa(i-1,j))*0.5d0*nUn
        !    compute limiter
!!$           D = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
        D = c/sqrt(9.d0*kap*kap + R*R)
        !    update rhs
        rhs(i,j) = rhs(i,j) - dxi*dxi*D*Ed
        
        !--------------
        ! x-direction, upper face
        !    compute gradients of E0, E
        E0d   = E0(i+1,j) - E0(i,j)
        Ed    = E(i+1,j) - E(i,j)
        E0avg = (E0(i+1,j) + E0(i,j))*0.5d0
        !    compute R for limiters
        R = max(dxi*abs(E0d)/E0avg, Rmin)
        !    compute average opacity over face
        kap = (kappa(i,j) + kappa(i+1,j))*0.5d0*nUn
        !    compute limiter
!!$           D = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
        D = c/sqrt(9.d0*kap*kap + R*R)
        !    update rhs
        rhs(i,j) = rhs(i,j) + dxi*dxi*D*Ed
        
        !--------------
        ! y-direction, upper face
        !    compute gradients of E0, E
        E0d   = E0(i,j+1) - E0(i,j)
        Ed    = E(i,j+1) - E(i,j)
        E0avg = (E0(i,j+1) + E0(i,j))*0.5d0
        !    compute R for limiters
        R = max(dyi*abs(E0d)/E0avg, Rmin)
        !    compute average opacity over face
        kap = (kappa(i,j) + kappa(i,j+1))*0.5d0*nUn
        !    compute limiter
!!$           D = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
        D = c/sqrt(9.d0*kap*kap + R*R)
        !    update rhs
        rhs(i,j) = rhs(i,j) + dyi*dyi*D*Ed
        
     enddo
  enddo
  
  return
end subroutine gFLDSplit_RadRHS2D
!=======================================================================






subroutine gFLDSplit_RadRHS1D(rhs, E0, E, Temp, kappa, src, a, adot, &
     ESpec, aUn, lUn, rUn, nUn, dx, Nx, NGxl, NGxr, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       June 2010
!  modified:   
!
!  PURPOSE: 1D version of the routine
!=======================================================================
#include "fortran.def"
  implicit none
  
  !--------------
  ! argument declarations
  integer, intent(in) :: ESpec
  integer, intent(in) :: Nx, NGxl, NGxr
  REALSUB, intent(in) :: a, adot
  real,    intent(in) :: dx
  real,    intent(in) :: aUn, lUn, rUn, nUn
  real, dimension(1-NGxl:Nx+NGxr), intent(in) :: E0, E, kappa, src, Temp
  real,    intent(out) :: rhs(1-NGxl:Nx+NGxr) 
  integer, intent(out) :: ier

  !--------------
  ! locals
  integer :: i
  real*8  :: kap, eta, c, pi, StBz, dxi, afac
  real*8  :: D, E0d, Ed, E0avg, R, Rmin

!=======================================================================
  
  ! initialize outputs to zero, flag to success
  rhs = 0.d0
  ier = 1

  ! set shortcut values
  if (ESpec == -1) then
     afac = 0.d0
  else
     afac = adot/a
  endif
  dxi   = a/dx/lUn
  c     = 2.99792458d10     ! speed of light [cm/s]
  pi    = 4.d0*atan(1.d0)
  Rmin  = 1.d-20
  StBz  = 5.6704d-5         ! Stefan-Boltzmann constant [ergs/(s cm^2 K^4)]


  ! iterate over the active domain
  do i=1,Nx,1
     
     !--------------
     ! this cell
     !    opacity values
     kap = kappa(i)*nUn
     !    black-body radiation (if applicable; otherwise Temp=0)
     eta = 4.d0*kap*StBz/rUn*Temp(i)**4
     !    initialize the rhs value
     rhs(i) = src(i)/rUn + eta - (afac + c*kap)*E(i) 
     
     !--------------
     ! x-direction, lower face
     !    compute gradients of E0, E
     E0d   = E0(i) - E0(i-1)
     Ed    = E(i) - E(i-1)
     E0avg = (E0(i) + E0(i-1))*0.5d0
     !    compute R for limiters
     R = max(dxi*abs(E0d)/E0avg, Rmin)
     !    compute average opacity over face
     kap = (kappa(i) + kappa(i-1))*0.5d0*nUn
     !    compute limiter
!!$           D = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
     D = c/sqrt(9.d0*kap*kap + R*R)
     !    update rhs
     rhs(i) = rhs(i) - dxi*dxi*D*Ed
     
     !--------------
     ! x-direction, upper face
     !    compute gradients of E0, E
     E0d   = E0(i+1) - E0(i)
     Ed    = E(i+1) - E(i)
     E0avg = (E0(i+1) + E0(i))*0.5d0
     !    compute R for limiters
     R = max(dxi*abs(E0d)/E0avg, Rmin)
     !    compute average opacity over face
     kap = (kappa(i) + kappa(i+1))*0.5d0*nUn
     !    compute limiter
!!$           D = c*(2.d0*kap+R)/(6.d0*kap*kap+3.d0*kap*R+R*R)
     D = c/sqrt(9.d0*kap*kap + R*R)
     !    update rhs
     rhs(i) = rhs(i) + dxi*dxi*D*Ed
        
  enddo
  
  return
end subroutine gFLDSplit_RadRHS1D
!=======================================================================
