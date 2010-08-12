!=======================================================================
!
! Copyright 2009 Daniel R. Reynolds
!
! This software is released under the terms of the "Enzo Public License"
! in the accompanying LICENSE file.
!
!=======================================================================
subroutine MFProb_RadResid(res_E1, res_E2, res_E3, E1, E10, E2, E20, E3, &
     E30, HI, HI0, HeI, HeI0, HeII, HeII0, src_E1, src_E2, src_E3, LTyp, &
     LImp, dt, theta, a, a0, adot, adot0, aUn, lUn, lUn0, rUn, rUn0,     &
     nUn, nUn0, dx, dy, dz, Nchem, Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr,   &
     NGzl, NGzr, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       August 2009
!  modified1:  
!
!  PURPOSE: Computes the nonlinear residual for the multi-frequency 
!           radiation problem,
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
!           Each of the above three forms has relative merits:
!             - The original LP formulation has been well-tested in 
!               the community; however it involves subtracting two
!               large numbers when R is close to 0, allowing for 
!               possibly large roundoff errors.  Additionally, the 
!               coth = sinh/cosh function may be costly to perform 
!               repeatedly.
!             - The rational approximation alleviates the cost of 
!               intrinsic functions and catastrophic floating-point 
!               cancellation errors, but introduces a 7.2% error 
!               away from the original limiter.
!             - The new approximation also alleviates any catastrophic 
!               floating-point cancellation errors, and introduces 
!               only a 4.8% error away from the original limiter; 
!               however it still involves the use of possibly-expensive 
!               intrinsic functions.
!
!  INPUTS:
!     E*,E*0     - Radiation energy density at each frequency (new & old times)
!     HI,HI0     - Hydrogen I density (new & old times)
!     HeI,HeI0   - Helium I density (new & old times)
!     HeII,HeII0 - Helium II density (new & old times)
!     src_E*     - emissivity sources for each frequency
!     LTyp       - integer flag denoting type of flux limiter:
!                       0 -> standard Levermore-Pomraning lim. (LP, 1981)
!                       1 -> rational approx. to LP lim. (LP, 1981)
!                       2 -> Reynolds approx to LP lim.
!                       3 -> turns off limiter (constant of 1/3)
!                       4 -> Zeus limiter
!     LImp       - integer flag denoting implicitness of flux limiter:
!                       0 -> fully lagged to previous time step
!                       1 -> fully lagged to previous newton iterate
!                       2 -> lag only temperature dependence
!     dt         - time step size
!     theta      - time integration parameter
!     a,a0       - cosmological expansion parameter (new & old times)
!     adot,adot0 - da/dt (new & old times)
!     *Un,*Un0   - variable scaling constants (new & old times)
!     dx,dy,dz   - mesh spacing in each direction
!     Nx,Ny,Nz   - active mesh size in each direction
!     NG*l/NG*r  - left/right ghost cells in each direction
!
!     Note: the vector inputs are of size (Nx + NGxl + NGxr) in 
!     the x-direction, others are similar.
!
!  OUTPUT ARGUMENTS: 
!     res_E1     - residual for frequency 1 (in internal Enzo units)
!     res_E2     - residual for frequency 2 (in internal Enzo units)
!     res_E3     - residual for frequency 3 (in internal Enzo units)
!     ier        - success/failure output flag (0->failure, 1->success)
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
  integer, intent(in)  :: LTyp, LImp, Nchem
  integer, intent(in)  :: Nx, NGxl, NGxr
  integer, intent(in)  :: Ny, NGyl, NGyr
  integer, intent(in)  :: Nz, NGzl, NGzr
  integer, intent(out) :: ier
  REALSUB, intent(in)  :: a, adot, a0, adot0
  real, intent(in) :: aUn, lUn, lUn0, rUn, rUn0, nUn, nUn0
  real, intent(in) :: dt, theta, dx, dy, dz
  real, intent(in),  dimension(*) :: E1, E2, E3, E10, E20, E30
  real, intent(in),  dimension(*) :: HI, HeI, HeII, HI0, HeI0, HeII0
  real, intent(in),  dimension(*) :: src_E1, src_E2, src_E3
  real, intent(out), dimension(*) :: res_E1, res_E2, res_E3

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
 

  ! supply arguments based on Nchem
  if (Nchem == 1) then

     !!!!!!! E1 evaluation !!!!!!!
     !   set shortcut variables
     call MFProb_HICrossSection(sHI, nu0_HI, 1)
     sHeI  = 0.d0
     sHeII = 0.d0

     ! call the appropriate dimension-specific routine
     if ((Nz == 1) .and. (Ny == 1)) then
        call MFProb_RadResid_1D(res_E1, E1, E10, HI, HI0, HI, HI0, &
             HI, HI0, src_E1, LTyp, LImp, dt, theta, sHI, sHeI,  &
             sHeII, a, a0, adot, adot0, aUn, lUn, lUn0, rUn, rUn0, nUn,  &
             nUn0, dx, Nx, NGxl, NGxr, ier)
     elseif (Nz == 1) then
        call MFProb_RadResid_2D(res_E1, E1, E10, HI, HI0, HI, HI0, &
             HI, HI0, src_E1, LTyp, LImp, dt, theta, sHI, sHeI,  &
             sHeII, a, a0, adot, adot0, aUn, lUn, lUn0, rUn, rUn0, nUn,  &
             nUn0, dx, dy, Nx, Ny, NGxl, NGxr, NGyl, NGyr, ier)
     else
        call MFProb_RadResid_3D(res_E1, E1, E10, HI, HI0, HI, HI0, &
             HI, HI0, src_E1, LTyp, LImp, dt, theta, sHI, sHeI,  &
             sHeII, a, a0, adot, adot0, aUn, lUn, lUn0, rUn, rUn0, nUn,  &
             nUn0, dx, dy, dz, Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr, NGzl, &
             NGzr, ier)
     end if


     !!!!!!! E2 evaluation !!!!!!!
     !   set shortcut variables
     call MFProb_HICrossSection(sHI, nu0_HeI, 1)

     ! call the appropriate dimension-specific routine
     if ((Nz == 1) .and. (Ny == 1)) then
        call MFProb_RadResid_1D(res_E2, E2, E20, HI, HI0, HI, HI0, &
             HI, HI0, src_E2, LTyp, LImp, dt, theta, sHI, sHeI,  &
             sHeII, a, a0, adot, adot0, aUn, lUn, lUn0, rUn, rUn0, nUn,  &
             nUn0, dx, Nx, NGxl, NGxr, ier)
     elseif (Nz == 1) then
        call MFProb_RadResid_2D(res_E2, E2, E20, HI, HI0, HI, HI0, &
             HI, HI0, src_E2, LTyp, LImp, dt, theta, sHI, sHeI,  &
             sHeII, a, a0, adot, adot0, aUn, lUn, lUn0, rUn, rUn0, nUn,  &
             nUn0, dx, dy, Nx, Ny, NGxl, NGxr, NGyl, NGyr, ier)
     else
        call MFProb_RadResid_3D(res_E2, E2, E20, HI, HI0, HI, HI0, &
             HI, HI0, src_E2, LTyp, LImp, dt, theta, sHI, sHeI,  &
             sHeII, a, a0, adot, adot0, aUn, lUn, lUn0, rUn, rUn0, nUn,  &
             nUn0, dx, dy, dz, Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr, NGzl, &
             NGzr, ier)
     end if


     !!!!!!! E3 evaluation !!!!!!!
     !   set shortcut variables
     call MFProb_HICrossSection(sHI, nu0_HeII, 1)

     ! call the appropriate dimension-specific routine
     if ((Nz == 1) .and. (Ny == 1)) then
        call MFProb_RadResid_1D(res_E3, E3, E30, HI, HI0, HI, HI0, &
             HI, HI0, src_E3, LTyp, LImp, dt, theta, sHI, sHeI,  &
             sHeII, a, a0, adot, adot0, aUn, lUn, lUn0, rUn, rUn0, nUn,  &
             nUn0, dx, Nx, NGxl, NGxr, ier)
     elseif (Nz == 1) then
        call MFProb_RadResid_2D(res_E3, E3, E30, HI, HI0, HI, HI0, &
             HI, HI0, src_E3, LTyp, LImp, dt, theta, sHI, sHeI,  &
             sHeII, a, a0, adot, adot0, aUn, lUn, lUn0, rUn, rUn0, nUn,  &
             nUn0, dx, dy, Nx, Ny, NGxl, NGxr, NGyl, NGyr, ier)
     else
        call MFProb_RadResid_3D(res_E3, E3, E30, HI, HI0, HI, HI0, &
             HI, HI0, src_E3, LTyp, LImp, dt, theta, sHI, sHeI,  &
             sHeII, a, a0, adot, adot0, aUn, lUn, lUn0, rUn, rUn0, nUn,  &
             nUn0, dx, dy, dz, Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr, NGzl, &
             NGzr, ier)
     end if

  else  ! Nchem == 3

     !!!!!!! E1 evaluation !!!!!!!
     !   set shortcut variables
     call MFProb_HICrossSection(sHI, nu0_HI, 1)
     sHeI  = 0.d0
     sHeII = 0.d0

     ! call the appropriate dimension-specific routine
     if ((Nz == 1) .and. (Ny == 1)) then
        call MFProb_RadResid_1D(res_E1, E1, E10, HI, HI0, HeI, HeI0, &
             HeII, HeII0, src_E1, LTyp, LImp, dt, theta, sHI, sHeI,  &
             sHeII, a, a0, adot, adot0, aUn, lUn, lUn0, rUn, rUn0, nUn,  &
             nUn0, dx, Nx, NGxl, NGxr, ier)
     elseif (Nz == 1) then
        call MFProb_RadResid_2D(res_E1, E1, E10, HI, HI0, HeI, HeI0, &
             HeII, HeII0, src_E1, LTyp, LImp, dt, theta, sHI, sHeI,  &
             sHeII, a, a0, adot, adot0, aUn, lUn, lUn0, rUn, rUn0, nUn,  &
             nUn0, dx, dy, Nx, Ny, NGxl, NGxr, NGyl, NGyr, ier)
     else
        call MFProb_RadResid_3D(res_E1, E1, E10, HI, HI0, HeI, HeI0, &
             HeII, HeII0, src_E1, LTyp, LImp, dt, theta, sHI, sHeI,  &
             sHeII, a, a0, adot, adot0, aUn, lUn, lUn0, rUn, rUn0, nUn,  &
             nUn0, dx, dy, dz, Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr, NGzl, &
             NGzr, ier)
     end if


     !!!!!!! E2 evaluation !!!!!!!
     !   set shortcut variables
     call MFProb_HICrossSection(sHI, nu0_HeI, 1)
     call MFProb_HeICrossSection(sHeI, nu0_HeI, 1)

     ! call the appropriate dimension-specific routine
     if ((Nz == 1) .and. (Ny == 1)) then
        call MFProb_RadResid_1D(res_E2, E2, E20, HI, HI0, HeI, HeI0, &
             HeII, HeII0, src_E2, LTyp, LImp, dt, theta, sHI, sHeI,  &
             sHeII, a, a0, adot, adot0, aUn, lUn, lUn0, rUn, rUn0, nUn,  &
             nUn0, dx, Nx, NGxl, NGxr, ier)
     elseif (Nz == 1) then
        call MFProb_RadResid_2D(res_E2, E2, E20, HI, HI0, HeI, HeI0, &
             HeII, HeII0, src_E2, LTyp, LImp, dt, theta, sHI, sHeI,  &
             sHeII, a, a0, adot, adot0, aUn, lUn, lUn0, rUn, rUn0, nUn,  &
             nUn0, dx, dy, Nx, Ny, NGxl, NGxr, NGyl, NGyr, ier)
     else
        call MFProb_RadResid_3D(res_E2, E2, E20, HI, HI0, HeI, HeI0, &
             HeII, HeII0, src_E2, LTyp, LImp, dt, theta, sHI, sHeI,  &
             sHeII, a, a0, adot, adot0, aUn, lUn, lUn0, rUn, rUn0, nUn,  &
             nUn0, dx, dy, dz, Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr, NGzl, &
             NGzr, ier)
     end if


     !!!!!!! E3 evaluation !!!!!!!
     !   set shortcut variables
     call MFProb_HICrossSection(sHI, nu0_HeII, 1)
     call MFProb_HeICrossSection( sHeI,  nu0_HeII, 1)
     call MFProb_HeIICrossSection(sHeII, nu0_HeII, 1)

     ! call the appropriate dimension-specific routine
     if ((Nz == 1) .and. (Ny == 1)) then
        call MFProb_RadResid_1D(res_E3, E3, E30, HI, HI0, HeI, HeI0, &
             HeII, HeII0, src_E3, LTyp, LImp, dt, theta, sHI, sHeI,  &
             sHeII, a, a0, adot, adot0, aUn, lUn, lUn0, rUn, rUn0, nUn,  &
             nUn0, dx, Nx, NGxl, NGxr, ier)
     elseif (Nz == 1) then
        call MFProb_RadResid_2D(res_E3, E3, E30, HI, HI0, HeI, HeI0, &
             HeII, HeII0, src_E3, LTyp, LImp, dt, theta, sHI, sHeI,  &
             sHeII, a, a0, adot, adot0, aUn, lUn, lUn0, rUn, rUn0, nUn,  &
             nUn0, dx, dy, Nx, Ny, NGxl, NGxr, NGyl, NGyr, ier)
     else
        call MFProb_RadResid_3D(res_E3, E3, E30, HI, HI0, HeI, HeI0, &
             HeII, HeII0, src_E3, LTyp, LImp, dt, theta, sHI, sHeI,  &
             sHeII, a, a0, adot, adot0, aUn, lUn, lUn0, rUn, rUn0, nUn,  &
             nUn0, dx, dy, dz, Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr, NGzl, &
             NGzr, ier)
     end if

  endif  ! Nchem

  return
end subroutine MFProb_RadResid
!=======================================================================






subroutine MFProb_RadResid_3D(res, E, E0, HI, HI0, HeI, HeI0, HeII,    &
     HeII0, src, LTyp, LImp, dt, theta, sHI, sHeI, sHeII, a, a0, adot, &
     adot0, aUn, lUn, lUn0, rUn, rUn0, nUn, nUn0, dx, dy, dz, Nx, Ny,  &
     Nz, NGxl, NGxr, NGyl, NGyr, NGzl, NGzr, ier)
!=======================================================================
!  PURPOSE: 3D version of the routine
!=======================================================================
#include "fortran.def"
  implicit none

  !--------------
  ! argument declarations
  integer, intent(in)  :: LTyp, LImp
  integer, intent(in)  :: Nx, NGxl, NGxr
  integer, intent(in)  :: Ny, NGyl, NGyr
  integer, intent(in)  :: Nz, NGzl, NGzr
  integer, intent(out) :: ier
  REALSUB, intent(in)  :: a, adot, a0, adot0
  real, intent(in) :: aUn, lUn, lUn0, rUn, rUn0, nUn, nUn0, sHI, sHeI, sHeII
  real, intent(in) :: dt, theta, dx, dy, dz
  real, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr), &
       intent(in) :: E, E0, src
  real, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr), &
       intent(in) :: HI, HeI, HeII, HI0, HeI0, HeII0
  real, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr), &
       intent(out) :: res

  !--------------
  ! locals
  integer  :: i, j, k
  real :: dtfac, dtfac0, etafac, c, pi
  real :: dxi, dxi0, dyi, dyi0, dzi, dzi0
  real :: D_xl, D0_xl, D_xr, D0_xr, Ed_xl, E0d_xl, Ed_xr, E0d_xr
  real :: D_yl, D0_yl, D_yr, D0_yr, Ed_yl, E0d_yl, Ed_yr, E0d_yr
  real :: D_zl, D0_zl, D_zr, D0_zr, Ed_zl, E0d_zl, Ed_zr, E0d_zr
  real :: Eavg, E0avg, R, R0, kap, kap0, lag, Rmin
  real :: kap_xl, kap0_xl, kap_xr, kap0_xr
  real :: kap_yl, kap0_yl, kap_yr, kap0_yr
  real :: kap_zl, kap0_zl, kap_zr, kap0_zr

  !=======================================================================

  ! initialize output to zero, flag to success
  res = 0.d0
  ier = 1

  ! set shortcut values 
  c = 2.99792458d10             ! speed of light [cm/s]
  pi = 4.d0*atan(1.d0)
  Rmin = 1.0d-20
  dtfac  = dt*theta
  dtfac0 = dt*(1.d0-theta)
  etafac = dt*8.d0*pi/(rUn+rUn0)
  dxi    = a/dx/lUn
  dyi    = a/dy/lUn
  dzi    = a/dz/lUn
  dxi0   = a0/dx/lUn0
  dyi0   = a0/dy/lUn0
  dzi0   = a0/dz/lUn0

  ! set lag constant to account for lagging E in limiter evaluation
  lag = 0.d0
  if (LImp == 0)  lag = 1.d0

  ! compute radiation energy gradient over domain
  do k=1,Nz
     do j=1,Ny
        do i=1,Nx

           ! compute old opacities at faces and cell center
           kap0_zl = (sHI*(  HI0(i,j,k)  +HI0(i,j,k-1))  &
                    + sHeI*( HeI0(i,j,k) +HeI0(i,j,k-1)) &
                    + sHeII*(HeII0(i,j,k)+HeII0(i,j,k-1)))*0.5d0*nUn0
           kap0_yl = (sHI*(  HI0(i,j,k)  +HI0(i,j-1,k))  &
                    + sHeI*( HeI0(i,j,k) +HeI0(i,j-1,k)) &
                    + sHeII*(HeII0(i,j,k)+HeII0(i,j-1,k)))*0.5d0*nUn0
           kap0_xl = (sHI*(  HI0(i,j,k)  +HI0(i-1,j,k))  &
                    + sHeI*( HeI0(i,j,k) +HeI0(i-1,j,k)) &
                    + sHeII*(HeII0(i,j,k)+HeII0(i-1,j,k)))*0.5d0*nUn0
           kap0_xr = (sHI*(  HI0(i,j,k)  +HI0(i+1,j,k))  &
                    + sHeI*( HeI0(i,j,k) +HeI0(i+1,j,k)) &
                    + sHeII*(HeII0(i,j,k)+HeII0(i+1,j,k)))*0.5d0*nUn0
           kap0_yr = (sHI*(  HI0(i,j,k)  +HI0(i,j+1,k))  &
                    + sHeI*( HeI0(i,j,k) +HeI0(i,j+1,k)) &
                    + sHeII*(HeII0(i,j,k)+HeII0(i,j+1,k)))*0.5d0*nUn0
           kap0_zr = (sHI*(  HI0(i,j,k)  +HI0(i,j+1,k))  &
                    + sHeI*( HeI0(i,j,k) +HeI0(i,j+1,k)) &
                    + sHeII*(HeII0(i,j,k)+HeII0(i,j+1,k)))*0.5d0*nUn0
           kap0 = (sHI*HI0(i,j,k) + sHeI*HeI0(i,j,k) + sHeII*HeII0(i,j,k))*nUn0

           ! compute new opacities at faces and cell center
           if (lag == 1) then
              kap_zl = kap0_zl
              kap_yl = kap0_yl
              kap_xl = kap0_xl
              kap_xr = kap0_xr
              kap_yr = kap0_yr
              kap_zr = kap0_zr
           else
              kap_zl = (sHI*(  HI(i,j,k)  +HI(i,j,k-1))  &
                      + sHeI*( HeI(i,j,k) +HeI(i,j,k-1)) &
                      + sHeII*(HeII(i,j,k)+HeII(i,j,k-1)))*0.5d0*nUn
              kap_yl = (sHI*(  HI(i,j,k)  +HI(i,j-1,k))  &
                      + sHeI*( HeI(i,j,k) +HeI(i,j-1,k)) &
                      + sHeII*(HeII(i,j,k)+HeII(i,j-1,k)))*0.5d0*nUn
              kap_xl = (sHI*(  HI(i,j,k)  +HI(i-1,j,k))  &
                      + sHeI*( HeI(i,j,k) +HeI(i-1,j,k)) &
                      + sHeII*(HeII(i,j,k)+HeII(i-1,j,k)))*0.5d0*nUn
              kap_xr = (sHI*(  HI(i,j,k)  +HI(i+1,j,k))  &
                      + sHeI*( HeI(i,j,k) +HeI(i+1,j,k)) &
                      + sHeII*(HeII(i,j,k)+HeII(i+1,j,k)))*0.5d0*nUn
              kap_yr = (sHI*(  HI(i,j,k)  +HI(i,j+1,k))  &
                      + sHeI*( HeI(i,j,k) +HeI(i,j+1,k)) &
                      + sHeII*(HeII(i,j,k)+HeII(i,j+1,k)))*0.5d0*nUn
              kap_zr = (sHI*(  HI(i,j,k)  +HI(i,j+1,k))  &
                      + sHeI*( HeI(i,j,k) +HeI(i,j+1,k)) &
                      + sHeII*(HeII(i,j,k)+HeII(i,j+1,k)))*0.5d0*nUn
           endif
           kap = (sHI*HI(i,j,k) + sHeI*HeI(i,j,k) + sHeII*HeII(i,j,k))*nUn


           !--------------
           ! z-directional limiter, lower face
           !    compute gradients, face values
           Ed_zl  = (E(i,j,k)  - E(i,j,k-1))
           E0d_zl = (E0(i,j,k) - E0(i,j,k-1))
           Eavg   = (E(i,j,k)  + E(i,j,k-1) )*0.5d0
           E0avg  = (E0(i,j,k) + E0(i,j,k-1))*0.5d0

           !    compute R for limiters
           R  = max(dzi *abs(Ed_zl )/Eavg,  Rmin)
           R0 = max(dzi0*abs(E0d_zl)/E0avg, Rmin)
           R  = lag*R0 + (1.d0-lag)*R

           !    compute limiter
           if (LTyp == 1) then       ! rational approx. to LP lim. (LP, 1981)
              D_zl  = c*(2.d0*kap_zl+R)/(6.d0*kap_zl*kap_zl+3.d0*kap_zl*R+R*R)
              D0_zl = c*(2.d0*kap0_zl+R0)/(6.d0*kap0_zl*kap0_zl+3.d0*kap0_zl*R0+R0*R0)
           else if (LTyp == 2) then  ! Reynolds approx to LP lim.
              D_zl  = 2.d0*c/pi*datan(R*pi/6.d0/kap_zl)/R
              D0_zl = 2.d0*c/pi*datan(R0*pi/6.d0/kap0_zl)/R0
           else if (LTyp == 3) then  ! no limiter
              D_zl  = c/kap_zl/3.d0
              D0_zl = c/kap0_zl/3.d0
           else if (LTyp == 4) then  ! Zeus limiter
              D_zl  = c*(2.d0*kap_zl+R)/(6.d0*kap_zl*kap_zl+3.d0*kap_zl*R+R*R)
              D0_zl = c*(2.d0*kap0_zl+R0)/(6.d0*kap0_zl*kap0_zl+3.d0*kap0_zl*R0+R0*R0)
           else                       ! standard Levermore-Pomraning (LP, 1981)
              D_zl  = c*(cosh(R/kap_zl)/sinh(R/kap_zl)-kap_zl/R)/R
              D0_zl = c*(cosh(R0/kap0_zl)/sinh(R0/kap0_zl)-kap0_zl/R0)/R0
           endif

           !--------------
           ! y-directional limiter, lower face
           !    compute gradients, face values
           Ed_yl  = (E(i,j,k)  - E(i,j-1,k))
           E0d_yl = (E0(i,j,k) - E0(i,j-1,k))
           Eavg   = (E(i,j,k)  + E(i,j-1,k) )*0.5d0
           E0avg  = (E0(i,j,k) + E0(i,j-1,k))*0.5d0

           !    compute R for limiters
           R  = max(dyi *abs(Ed_yl )/Eavg,  Rmin)
           R0 = max(dyi0*abs(E0d_yl)/E0avg, Rmin)
           R  = lag*R0 + (1.d0-lag)*R

           !    compute limiter
           if (LTyp == 1) then       ! rational approx. to LP lim. (LP, 1981)
              D_yl  = c*(2.d0*kap_yl+R)/(6.d0*kap_yl*kap_yl+3.d0*kap_yl*R+R*R)
              D0_yl = c*(2.d0*kap0_yl+R0)/(6.d0*kap0_yl*kap0_yl+3.d0*kap0_yl*R0+R0*R0)
           else if (LTyp == 2) then  ! Reynolds approx to LP lim.
              D_yl  = 2.d0*c/pi*datan(R*pi/6.d0/kap_yl)/R
              D0_yl = 2.d0*c/pi*datan(R0*pi/6.d0/kap0_yl)/R0
           else if (LTyp == 3) then  ! no limiter
              D_yl  = c/kap_yl/3.d0
              D0_yl = c/kap0_yl/3.d0
           else if (LTyp == 4) then  ! Zeus limiter
              D_yl  = c*(2.d0*kap_yl+R)/(6.d0*kap_yl*kap_yl+3.d0*kap_yl*R+R*R)
              D0_yl = c*(2.d0*kap0_yl+R0)/(6.d0*kap0_yl*kap0_yl+3.d0*kap0_yl*R0+R0*R0)
           else                       ! standard Levermore-Pomraning (LP, 1981)
              D_yl  = c*(cosh(R/kap_yl)/sinh(R/kap_yl)-kap_yl/R)/R
              D0_yl = c*(cosh(R0/kap0_yl)/sinh(R0/kap0_yl)-kap0_yl/R0)/R0
           endif

           !--------------
           ! x-directional limiter, lower face
           !    compute gradients, face values
           Ed_xl  = (E(i,j,k)  - E(i-1,j,k))
           E0d_xl = (E0(i,j,k) - E0(i-1,j,k))
           Eavg   = (E(i,j,k)  + E(i-1,j,k) )*0.5d0
           E0avg  = (E0(i,j,k) + E0(i-1,j,k))*0.5d0

           !    compute R for limiters
           R  = max(dxi *abs(Ed_xl )/Eavg,  Rmin)
           R0 = max(dxi0*abs(E0d_xl)/E0avg, Rmin)
           R  = lag*R0 + (1.d0-lag)*R

           !    compute limiter
           if (LTyp == 1) then       ! rational approx. to LP lim. (LP, 1981)
              D_xl  = c*(2.d0*kap_xl+R)/(6.d0*kap_xl*kap_xl+3.d0*kap_xl*R+R*R)
              D0_xl = c*(2.d0*kap0_xl+R0)/(6.d0*kap0_xl*kap0_xl+3.d0*kap0_xl*R0+R0*R0)
           else if (LTyp == 2) then  ! Reynolds approx to LP lim.
              D_xl  = 2.d0*c/pi*datan(R*pi/6.d0/kap_xl)/R
              D0_xl = 2.d0*c/pi*datan(R0*pi/6.d0/kap0_xl)/R0
           else if (LTyp == 3) then  ! no limiter
              D_xl  = c/kap_xl/3.d0
              D0_xl = c/kap0_xl/3.d0
           else if (LTyp == 4) then  ! Zeus limiter
              D_xl  = c*(2.d0*kap_xl+R)/(6.d0*kap_xl*kap_xl+3.d0*kap_xl*R+R*R)
              D0_xl = c*(2.d0*kap0_xl+R0)/(6.d0*kap0_xl*kap0_xl+3.d0*kap0_xl*R0+R0*R0)
           else                       ! standard Levermore-Pomraning (LP, 1981)
              D_xl  = c*(cosh(R/kap_xl)/sinh(R/kap_xl)-kap_xl/R)/R
              D0_xl = c*(cosh(R0/kap0_xl)/sinh(R0/kap0_xl)-kap0_xl/R0)/R0
           endif

           !--------------
           ! x-directional limiter, upper face
           !    compute gradients, face values
           Ed_xr  = (E(i,j,k)  - E(i+1,j,k))
           E0d_xr = (E0(i,j,k) - E0(i+1,j,k))
           Eavg   = (E(i,j,k)  + E(i+1,j,k) )*0.5d0
           E0avg  = (E0(i,j,k) + E0(i+1,j,k))*0.5d0

           !    compute R for limiters
           R  = max(dxi *abs(Ed_xr )/Eavg,  Rmin)
           R0 = max(dxi0*abs(E0d_xr)/E0avg, Rmin)
           R  = lag*R0 + (1.d0-lag)*R

           !    compute limiter
           if (LTyp == 1) then       ! rational approx. to LP lim. (LP, 1981)
              D_xr  = c*(2.d0*kap_xr+R)/(6.d0*kap_xr*kap_xr+3.d0*kap_xr*R+R*R)
              D0_xr = c*(2.d0*kap0_xr+R0)/(6.d0*kap0_xr*kap0_xr+3.d0*kap0_xr*R0+R0*R0)
           else if (LTyp == 2) then  ! Reynolds approx to LP lim.
              D_xr  = 2.d0*c/pi*datan(R*pi/6.d0/kap_xr)/R
              D0_xr = 2.d0*c/pi*datan(R0*pi/6.d0/kap0_xr)/R0
           else if (LTyp == 3) then  ! no limiter
              D_xr  = c/kap_xr/3.d0
              D0_xr = c/kap0_xr/3.d0
           else if (LTyp == 4) then  ! Zeus limiter
              D_xr  = c*(2.d0*kap_xr+R)/(6.d0*kap_xr*kap_xr+3.d0*kap_xr*R+R*R)
              D0_xr = c*(2.d0*kap0_xr+R0)/(6.d0*kap0_xr*kap0_xr+3.d0*kap0_xr*R0+R0*R0)
           else                       ! standard Levermore-Pomraning (LP, 1981)
              D_xr  = c*(cosh(R/kap_xr)/sinh(R/kap_xr)-kap_xr/R)/R
              D0_xr = c*(cosh(R0/kap0_xr)/sinh(R0/kap0_xr)-kap0_xr/R0)/R0
           endif

           !--------------
           ! y-directional limiter, upper face
           !    compute gradients, face values
           Ed_yr  = (E(i,j,k)  - E(i,j+1,k))
           E0d_yr = (E0(i,j,k) - E0(i,j+1,k))
           Eavg   = (E(i,j,k)  + E(i,j+1,k) )*0.5d0
           E0avg  = (E0(i,j,k) + E0(i,j+1,k))*0.5d0

           !    compute R for limiters
           R  = max(dyi *abs(Ed_yr )/Eavg,  Rmin)
           R0 = max(dyi0*abs(E0d_yr)/E0avg, Rmin)
           R  = lag*R0 + (1.d0-lag)*R

           !    compute limiter
           if (LTyp == 1) then       ! rational approx. to LP lim. (LP, 1981)
              D_yr  = c*(2.d0*kap_yr+R)/(6.d0*kap_yr*kap_yr+3.d0*kap_yr*R+R*R)
              D0_yr = c*(2.d0*kap0_yr+R0)/(6.d0*kap0_yr*kap0_yr+3.d0*kap0_yr*R0+R0*R0)
           else if (LTyp == 2) then  ! Reynolds approx to LP lim.
              D_yr  = 2.d0*c/pi*datan(R*pi/6.d0/kap_yr)/R
              D0_yr = 2.d0*c/pi*datan(R0*pi/6.d0/kap0_yr)/R0
           else if (LTyp == 3) then  ! no limiter
              D_yr  = c/kap_yr/3.d0
              D0_yr = c/kap0_yr/3.d0
           else if (LTyp == 4) then  ! Zeus limiter
              D_yr  = c*(2.d0*kap_yr+R)/(6.d0*kap_yr*kap_yr+3.d0*kap_yr*R+R*R)
              D0_yr = c*(2.d0*kap0_yr+R0)/(6.d0*kap0_yr*kap0_yr+3.d0*kap0_yr*R0+R0*R0)
           else                       ! standard Levermore-Pomraning (LP, 1981)
              D_yr  = c*(cosh(R/kap_yr)/sinh(R/kap_yr)-kap_yr/R)/R
              D0_yr = c*(cosh(R0/kap0_yr)/sinh(R0/kap0_yr)-kap0_yr/R0)/R0
           endif

           !--------------
           ! z-directional limiter, upper face
           !    compute gradients, face values
           Ed_zr  = (E(i,j,k)  - E(i,j,k+1))
           E0d_zr = (E0(i,j,k) - E0(i,j,k+1))
           Eavg   = (E(i,j,k)  + E(i,j,k+1) )*0.5d0
           E0avg  = (E0(i,j,k) + E0(i,j,k+1))*0.5d0

           !    compute R for limiters
           R  = max(dzi *abs(Ed_zr )/Eavg,  Rmin)
           R0 = max(dzi0*abs(E0d_zr)/E0avg, Rmin)
           R  = lag*R0 + (1.d0-lag)*R

           !    compute limiter
           if (LTyp == 1) then       ! rational approx. to LP lim. (LP, 1981)
              D_zr  = c*(2.d0*kap_zr+R)/(6.d0*kap_zr*kap_zr+3.d0*kap_zr*R+R*R)
              D0_zr = c*(2.d0*kap0_zr+R0)/(6.d0*kap0_zr*kap0_zr+3.d0*kap0_zr*R0+R0*R0)
           else if (LTyp == 2) then  ! Reynolds approx to LP lim.
              D_zr  = 2.d0*c/pi*datan(R*pi/6.d0/kap_zr)/R
              D0_zr = 2.d0*c/pi*datan(R0*pi/6.d0/kap0_zr)/R0
           else if (LTyp == 3) then  ! no limiter
              D_zr  = c/kap_zr/3.d0
              D0_zr = c/kap0_zr/3.d0
           else if (LTyp == 4) then  ! Zeus limiter
              D_zr  = c*(2.d0*kap_zr+R)/(6.d0*kap_zr*kap_zr+3.d0*kap_zr*R+R*R)
              D0_zr = c*(2.d0*kap0_zr+R0)/(6.d0*kap0_zr*kap0_zr+3.d0*kap0_zr*R0+R0*R0)
           else                       ! standard Levermore-Pomraning (LP, 1981)
              D_zr  = c*(cosh(R/kap_zr)/sinh(R/kap_zr)-kap_zr/R)/R
              D0_zr = c*(cosh(R0/kap0_zr)/sinh(R0/kap0_zr)-kap0_zr/R0)/R0
           endif

           ! put it all together
           !    Edot - Div(D(E)*Grad(E)) - 4*pi*src + c*kap*E = 0,
           res(i,j,k) = -src(i,j,k)*etafac                        &
                + (1.d0 + dtfac*4.d0*kap)*E(i,j,k)                &
                  - dtfac*dxi*dxi*(D_xr*Ed_xr-D_xl*Ed_xl)         &
                  - dtfac*dyi*dyi*(D_yr*Ed_yr-D_yl*Ed_yl)         &
                  - dtfac*dzi*dzi*(D_zr*Ed_zr-D_zr*Ed_zl)         &
                - (1.d0 - dtfac0*4.d0*kap0)*E0(i,j,k)             &
                  - dtfac0*dxi0*dxi0*(D0_xr*E0d_xr-D0_xl*E0d_xl)  &
                  - dtfac0*dyi0*dyi0*(D0_yr*E0d_yr-D0_yl*E0d_yl)  &
                  - dtfac0*dzi0*dzi0*(D0_zr*E0d_zr-D0_zr*E0d_zl)

        enddo
     enddo
  enddo

  return
end subroutine MFProb_RadResid_3D
!=======================================================================






subroutine MFProb_RadResid_2D(res, E, E0, HI, HI0, HeI, HeI0, HeII,    &
     HeII0, src, LTyp, LImp, dt, theta, sHI, sHeI, sHeII, a, a0, adot, &
     adot0, aUn, lUn, lUn0, rUn, rUn0, nUn, nUn0, dx, dy, Nx, Ny,      &
     NGxl, NGxr, NGyl, NGyr, ier)
!=======================================================================
!  PURPOSE: 2D version of the routine
!=======================================================================
#include "fortran.def"
  implicit none

  !--------------
  ! argument declarations
  integer, intent(in)  :: LTyp, LImp
  integer, intent(in)  :: Nx, NGxl, NGxr
  integer, intent(in)  :: Ny, NGyl, NGyr
  integer, intent(out) :: ier
  REALSUB, intent(in)  :: a, adot, a0, adot0
  real, intent(in) :: aUn, lUn, lUn0, rUn, rUn0, nUn, nUn0, sHI, sHeI, sHeII
  real, intent(in) :: dt, theta, dx, dy
  real, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr), intent(in) :: &
       E, E0, src
  real, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr), intent(in) :: &
       HI, HeI, HeII, HI0, HeI0, HeII0
  real, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr), intent(out) :: res

  !--------------
  ! locals
  integer  :: i, j
  real :: dtfac, dtfac0, etafac, c, pi
  real :: dxi, dxi0, dyi, dyi0
  real :: D_xl, D0_xl, D_xr, D0_xr, Ed_xl, E0d_xl, Ed_xr, E0d_xr
  real :: D_yl, D0_yl, D_yr, D0_yr, Ed_yl, E0d_yl, Ed_yr, E0d_yr
  real :: Eavg, E0avg, R, R0, kap, kap0, lag, Rmin
  real :: kap_xl, kap0_xl, kap_xr, kap0_xr
  real :: kap_yl, kap0_yl, kap_yr, kap0_yr

  !=======================================================================

  ! initialize output to zero, flag to success
  res = 0.d0
  ier = 1

  ! set shortcut values 
  c = 2.99792458d10             ! speed of light [cm/s]
  pi = 4.d0*atan(1.d0)
  Rmin = 1.0d-20
  dtfac  = dt*theta
  dtfac0 = dt*(1.d0-theta)
  etafac = dt*8.d0*pi/(rUn+rUn0)
  dxi    = a/dx/lUn
  dyi    = a/dy/lUn
  dxi0   = a0/dx/lUn0
  dyi0   = a0/dy/lUn0

  ! set lag constant to account for lagging E in limiter evaluation
  lag = 0.d0
  if (LImp == 0)  lag = 1.d0

  ! compute radiation energy gradient over domain
  do j=1,Ny
     do i=1,Nx
        
        ! compute old opacities at faces and cell center
        kap0_yl = (sHI*(  HI0(i,j)  +HI0(i,j-1))  &
                 + sHeI*( HeI0(i,j) +HeI0(i,j-1)) &
                 + sHeII*(HeII0(i,j)+HeII0(i,j-1)))*0.5d0*nUn0
        kap0_xl = (sHI*(  HI0(i,j)  +HI0(i-1,j))  &
                 + sHeI*( HeI0(i,j) +HeI0(i-1,j)) &
                 + sHeII*(HeII0(i,j)+HeII0(i-1,j)))*0.5d0*nUn0
        kap0_xr = (sHI*(  HI0(i,j)  +HI0(i+1,j))  &
                 + sHeI*( HeI0(i,j) +HeI0(i+1,j)) &
                 + sHeII*(HeII0(i,j)+HeII0(i+1,j)))*0.5d0*nUn0
        kap0_yr = (sHI*(  HI0(i,j)  +HI0(i,j+1))  &
                 + sHeI*( HeI0(i,j) +HeI0(i,j+1)) &
                 + sHeII*(HeII0(i,j)+HeII0(i,j+1)))*0.5d0*nUn0
        kap0 = (sHI*HI0(i,j) + sHeI*HeI0(i,j) + sHeII*HeII0(i,j))*nUn0

        ! compute new opacities at faces and cell center
        if (lag == 1) then
           kap_yl = kap0_yl
           kap_xl = kap0_xl
           kap_xr = kap0_xr
           kap_yr = kap0_yr
        else
           kap_yl = (sHI*(  HI(i,j)  +HI(i,j-1))  &
                   + sHeI*( HeI(i,j) +HeI(i,j-1)) &
                   + sHeII*(HeII(i,j)+HeII(i,j-1)))*0.5d0*nUn
           kap_xl = (sHI*(  HI(i,j)  +HI(i-1,j))  &
                   + sHeI*( HeI(i,j) +HeI(i-1,j)) &
                   + sHeII*(HeII(i,j)+HeII(i-1,j)))*0.5d0*nUn
           kap_xr = (sHI*(  HI(i,j)  +HI(i+1,j))  &
                   + sHeI*( HeI(i,j) +HeI(i+1,j)) &
                   + sHeII*(HeII(i,j)+HeII(i+1,j)))*0.5d0*nUn
           kap_yr = (sHI*(  HI(i,j)  +HI(i,j+1))  &
                   + sHeI*( HeI(i,j) +HeI(i,j+1)) &
                   + sHeII*(HeII(i,j)+HeII(i,j+1)))*0.5d0*nUn
        endif
        kap = (sHI*HI(i,j) + sHeI*HeI(i,j) + sHeII*HeII(i,j))*nUn


        !--------------
        ! y-directional limiter, lower face
        !    compute gradients, face values
        Ed_yl  = (E(i,j)  - E(i,j-1))
        E0d_yl = (E0(i,j) - E0(i,j-1))
        Eavg   = (E(i,j)  + E(i,j-1) )*0.5d0
        E0avg  = (E0(i,j) + E0(i,j-1))*0.5d0
        
        !    compute R for limiters
        R  = max(dyi *abs(Ed_yl )/Eavg,  Rmin)
        R0 = max(dyi0*abs(E0d_yl)/E0avg, Rmin)
        R  = lag*R0 + (1.d0-lag)*R
        
        !    compute limiter
        if (LTyp == 1) then       ! rational approx. to LP lim. (LP, 1981)
           D_yl  = c*(2.d0*kap_yl+R)/(6.d0*kap_yl*kap_yl+3.d0*kap_yl*R+R*R)
           D0_yl = c*(2.d0*kap0_yl+R0)/(6.d0*kap0_yl*kap0_yl+3.d0*kap0_yl*R0+R0*R0)
        else if (LTyp == 2) then  ! Reynolds approx to LP lim.
           D_yl  = 2.d0*c/pi*datan(R*pi/6.d0/kap_yl)/R
           D0_yl = 2.d0*c/pi*datan(R0*pi/6.d0/kap0_yl)/R0
        else if (LTyp == 3) then  ! no limiter
           D_yl  = c/kap_yl/3.d0
           D0_yl = c/kap0_yl/3.d0
        else if (LTyp == 4) then  ! Zeus limiter
           D_yl  = c*(2.d0*kap_yl+R)/(6.d0*kap_yl*kap_yl+3.d0*kap_yl*R+R*R)
           D0_yl = c*(2.d0*kap0_yl+R0)/(6.d0*kap0_yl*kap0_yl+3.d0*kap0_yl*R0+R0*R0)
        else                       ! standard Levermore-Pomraning (LP, 1981)
           D_yl  = c*(cosh(R/kap_yl)/sinh(R/kap_yl)-kap_yl/R)/R
           D0_yl = c*(cosh(R0/kap0_yl)/sinh(R0/kap0_yl)-kap0_yl/R0)/R0
        endif
        
        !--------------
        ! x-directional limiter, lower face
        !    compute gradients, face values
        Ed_xl  = (E(i,j)  - E(i-1,j))
        E0d_xl = (E0(i,j) - E0(i-1,j))
        Eavg   = (E(i,j)  + E(i-1,j) )*0.5d0
        E0avg  = (E0(i,j) + E0(i-1,j))*0.5d0
        
        !    compute R for limiters
        R  = max(dxi *abs(Ed_xl )/Eavg,  Rmin)
        R0 = max(dxi0*abs(E0d_xl)/E0avg, Rmin)
        R  = lag*R0 + (1.d0-lag)*R
        
        !    compute limiter
        if (LTyp == 1) then       ! rational approx. to LP lim. (LP, 1981)
           D_xl  = c*(2.d0*kap_xl+R)/(6.d0*kap_xl*kap_xl+3.d0*kap_xl*R+R*R)
           D0_xl = c*(2.d0*kap0_xl+R0)/(6.d0*kap0_xl*kap0_xl+3.d0*kap0_xl*R0+R0*R0)
        else if (LTyp == 2) then  ! Reynolds approx to LP lim.
           D_xl  = 2.d0*c/pi*datan(R*pi/6.d0/kap_xl)/R
           D0_xl = 2.d0*c/pi*datan(R0*pi/6.d0/kap0_xl)/R0
        else if (LTyp == 3) then  ! no limiter
           D_xl  = c/kap_xl/3.d0
           D0_xl = c/kap0_xl/3.d0
        else if (LTyp == 4) then  ! Zeus limiter
           D_xl  = c*(2.d0*kap_xl+R)/(6.d0*kap_xl*kap_xl+3.d0*kap_xl*R+R*R)
           D0_xl = c*(2.d0*kap0_xl+R0)/(6.d0*kap0_xl*kap0_xl+3.d0*kap0_xl*R0+R0*R0)
        else                       ! standard Levermore-Pomraning (LP, 1981)
           D_xl  = c*(cosh(R/kap_xl)/sinh(R/kap_xl)-kap_xl/R)/R
           D0_xl = c*(cosh(R0/kap0_xl)/sinh(R0/kap0_xl)-kap0_xl/R0)/R0
        endif
        
        !--------------
        ! x-directional limiter, upper face
        !    compute gradients, face values
        Ed_xr  = (E(i,j)  - E(i+1,j))
        E0d_xr = (E0(i,j) - E0(i+1,j))
        Eavg   = (E(i,j)  + E(i+1,j) )*0.5d0
        E0avg  = (E0(i,j) + E0(i+1,j))*0.5d0
        
        !    compute R for limiters
        R  = max(dxi *abs(Ed_xr )/Eavg,  Rmin)
        R0 = max(dxi0*abs(E0d_xr)/E0avg, Rmin)
        R  = lag*R0 + (1.d0-lag)*R
        
        !    compute limiter
        if (LTyp == 1) then       ! rational approx. to LP lim. (LP, 1981)
           D_xr  = c*(2.d0*kap_xr+R)/(6.d0*kap_xr*kap_xr+3.d0*kap_xr*R+R*R)
           D0_xr = c*(2.d0*kap0_xr+R0)/(6.d0*kap0_xr*kap0_xr+3.d0*kap0_xr*R0+R0*R0)
        else if (LTyp == 2) then  ! Reynolds approx to LP lim.
           D_xr  = 2.d0*c/pi*datan(R*pi/6.d0/kap_xr)/R
           D0_xr = 2.d0*c/pi*datan(R0*pi/6.d0/kap0_xr)/R0
        else if (LTyp == 3) then  ! no limiter
           D_xr  = c/kap_xr/3.d0
           D0_xr = c/kap0_xr/3.d0
        else if (LTyp == 4) then  ! Zeus limiter
           D_xr  = c*(2.d0*kap_xr+R)/(6.d0*kap_xr*kap_xr+3.d0*kap_xr*R+R*R)
           D0_xr = c*(2.d0*kap0_xr+R0)/(6.d0*kap0_xr*kap0_xr+3.d0*kap0_xr*R0+R0*R0)
        else                       ! standard Levermore-Pomraning (LP, 1981)
           D_xr  = c*(cosh(R/kap_xr)/sinh(R/kap_xr)-kap_xr/R)/R
           D0_xr = c*(cosh(R0/kap0_xr)/sinh(R0/kap0_xr)-kap0_xr/R0)/R0
        endif
        
        !--------------
        ! y-directional limiter, upper face
        !    compute gradients, face values
        Ed_yr  = (E(i,j)  - E(i,j+1))
        E0d_yr = (E0(i,j) - E0(i,j+1))
        Eavg   = (E(i,j)  + E(i,j+1) )*0.5d0
        E0avg  = (E0(i,j) + E0(i,j+1))*0.5d0
        
        !    compute R for limiters
        R  = max(dyi *abs(Ed_yr )/Eavg,  Rmin)
        R0 = max(dyi0*abs(E0d_yr)/E0avg, Rmin)
        R  = lag*R0 + (1.d0-lag)*R
        
        !    compute limiter
        if (LTyp == 1) then       ! rational approx. to LP lim. (LP, 1981)
           D_yr  = c*(2.d0*kap_yr+R)/(6.d0*kap_yr*kap_yr+3.d0*kap_yr*R+R*R)
           D0_yr = c*(2.d0*kap0_yr+R0)/(6.d0*kap0_yr*kap0_yr+3.d0*kap0_yr*R0+R0*R0)
        else if (LTyp == 2) then  ! Reynolds approx to LP lim.
           D_yr  = 2.d0*c/pi*datan(R*pi/6.d0/kap_yr)/R
           D0_yr = 2.d0*c/pi*datan(R0*pi/6.d0/kap0_yr)/R0
        else if (LTyp == 3) then  ! no limiter
           D_yr  = c/kap_yr/3.d0
           D0_yr = c/kap0_yr/3.d0
        else if (LTyp == 4) then  ! Zeus limiter
           D_yr  = c*(2.d0*kap_yr+R)/(6.d0*kap_yr*kap_yr+3.d0*kap_yr*R+R*R)
           D0_yr = c*(2.d0*kap0_yr+R0)/(6.d0*kap0_yr*kap0_yr+3.d0*kap0_yr*R0+R0*R0)
        else                       ! standard Levermore-Pomraning (LP, 1981)
           D_yr  = c*(cosh(R/kap_yr)/sinh(R/kap_yr)-kap_yr/R)/R
           D0_yr = c*(cosh(R0/kap0_yr)/sinh(R0/kap0_yr)-kap0_yr/R0)/R0
        endif
        
        ! put it all together
        !    Edot - Div(D(E)*Grad(E)) - 4*pi*src + c*kap*E = 0,
        res(i,j) = -src(i,j)*etafac                            &
             + (1.d0 + dtfac*4.d0*kap)*E(i,j)                  &
               - dtfac*dxi*dxi*(D_xr*Ed_xr-D_xl*Ed_xl)         &
               - dtfac*dyi*dyi*(D_yr*Ed_yr-D_yl*Ed_yl)         &
             - (1.d0 - dtfac0*4.d0*kap0)*E0(i,j)               &
               - dtfac0*dxi0*dxi0*(D0_xr*E0d_xr-D0_xl*E0d_xl)  &
               - dtfac0*dyi0*dyi0*(D0_yr*E0d_yr-D0_yl*E0d_yl)

     enddo
  enddo

  return
end subroutine MFProb_RadResid_2D
!=======================================================================






subroutine MFProb_RadResid_1D(res, E, E0, HI, HI0, HeI, HeI0, HeII,    &
     HeII0, src, LTyp, LImp, dt, theta, sHI, sHeI, sHeII, a, a0, adot, &
     adot0, aUn, lUn, lUn0, rUn, rUn0, nUn, nUn0, dx, Nx, NGxl, NGxr, ier)
!=======================================================================
!  PURPOSE: 1D version of the routine
!=======================================================================
#include "fortran.def"
  implicit none

  !--------------
  ! argument declarations
  integer, intent(in)  :: LTyp, LImp
  integer, intent(in)  :: Nx, NGxl, NGxr
  integer, intent(out) :: ier
  REALSUB, intent(in)  :: a, adot, a0, adot0
  real, intent(in) :: aUn, lUn, lUn0, rUn, rUn0, nUn, nUn0, sHI, sHeI, sHeII
  real, intent(in) :: dt, theta, dx
  real, dimension(1-NGxl:Nx+NGxr), intent(in) :: E, E0, src
  real, dimension(1-NGxl:Nx+NGxr), intent(in) :: &
       HI, HeI, HeII, HI0, HeI0, HeII0
  real, dimension(1-NGxl:Nx+NGxr), intent(out) :: res

  !--------------
  ! locals
  integer  :: i
  real :: dtfac, dtfac0, etafac, c, pi
  real :: dxi, dxi0
  real :: D_xl, D0_xl, D_xr, D0_xr, Ed_xl, E0d_xl, Ed_xr, E0d_xr
  real :: Eavg, E0avg, R, R0, kap, kap0, lag, Rmin
  real :: kap_xl, kap0_xl, kap_xr, kap0_xr

  !=======================================================================

  ! initialize output to zero, flag to success
  res = 0.d0
  ier = 1

  ! set shortcut values 
  c = 2.99792458d10             ! speed of light [cm/s]
  pi = 4.d0*atan(1.d0)
  Rmin = 1.0d-20
  dtfac  = dt*theta
  dtfac0 = dt*(1.d0-theta)
  etafac = dt*8.d0*pi/(rUn+rUn0)
  dxi    = a/dx/lUn
  dxi0   = a0/dx/lUn0

  ! set lag constant to account for lagging E in limiter evaluation
  lag = 0.d0
  if (LImp == 0)  lag = 1.d0

  ! compute radiation energy gradient over domain
  do i=1,Nx

     ! compute old opacities at faces and cell center
     kap0_xl = (sHI*(  HI0(i)  +HI0(i-1))  &
              + sHeI*( HeI0(i) +HeI0(i-1)) &
              + sHeII*(HeII0(i)+HeII0(i-1)))*0.5d0*nUn0
     kap0_xr = (sHI*(  HI0(i)  +HI0(i+1))  &
              + sHeI*( HeI0(i) +HeI0(i+1)) &
              + sHeII*(HeII0(i)+HeII0(i+1)))*0.5d0*nUn0
     kap0 = (sHI*HI0(i) + sHeI*HeI0(i) + sHeII*HeII0(i))*nUn0

     ! compute new opacities at faces and cell center
     if (lag == 1) then
        kap_xl = kap0_xl
        kap_xr = kap0_xr
     else
        kap_xl = (sHI*(  HI(i)  +HI(i-1))  &
                + sHeI*( HeI(i) +HeI(i-1)) &
                + sHeII*(HeII(i)+HeII(i-1)))*0.5d0*nUn
        kap_xr = (sHI*(  HI(i)  +HI(i+1))  &
                + sHeI*( HeI(i) +HeI(i+1)) &
                + sHeII*(HeII(i)+HeII(i+1)))*0.5d0*nUn
     endif
     kap = (sHI*HI(i) + sHeI*HeI(i) + sHeII*HeII(i))*nUn

     
     !--------------
     ! x-directional limiter, lower face
     !    compute gradients, face values
     Ed_xl  = (E(i)  - E(i-1))
     E0d_xl = (E0(i) - E0(i-1))
     Eavg   = (E(i)  + E(i-1) )*0.5d0
     E0avg  = (E0(i) + E0(i-1))*0.5d0
     
     !    compute R for limiters
     R  = max(dxi *abs(Ed_xl )/Eavg,  Rmin)
     R0 = max(dxi0*abs(E0d_xl)/E0avg, Rmin)
     R  = lag*R0 + (1.d0-lag)*R
     
     !    compute limiter
     if (LTyp == 1) then       ! rational approx. to LP lim. (LP, 1981)
        D_xl  = c*(2.d0*kap_xl+R)/(6.d0*kap_xl*kap_xl+3.d0*kap_xl*R+R*R)
        D0_xl = c*(2.d0*kap0_xl+R0)/(6.d0*kap0_xl*kap0_xl+3.d0*kap0_xl*R0+R0*R0)
     else if (LTyp == 2) then  ! Reynolds approx to LP lim.
        D_xl  = 2.d0*c/pi*datan(R*pi/6.d0/kap_xl)/R
        D0_xl = 2.d0*c/pi*datan(R0*pi/6.d0/kap0_xl)/R0
     else if (LTyp == 3) then  ! no limiter
        D_xl  = c/kap_xl/3.d0
        D0_xl = c/kap0_xl/3.d0
     else if (LTyp == 4) then  ! Zeus limiter
        D_xl  = c*(2.d0*kap_xl+R)/(6.d0*kap_xl*kap_xl+3.d0*kap_xl*R+R*R)
        D0_xl = c*(2.d0*kap0_xl+R0)/(6.d0*kap0_xl*kap0_xl+3.d0*kap0_xl*R0+R0*R0)
     else                       ! standard Levermore-Pomraning (LP, 1981)
        D_xl  = c*(cosh(R/kap_xl)/sinh(R/kap_xl)-kap_xl/R)/R
        D0_xl = c*(cosh(R0/kap0_xl)/sinh(R0/kap0_xl)-kap0_xl/R0)/R0
     endif
     
     !--------------
     ! x-directional limiter, upper face
     !    compute gradients, face values
     Ed_xr  = (E(i)  - E(i+1))
     E0d_xr = (E0(i) - E0(i+1))
     Eavg   = (E(i)  + E(i+1) )*0.5d0
     E0avg  = (E0(i) + E0(i+1))*0.5d0
     
     !    compute R for limiters
     R  = max(dxi *abs(Ed_xr )/Eavg,  Rmin)
     R0 = max(dxi0*abs(E0d_xr)/E0avg, Rmin)
     R  = lag*R0 + (1.d0-lag)*R
     
     !    compute limiter
     if (LTyp == 1) then       ! rational approx. to LP lim. (LP, 1981)
        D_xr  = c*(2.d0*kap_xr+R)/(6.d0*kap_xr*kap_xr+3.d0*kap_xr*R+R*R)
        D0_xr = c*(2.d0*kap0_xr+R0)/(6.d0*kap0_xr*kap0_xr+3.d0*kap0_xr*R0+R0*R0)
     else if (LTyp == 2) then  ! Reynolds approx to LP lim.
        D_xr  = 2.d0*c/pi*datan(R*pi/6.d0/kap_xr)/R
        D0_xr = 2.d0*c/pi*datan(R0*pi/6.d0/kap0_xr)/R0
     else if (LTyp == 3) then  ! no limiter
        D_xr  = c/kap_xr/3.d0
        D0_xr = c/kap0_xr/3.d0
     else if (LTyp == 4) then  ! Zeus limiter
        D_xr  = c*(2.d0*kap_xr+R)/(6.d0*kap_xr*kap_xr+3.d0*kap_xr*R+R*R)
        D0_xr = c*(2.d0*kap0_xr+R0)/(6.d0*kap0_xr*kap0_xr+3.d0*kap0_xr*R0+R0*R0)
     else                       ! standard Levermore-Pomraning (LP, 1981)
        D_xr  = c*(cosh(R/kap_xr)/sinh(R/kap_xr)-kap_xr/R)/R
        D0_xr = c*(cosh(R0/kap0_xr)/sinh(R0/kap0_xr)-kap0_xr/R0)/R0
     endif
     
     ! put it all together
     !    Edot - Div(D(E)*Grad(E)) - 4*pi*src + c*kap*E = 0,
     res(i) = -src(i)*etafac                             &
          + (1.d0 + dtfac*4.d0*kap)*E(i)                 &
            - dtfac*dxi*dxi*(D_xr*Ed_xr-D_xl*Ed_xl)      &
          - (1.d0 - dtfac0*4.d0*kap0)*E0(i)              &
            - dtfac0*dxi0*dxi0*(D0_xr*E0d_xr-D0_xl*E0d_xl)

  enddo
  
  return
end subroutine MFProb_RadResid_1D
!=======================================================================
