!=======================================================================
!
! Copyright 2009 Daniel R. Reynolds
!
! This software is released under the terms of the "Enzo Public License"
! in the accompanying LICENSE file.
!
!=======================================================================
subroutine MFProb_ComputeRadiationIntegrals(piHI, piHeI, piHeII, GHI, &
     GHeI, GHeII, Ef, E1, E2, E3, E1old, E2old, E3old, Nchem,         &
     ESpectrum, fsUn, rUn, rUn0, Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr,  &
     NGzl, NGzr, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       August 2009
!
!  PURPOSE: Computes numerical quadratures (frequency space) for 
!           radiation integrals to obtain photo-ionization and 
!           photo-heating terms.
!
!           These are computed using a simple 4th-order accurate 
!           numerical quadrature rule (composite Simpson's), on the 
!           definite integrals
!                  int_{nu0}^{nu1} f(nu) dnu
!                  int_{nu1}^{nu2} f(nu) dnu
!           and the re-mapped indefinite integral
!                  int_{nu2}^{inf} f(nu) dnu
!                = int_0^1 nu2/(x^2)*f(nu2/x) dx
!           The composite Simpson's rule computes each integral as 
!           the sum of small quadrature intervals [xl,xr]_i, where 
!                  [0,1] = Union_i [xl,xr]_i,
!           and xm = (xl+xr)/2, via the quadrature formula
!                  int_{xl}^{xr} g(x) dx 
!                = (xr-xl)/6*(g(xl) + 4*g(xm) + g(xr))
!           Note: since 1/0 = infty, we begin the mapped quadrature 
!           integration at a small positive value.
!
!  INPUTS:
!     Ef         - free-streaming radiation (new time)
!     E1         - monochromatic radiation at HI ionization threshold
!     E2         - monochromatic radiation at HeI ionization threshold
!     E3         - monochromatic radiation at HeII ionization threshold
!     E1old      - E1 at old time step
!     E2old      - E2 at old time step
!     E3old      - E3 at old time step
!     Nchem      - number of chemical species
!     ESpectrum  - radiation spectrum type
!     fsUn       - free-streaming radiation unit scaling
!     rUn,rUn0   - radiation unit scaling (new and old time steps)
!     Nx,Ny,Nz   - active mesh size in each direction
!     NG*l/NG*r  - left/right ghost cells in each direction
!
!     Note: the vector inputs are of size (Nx + NGxl + NGxr) in 
!           the x-direction, others are similar.
!     Note: the photo-heating coefficients must still be multiplied by 
!           n_i/rho and summed up (over i) for the photo-heating rate G.
!
!  OUTPUT ARGUMENTS: 
!     piHI       - array of HI   photo-ionization coefficients
!     piHeI      - array of HeI  photo-ionization coefficients
!     piHeII     - array of HeII photo-ionization coefficients
!     GHI        - array of HI   photo-heating coefficients
!     GHeI       - array of HeI  photo-heating coefficients
!     GHeII      - array of HeII photo-heating coefficients
!     ier        - success/failure flag (1->success, 0->failure)
!
!=======================================================================
  implicit none
#include "fortran.def"

  !--------------
  ! argument declarations
  integer, intent(in) :: Nchem, ESpectrum
  integer, intent(in) :: Nx, NGxl, NGxr
  integer, intent(in) :: Ny, NGyl, NGyr
  integer, intent(in) :: Nz, NGzl, NGzr
  integer, intent(out) :: ier
  real, intent(in) :: fsUn, rUn, rUn0
  real, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr), &
       intent(in) :: Ef, E1, E2, E3, E1old, E2old, E3old
  real, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr), &
       intent(out) :: piHI, piHeI, piHeII, GHI, GHeI, GHeII
  
  !--------------
  ! locals
  integer :: i, j, k, l, nbins
!  integer, parameter :: nnus=21     ! must be odd
  integer, parameter :: nnus=13     ! must be odd
  real :: hp, ev2erg, c, epsilon, Llimit, Ulimit
  real :: nu0_HI, nu0_HeI, nu0_HeII
  real :: Efval, E1val, E2val, E3val, Eval, chibar
  real :: s1n1, s1n2, s1n3, s2n2, s2n3, s3n3, chin1, chin2, chin3
  real, dimension(nnus) :: nusA, nusB, nusC, etasC
  real, dimension(nnus) :: fHIA, fHIB, fHIC
  real, dimension(nnus) :: fHeIB, fHeIC, fHeIIC
  real, dimension(nnus) :: fnuA, fnuB, fnuC
  real, dimension(nnus) :: chiA, chiB, chiC
  real, dimension(nnus) :: S1A, S1B, S1C
  real, dimension(nnus) :: S2B, S2C, S3C
  real, dimension(nnus) :: fE
  real, dimension(nnus) :: wtsA, wtsB, wtsC

  !=======================================================================

!!$  print *, '  '
!!$  print *, 'Entering MFProb_ComputeRadiationIntegrals'
!!$  print *, '   fsUn =',fsUn
!!$  print *, '   rUn* =',rUn,rUn0

  ! initialize integrals, return flag
  ier = 1
  piHI = 0.d0
  GHI = 0.d0
  if (Nchem > 1) then
     piHeI  = 0.d0
     piHeII = 0.d0
     GHeI   = 0.d0
     GHeII  = 0.d0
  endif

  ! set constants
  hp = 6.6260693d-27            ! Planck's constant (ergs*s)
  ev2erg = 1.60217653d-12       ! conversion constant from eV to ergs
  c = 2.99792458d10             ! speed of light (cm/s)
  nu0_HI   = 13.6d0*ev2erg/hp   ! ionization threshold of HI (hz)
  nu0_HeI  = 24.6d0*ev2erg/hp   ! ionization threshold of HeI (hz)
  nu0_HeII = 54.4d0*ev2erg/hp   ! ionization threshold of HeII (hz)
  epsilon = 1.d0                ! floating point roundoff
  do while ((1.d0 + epsilon) > 1.d0)  
     epsilon = epsilon*0.5d0
  end do
  Llimit = 1.d-1                ! lower limit of int (shift away from 0)
  Ulimit = 1.d0 - sqrt(epsilon) ! upper limit of int (shift away from 1)

  ! set frequencies for evaluation
  nbins = (nnus-1)/2
  do l=1,nnus
     ! nusA = linspace(nu0_HI,nu0_HeI,nnus)
     nusA(l) = nu0_HI + (l-1)*(0.9999999999d0*nu0_HeI-nu0_HI)/(nnus-1)
     ! nusB = linspace(nu0_HeI,nu0_HeII,nnus)
     nusB(l) = nu0_HeI + (l-1)*(0.9999999999d0*nu0_HeII-nu0_HeI)/(nnus-1)
     ! etasC = linspace(Llimit,Ulimit,nnus)
     etasC(l) = Llimit + (l-1)*(Ulimit-Llimit)/(nnus-1)
     ! nusC = nu0_HeII./etasC
     nusC(l) = nu0_HeII/etasC(l)
  enddo

!!$  print '(A,21(es9.1,1x))','  nusA =',nusA
!!$  print '(A,21(es9.1,1x))','  nusB =',nusB
!!$  print '(A,21(es9.1,1x))','  etasC =',etasC
!!$  print '(A,21(es9.1,1x))','  nusC =',nusC

  ! evaluate frequency-dependent functions at nu values
  !    species cross-sections at ionization thresholds
  call MFProb_HICrossSection(s1n1, nu0_HI, 1)
  call MFProb_HICrossSection(s1n2, nu0_HeI,1)
  call MFProb_HICrossSection(s1n3, nu0_HeII,1)
  call MFProb_HeICrossSection(s2n2, nu0_HeI, 1)
  call MFProb_HeICrossSection(s2n3, nu0_HeII, 1)
  call MFProb_HeIICrossSection(s3n3, nu0_HeII, 1)
  s1n2 = s1n2/s1n1
  s1n3 = s1n3/s1n1
  s2n3 = s2n3/s2n2
!!$  print '(A,6(es9.1))', '  s*n* =', s1n1, s1n2, s1n3, s2n2, s2n3, s3n3
  
  !    species cross sections across intervals
  call MFProb_HICrossSection(fHIA, nusA, nnus)
  call MFProb_HICrossSection(fHIB, nusB, nnus)
  call MFProb_HICrossSection(fHIC, nusC, nnus)
  call MFProb_HeICrossSection(fHeIB, nusB, nnus)
  call MFProb_HeICrossSection(fHeIC, nusC, nnus)
  call MFProb_HeIICrossSection(fHeIIC, nusC, nnus)
!!$  print '(A,6(es9.1))', '  fH* =', sum(fHIA), sum(fHIB), sum(fHIC), &
!!$       sum(fHeIB), sum(fHeIC), sum(fHeIIC)
  !    radiation spectrum across intervals
  call MFProb_RadiationSpectrum(chiA, nusA, ESpectrum, nnus)
  call MFProb_RadiationSpectrum(chiB, nusB, ESpectrum, nnus)
  call MFProb_RadiationSpectrum(chiC, nusC, ESpectrum, nnus)
!!$  print '(A,3(es9.1))', '  chi* =', sum(chiA), sum(chiB), sum(chiC)
  !    exponents for radiation terms
  S3C = fHeIIC/s3n3
  S2B = fHeIB/s2n2
  S1A = fHIA/s1n1
  S2C = fHeIC/s2n2 - s2n3*S3C
  S1B = fHIB/s1n1  - s1n2*S2B
  S1C = fHIC/s1n1  - s1n3*S3C - s1n2*S2C
!!$  print '(A,6(es9.1))', '  S* =', sum(S3C), sum(S2B), sum(S1A), &
!!$       sum(S2C), sum(S1B), sum(S1C)

  ! compute the quadrature weights across the intervals
  wtsA = 0.d0
  wtsB = 0.d0
  wtsC = 0.d0
  do l=1,nbins
     i = 2*l-1
     j = 2*l
     k = 2*l+1
     wtsA(i) = wtsA(i) + (nusA(k)-nusA(i))/6.d0
     wtsA(j) = wtsA(j) + (nusA(k)-nusA(i))*2.d0/3.d0
     wtsA(k) = wtsA(k) + (nusA(k)-nusA(i))/6.d0
     wtsB(i) = wtsB(i) + (nusB(k)-nusB(i))/6.d0
     wtsB(j) = wtsB(j) + (nusB(k)-nusB(i))*2.d0/3.d0
     wtsB(k) = wtsB(k) + (nusB(k)-nusB(i))/6.d0
     wtsC(i) = wtsC(i) + (etasC(k)-etasC(i))/6.d0*nu0_HeII/etasC(i)**2
     wtsC(j) = wtsC(j) + (etasC(k)-etasC(i))*2.d0/3.d0*nu0_HeII/etasC(j)**2
     wtsC(k) = wtsC(k) + (etasC(k)-etasC(i))/6.d0*nu0_HeII/etasC(k)**2
  enddo
!!$  print '(A,3(es9.1,1x))', '  wts* =', sum(wtsA), sum(wtsB), sum(wtsC)

  ! compute the integrated radiation spectrum
  chibar = sum(wtsA*chiA) + sum(wtsB*chiB) + sum(wtsC*chiC)
!!$  print '(A,1(es9.1))', '  chibar =', chibar

  ! rescale weights with common factors
  wtsA = wtsA*c
  wtsB = wtsB*c
  wtsC = wtsC*c

  !    radiation spectrum at ionization thresholds
  call MFProb_RadiationSpectrum(chin1, nu0_HI,   ESpectrum, 1)
  call MFProb_RadiationSpectrum(chin2, nu0_HeI,  ESpectrum, 1)
  call MFProb_RadiationSpectrum(chin3, nu0_HeII, ESpectrum, 1)
!!$  print '(A,3(es9.1))', '  chin* =', chin1, chin2, chin3

  !    rescale radiation spectra
  chiA  = chiA/chibar
  chiB  = chiB/chibar
  chiC  = chiC/chibar
  chin1 = chin1/chibar
  chin2 = chin2/chibar
  chin3 = chin3/chibar
!!$  print '(A,3(es9.1))', '  chin* =', chin1, chin2, chin3

  ! single-species problem; compute everything for hydrogen only
  if (Nchem == 1) then

     ! iterate over the domain
     do k=1,Nz
        do j=1,Ny
           do i=1,Nx

              ! get time-centered values for radiation species
              Efval = Ef(i,j,k)*fsUn
              E1val = 0.5d0*(E1(i,j,k)*rUn + E1old(i,j,k)*rUn0)/Efval/chin1
              E2val = 0.5d0*(E2(i,j,k)*rUn + E2old(i,j,k)*rUn0)/Efval/chin2
              E3val = 0.5d0*(E3(i,j,k)*rUn + E3old(i,j,k)*rUn0)/Efval/chin3

              ! initialize the outputs
              piHI(i,j,k) = 0.d0
              GHI(i,j,k)  = 0.d0

              ! get E values at prescribed frequencies
              fE = wtsA * chiA * E1val**S1A
                 
              ! compute the integrated photo-ionization and photo-heating terms
              do l=1,nnus
                 piHI(i,j,k) = piHI(i,j,k) + fHIA(l)*fE(l)/nusA(l)
                 GHI(i,j,k)  = GHI(i,j,k)  + fHIA(l)*fE(l)*(1.d0-nu0_HI/nusA(l))
              enddo

              ! get E values at prescribed frequencies
              fE = wtsB * chiB * E1val**S1B * E2val**S2B
                 
              ! compute the integrated photo-ionization and photo-heating terms
              do l=1,nnus
                 piHI(i,j,k) = piHI(i,j,k) + fHIB(l)*fE(l)/nusB(l)
                 GHI(i,j,k)  = GHI(i,j,k)  + fHIB(l)*fE(l)*(1.d0-nu0_HI/nusB(l))
              enddo

              ! get E values at prescribed frequencies
              fE = wtsC * chiC * E1val**S1C * E2val**S2C * E3val**S3C
                 
              ! compute the integrated photo-ionization and photo-heating terms
              do l=1,nnus
                 piHI(i,j,k) = piHI(i,j,k) + fHIC(l)*fE(l)/nusC(l)
                 GHI(i,j,k)  = GHI(i,j,k)  + fHIC(l)*fE(l)*(1.d0-nu0_HI/nusC(l))
              enddo 

              ! add on final coefficients
              piHI(i,j,k) = piHI(i,j,k)*Efval/hp
              GHI(i,j,k)  = GHI(i,j,k)*Efval

           enddo  ! (i loop)
        enddo  ! (j loop)
     enddo  ! (k loop)

  ! multi-species problem; compute everything for hydrogen & helium
  else

     ! iterate over the domain
     do k=1,Nz
        do j=1,Ny
           do i=1,Nx
              
              ! get time-centered values for radiation species
              Efval = Ef(i,j,k)*fsUn
              E1val = 0.5d0/Efval/chin1*(E1(i,j,k)*rUn + E1old(i,j,k)*rUn0)
              E2val = 0.5d0/Efval/chin2*(E2(i,j,k)*rUn + E2old(i,j,k)*rUn0)
              E3val = 0.5d0/Efval/chin3*(E3(i,j,k)*rUn + E3old(i,j,k)*rUn0)

              ! get E values at prescribed frequencies
              fE = wtsA * chiA * E1val**S1A

              ! compute the integrated photo-ionization and photo-heating terms
              piHI(i,j,k) = sum(fHIA*fE/nusA)
              GHI(i,j,k)  = sum(fHIA*fE*(1.d0-nu0_HI/nusA))
                
              ! get E values at prescribed frequencies
              fE = wtsB * chiB * E1val**S1B * E2val**S2B

              ! compute the integrated photo-ionization and photo-heating terms
              piHI(i,j,k)  = piHI(i,j,k) + sum(fHIB*fE/nusB)
              GHI(i,j,k)   = GHI(i,j,k)  + sum(fHIB*fE*(1.d0-nu0_HI/nusB))
              piHeI(i,j,k) = sum(fHeIB*fE/nusB)
              GHeI(i,j,k)  = sum(fHeIB*fE*(1.d0-nu0_HeI/nusB))
                
              ! get E values at prescribed frequencies
              fE = wtsC * chiC * E1val**S1C * E2val**S2C * E3val**S3C

              ! compute the integrated photo-ionization and photo-heating terms,
              ! (including the final frequency-independent coefficients)
              piHI(i,j,k)   = (piHI(i,j,k)  + sum(fHIC*fE/nusC))*Efval/hp
              GHI(i,j,k)    = (GHI(i,j,k)   + sum(fHIC*fE*(1.d0-nu0_HI/nusC)))*Efval
              piHeI(i,j,k)  = (piHeI(i,j,k) + sum(fHeIC*fE/nusC))*Efval/hp
              GHeI(i,j,k)   = (GHeI(i,j,k)  + sum(fHeIC*fE*(1.d0-nu0_HeI/nusC)))*Efval
              piHeII(i,j,k) = sum(fHeIIC*fE/nusC)*Efval/hp
              GHeII(i,j,k)  = sum(fHeIIC*fE*(1.d0-nu0_HeII/nusC))*Efval

           enddo  ! (i loop)
        enddo  ! (j loop)
     enddo  ! (k loop)

  endif  ! Nchem


!!$  print '(A,1(es9.1))', '  piHI =', sum(piHI)
!!$  print '(A,1(es9.1))', '  piHeI =', sum(piHeI)
!!$  print '(A,1(es9.1))', '  piHeII =', sum(piHeII)
!!$  print '(A,1(es9.1))', '  GHI =', sum(GHI)
!!$  print '(A,1(es9.1))', '  GHeI =', sum(GHeI)
!!$  print '(A,1(es9.1))', '  GHeII =', sum(GHeII)
!!$  print *, '  '

  return
end subroutine MFProb_ComputeRadiationIntegrals
!=======================================================================




subroutine MFProb_RadiationSpectrum(chi, nu, ESpectrum, nnus)
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
  integer, intent(in)   :: ESpectrum, nnus
  real, intent(in)  :: nu(nnus)
  real, intent(out) :: chi(nnus)
  
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
end subroutine MFProb_RadiationSpectrum
!=======================================================================




subroutine MFProb_EnforceRadiationBounds(Ef, E1, E2, E3, fsUn, rUn, &
     ESpectrum, Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr, NGzl, NGzr, ier)
!=======================================================================
!  PURPOSE: Enforces the constraints on the monochromatic radiation 
!           density values: 
!                0 < Ei < Eot(nu_i) = Ef*chi(nu_i)/chiint
!  INPUTS:
!     Ef        - free-streaming radiation
!     E1        - monochromatic radiation density at nu_1
!     E2        - monochromatic radiation density at nu_2
!     E3        - monochromatic radiation density at nu_3
!     ESpectrum - radiation spectrum type
!     *Un       - unit conversion factors
!     N*        - grid information
!  OUTPUTS:
!     ier - success/failure flag (1->success, 0->failure)
!=======================================================================
  implicit none
#include "fortran.def"

  !--------------
  ! argument declarations
  real, intent(in) :: Ef(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr)
  real, dimension(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr) &
       :: E1, E2, E3
  real, intent(in) :: fsUn, rUn
  integer, intent(in)  :: ESpectrum
  integer, intent(in)  :: Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr, NGzl, NGzr
  integer, intent(out) :: ier
  
  !--------------
  ! locals
  integer :: i, j, k, l, nbins
!  integer, parameter :: nnus=21     ! must be odd
  integer, parameter :: nnus=13     ! must be odd
  real :: hp, ev2erg, epsilon, Llimit, Ulimit
  real :: nu0_HI, nu0_HeI, nu0_HeII
  real :: chibar, chin1, chin2, chin3
  real, dimension(nnus) :: nusA, nusB, nusC, etasC
  real, dimension(nnus) :: chiA, chiB, chiC
  real, dimension(nnus) :: wtsA, wtsB, wtsC

  !=======================================================================

  ! initialize return flag
  ier = 1

  ! set constants
  hp = 6.6260693d-27            ! Planck's constant (ergs*s)
  ev2erg = 1.60217653d-12       ! conversion constant from eV to ergs
  nu0_HI   = 13.6d0*ev2erg/hp   ! ionization threshold of HI (hz)
  nu0_HeI  = 24.6d0*ev2erg/hp   ! ionization threshold of HeI (hz)
  nu0_HeII = 54.4d0*ev2erg/hp   ! ionization threshold of HeII (hz)
  epsilon = 1.d0                ! floating point roundoff
  do while ((1.d0 + epsilon) > 1.d0)  
     epsilon = epsilon*0.5d0
  end do
  Llimit = 1.d-1                ! lower limit of int (shift away from 0)
  Ulimit = 1.d0 - sqrt(epsilon) ! upper limit of int (shift away from 1)

  ! set frequencies for evaluation
  nbins = (nnus-1)/2
  do l=1,nnus
     ! nusA = linspace(nu0_HI,nu0_HeI,nnus)
     nusA(l) = nu0_HI + (l-1)*(0.9999999999d0*nu0_HeI-nu0_HI)/(nnus-1)
     ! nusB = linspace(nu0_HeI,nu0_HeII,nnus)
     nusB(l) = nu0_HeI + (l-1)*(0.9999999999d0*nu0_HeII-nu0_HeI)/(nnus-1)
     ! etasC = linspace(Llimit,Ulimit,nnus)
     etasC(l) = Llimit + (l-1)*(Ulimit-Llimit)/(nnus-1)
     ! nusC = nu0_HeII./etasC
     nusC(l) = nu0_HeII/etasC(l)
  enddo

  ! compute the quadrature weights across the intervals
  wtsA = 0.d0
  wtsB = 0.d0
  wtsC = 0.d0
  do l=1,nbins
     i = 2*l-1
     j = 2*l
     k = 2*l+1
     wtsA(i) = wtsA(i) + (nusA(k)-nusA(i))/6.d0
     wtsA(j) = wtsA(j) + (nusA(k)-nusA(i))*2.d0/3.d0
     wtsA(k) = wtsA(k) + (nusA(k)-nusA(i))/6.d0
     wtsB(i) = wtsB(i) + (nusB(k)-nusB(i))/6.d0
     wtsB(j) = wtsB(j) + (nusB(k)-nusB(i))*2.d0/3.d0
     wtsB(k) = wtsB(k) + (nusB(k)-nusB(i))/6.d0
     wtsC(i) = wtsC(i) + (etasC(k)-etasC(i))/6.d0*nu0_HeII/etasC(i)**2
     wtsC(j) = wtsC(j) + (etasC(k)-etasC(i))*2.d0/3.d0*nu0_HeII/etasC(j)**2
     wtsC(k) = wtsC(k) + (etasC(k)-etasC(i))/6.d0*nu0_HeII/etasC(k)**2
  enddo

  !    radiation spectrum across intervals
  call MFProb_RadiationSpectrum(chiA, nusA, ESpectrum, nnus)
  call MFProb_RadiationSpectrum(chiB, nusB, ESpectrum, nnus)
  call MFProb_RadiationSpectrum(chiC, nusC, ESpectrum, nnus)

  ! compute the integrated radiation spectrum
  chibar = sum(wtsA*chiA) + sum(wtsB*chiB) + sum(wtsC*chiC)

  !    radiation spectrum at ionization thresholds
  call MFProb_RadiationSpectrum(chin1, nu0_HI,   ESpectrum, 1)
  call MFProb_RadiationSpectrum(chin2, nu0_HeI,  ESpectrum, 1)
  call MFProb_RadiationSpectrum(chin3, nu0_HeII, ESpectrum, 1)

  !    rescale radiation spectra
  chin1 = chin1/chibar
  chin2 = chin2/chibar
  chin3 = chin3/chibar
  
  ! enforce bounds on monochromatic radiation fields
  E1 = max(epsilon,min(E1,Ef*fsUn*chin1/rUn))
  E2 = max(epsilon,min(E2,Ef*fsUn*chin2/rUn))
  E3 = max(epsilon,min(E3,Ef*fsUn*chin3/rUn))

  return
end subroutine MFProb_EnforceRadiationBounds
!=======================================================================




subroutine MFProb_HICrossSection(sig, nu, nnus)
!=======================================================================
!  PURPOSE: Returns the HI cross section sig(nu). 
!  INPUTS:    nu - frequencies, nnus - length of nu
!  OUTPUT:   sig - spectrum values
!  *** assumes that nu >= nu0_HI ***
!=======================================================================
  implicit none
#include "fortran.def"

  !--------------
  ! argument declarations
  integer,  intent(in)  :: nnus
  real, intent(in)  :: nu(nnus)
  real, intent(out) :: sig(nnus)
  
  !--------------
  ! locals
  integer :: i
  real*8, parameter :: pi=3.141592653589793238d0
  real*8, parameter :: nu0_HI=3.28846557762383d15    ! HI ionization (hz)
  real*8 :: eps

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
end subroutine MFProb_HICrossSection
!=======================================================================




subroutine MFProb_HeICrossSection(sig, nu, nnus)
!=======================================================================
!  PURPOSE: Returns the HeI cross section sig(nu). 
!  INPUTS:    nu - frequencies, nnus - length of nu
!  OUTPUT:   sig - spectrum values
!  *** assumes that nu >= nu0_HeI ***
!=======================================================================
  implicit none
#include "fortran.def"

  !--------------
  ! argument declarations
  integer,  intent(in)  :: nnus
  real, intent(in)  :: nu(nnus)
  real, intent(out) :: sig(nnus)
  
  !--------------
  ! locals
  integer :: i
  real*8, parameter :: pi=3.141592653589793238d0
  real*8, parameter :: nu0_HeI=5.94825391246663d15   ! HeI ionization (hz)
  real*8 :: eps

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
end subroutine MFProb_HeICrossSection
!=======================================================================




subroutine MFProb_HeIICrossSection(sig, nu, nnus)
!=======================================================================
!  PURPOSE: Returns the HeII cross section sig(nu). 
!  INPUTS:    nu - frequencies, nnus - length of nu
!  OUTPUT:   sig - spectrum values
!  *** assumes that nu >= nu0_HeI ***
!=======================================================================
  implicit none
#include "fortran.def"

  !--------------
  ! argument declarations
  integer,  intent(in)  :: nnus
  real, intent(in)  :: nu(nnus)
  real, intent(out) :: sig(nnus)
  
  !--------------
  ! locals
  integer :: i
  real*8, parameter :: pi=3.141592653589793238d0
  real*8, parameter :: nu0_HeII=1.31538623104953d16  ! HeII ionization (hz)
  real*8 :: eps

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
end subroutine MFProb_HeIICrossSection
!=======================================================================
