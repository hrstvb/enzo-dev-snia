!=======================================================================
!
! Copyright 2009 Daniel R. Reynolds
!
! This software is released under the terms of the "Enzo Public License"
! in the accompanying LICENSE file.
!
!=======================================================================
subroutine MFSplit_ComputeRadiationIntegrals(piHI, piHeI, piHeII, GHI, &
     GHeI, GHeII, Ef, E1, E2, E3, E1old, E2old, E3old, Nchem,          &
     ESpectrum, chibar, fsUn, Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr,      &
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
  real, intent(in) :: fsUn, chibar
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
  real :: Efval, E1val, E2val, E3val
  real :: s1n1, s1n2, s1n3, s2n2, s2n3, s3n3
  real, dimension(nnus) :: nusA, nusB, nusC, etasC
  real, dimension(nnus) :: sHIA, sHIB, sHIC
  real, dimension(nnus) :: sHeIB, sHeIC, sHeIIC
  real, dimension(nnus) :: fnuA, fnuB, fnuC
  real, dimension(nnus) :: chiA, chiB, chiC
  real, dimension(nnus) :: S1A, S1B, S1C
  real, dimension(nnus) :: S2B, S2C, S3C
  real, dimension(nnus) :: fE
  real, dimension(nnus) :: wtsA, wtsB, wtsC

  !=======================================================================

!!$  print *, '  '
!!$  print *, 'Entering MFSplit_ComputeRadiationIntegrals'

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

  ! evaluate frequency-dependent functions at nu values
  !    species cross-sections at ionization thresholds
  call MFSplit_HICrossSection(s1n1, nu0_HI, 1)
  call MFSplit_HICrossSection(s1n2, nu0_HeI,1)
  call MFSplit_HICrossSection(s1n3, nu0_HeII,1)
  call MFSplit_HeICrossSection(s2n2, nu0_HeI, 1)
  call MFSplit_HeICrossSection(s2n3, nu0_HeII, 1)
  call MFSplit_HeIICrossSection(s3n3, nu0_HeII, 1)
  
  !    species cross sections across intervals
  call MFSplit_HICrossSection(sHIA, nusA, nnus)
  call MFSplit_HICrossSection(sHIB, nusB, nnus)
  call MFSplit_HICrossSection(sHIC, nusC, nnus)
  call MFSplit_HeICrossSection(sHeIB, nusB, nnus)
  call MFSplit_HeICrossSection(sHeIC, nusC, nnus)
  call MFSplit_HeIICrossSection(sHeIIC, nusC, nnus)

  !    radiation spectrum across intervals
  call MFSplit_RadiationSpectrum(chiA, nusA, ESpectrum, nnus)
  call MFSplit_RadiationSpectrum(chiB, nusB, ESpectrum, nnus)
  call MFSplit_RadiationSpectrum(chiC, nusC, ESpectrum, nnus)

  !    exponents for radiation terms
  S3C = sHeIIC/s3n3
  S2B = sHeIB/s2n2
  S2C = sHeIC/s2n2 - s2n3/s2n2*S3C
  S1A = sHIA/s1n1
  S1B = sHIB/s1n1  - s1n2/s1n1*S2B
  S1C = sHIC/s1n1  - s1n2/s1n1*S2C - s1n3/s1n1*S3C

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

  ! single-species problem; compute everything for hydrogen only
  if (Nchem == 1) then

     ! iterate over the domain
     do k=1,Nz
        do j=1,Ny
           do i=1,Nx

              ! get time-centered values for radiation species
              ! note: E1,E2,E3 units already take the chiX/chibar into account
              Efval = Ef(i,j,k)*fsUn
              E1val = 0.5d0*(E1(i,j,k) + E1old(i,j,k))/Ef(i,j,k)
              E2val = 0.5d0*(E2(i,j,k) + E2old(i,j,k))/Ef(i,j,k)
              E3val = 0.5d0*(E3(i,j,k) + E3old(i,j,k))/Ef(i,j,k)

              ! initialize the outputs
              piHI(i,j,k) = 0.d0
              GHI(i,j,k)  = 0.d0

              ! get E values at prescribed frequencies
              fE = wtsA * chiA * E1val**S1A
                 
              ! compute the integrated photo-ionization and photo-heating terms
              do l=1,nnus
                 piHI(i,j,k) = piHI(i,j,k) + sHIA(l)*fE(l)/nusA(l)
                 GHI(i,j,k)  = GHI(i,j,k)  + sHIA(l)*fE(l)*(1.d0-nu0_HI/nusA(l))
              enddo

              ! get E values at prescribed frequencies
              fE = wtsB * chiB * E1val**S1B * E2val**S2B
                 
              ! compute the integrated photo-ionization and photo-heating terms
              do l=1,nnus
                 piHI(i,j,k) = piHI(i,j,k) + sHIB(l)*fE(l)/nusB(l)
                 GHI(i,j,k)  = GHI(i,j,k)  + sHIB(l)*fE(l)*(1.d0-nu0_HI/nusB(l))
              enddo

              ! get E values at prescribed frequencies
              fE = wtsC * chiC * E1val**S1C * E2val**S2C * E3val**S3C
                 
              ! compute the integrated photo-ionization and photo-heating terms
              do l=1,nnus
                 piHI(i,j,k) = piHI(i,j,k) + sHIC(l)*fE(l)/nusC(l)
                 GHI(i,j,k)  = GHI(i,j,k)  + sHIC(l)*fE(l)*(1.d0-nu0_HI/nusC(l))
              enddo 

              ! add on final coefficients
              piHI(i,j,k) = piHI(i,j,k)*c*Efval/hp/chibar
              GHI(i,j,k)  = GHI(i,j,k)*c*Efval/chibar


              if ((i==int(Nx/2)) .and. (j==int(Ny/2)) .and. (k==int(Nz/2))) then
                 print *, ' Far from ionization source:'
                 print *, '    Ef =',Ef(i,j,k)
                 print *, '    E1 =',E1(i,j,k)
                 print *, '    E2 =',E2(i,j,k)
                 print *, '    E3 =',E3(i,j,k)
                 print *, '    piHI =',piHI(i,j,k)
                 print *, '     GHI =',GHI(i,j,k)
              end if

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
              ! note: E1,E2,E3 units already take the chiX/chibar into account
              Efval = Ef(i,j,k)*fsUn
              E1val = 0.5d0*(E1(i,j,k) + E1old(i,j,k))/Ef(i,j,k)
              E2val = 0.5d0*(E2(i,j,k) + E2old(i,j,k))/Ef(i,j,k)
              E3val = 0.5d0*(E3(i,j,k) + E3old(i,j,k))/Ef(i,j,k)

              ! initialize the outputs
              piHI(i,j,k)   = 0.d0
              piHeI(i,j,k)  = 0.d0
              piHeII(i,j,k) = 0.d0
              GHI(i,j,k)    = 0.d0
              GHeI(i,j,k)   = 0.d0
              GHeII(i,j,k)  = 0.d0

              ! get E values at prescribed frequencies
              fE = wtsA * chiA * E1val**S1A

              ! compute the integrated photo-ionization and photo-heating terms
              piHI(i,j,k) = piHI(i,j,k) + sum(sHIA*fE/nusA)
              GHI(i,j,k)  = GHI(i,j,k)  + sum(sHIA*fE*(1.d0-nu0_HI/nusA))
                
              ! get E values at prescribed frequencies
              fE = wtsB * chiB * E1val**S1B * E2val**S2B

              ! compute the integrated photo-ionization and photo-heating terms
              piHI(i,j,k)  = piHI(i,j,k)  + sum(sHIB*fE/nusB)
              piHeI(i,j,k) = piHeI(i,j,k) + sum(sHeIB*fE/nusB)
              GHI(i,j,k)   = GHI(i,j,k)   + sum(sHIB*fE*(1.d0-nu0_HI/nusB))
              GHeI(i,j,k)  = GHeI(i,j,k)  + sum(sHeIB*fE*(1.d0-nu0_HeI/nusB))
                
              ! get E values at prescribed frequencies
              fE = wtsC * chiC * E1val**S1C * E2val**S2C * E3val**S3C

              ! compute the integrated photo-ionization and photo-heating terms,
              ! (including the final frequency-independent coefficients)
              piHI(i,j,k)   = piHI(i,j,k)   + sum(sHIC*fE/nusC)
              piHeI(i,j,k)  = piHeI(i,j,k)  + sum(sHeIC*fE/nusC)
              piHeII(i,j,k) = piHeII(i,j,k) + sum(sHeIIC*fE/nusC)
              GHI(i,j,k)    = GHI(i,j,k)    + sum(sHIC*fE*(1.d0-nu0_HI/nusC))
              GHeI(i,j,k)   = GHeI(i,j,k)   + sum(sHeIC*fE*(1.d0-nu0_HeI/nusC))
              GHeII(i,j,k)  = GHeII(i,j,k)  + sum(sHeIIC*fE*(1.d0-nu0_HeII/nusC))

              ! add on final coefficients
              piHI(i,j,k)   = c*Efval/hp/chibar*piHI(i,j,k)
              piHeI(i,j,k)  = c*Efval/hp/chibar*piHeI(i,j,k)
              piHeII(i,j,k) = c*Efval/hp/chibar*piHeII(i,j,k)
              GHI(i,j,k)    = c*Efval/chibar*GHI(i,j,k)
              GHeI(i,j,k)   = c*Efval/chibar*GHeI(i,j,k)
              GHeII(i,j,k)  = c*Efval/chibar*GHeII(i,j,k)

           enddo  ! (i loop)
        enddo  ! (j loop)
     enddo  ! (k loop)

  endif  ! Nchem


  return
end subroutine MFSplit_ComputeRadiationIntegrals
!=======================================================================
