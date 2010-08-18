#ifdef USE_MKL
#include "../fortran.def"
c=======================================================================
c///////////////////////  SUBROUTINE INTVAR  \\\\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine intvar(qslice, idim, i1, i2, isteep, steepen, iflatten, 
     &                  flatten, c1, c2, c3, c4, c5, c6, char1, char2,
     &                  c0, dq, ql, qr, q6, qla, qra, ql0, qr0)
c
c  COMPUTES LEFT AND RIGHT EULERIAN INTERFACE VALUES FOR RIEMANN SOLVER
c
c  written by: Greg Bryan
c  date:       March, 1996
c  modified1:
c
c  PURPOSE:  Uses piecewise parabolic interpolation to compute left-
c    and right interface values to be fed into Riemann solver during a
c    one dimensional sweeps.  This version computes the Eulerian corrections
c    to the left and right states described in section three of Colella &
c    Woodward (1984), JCP.  The routine works on a single variable in
c    one dimension.
c
c  INPUT:
c    qslice   - one dimensional field of quantity q (one of d,e,u,v...)
c    idim     - declared dimension of 1D fields
c    i1, i2   - start and end indexes of active region
c    isteep   - steepening flag (1 = on, 0 = off); only apply to density!
c    steepen    - steepening coefficients
c    iflatten - flattening flag (1 = on, 0 = off)
c    flatten  - flattening coefficients
c    c1-6     - precomputed grid coefficients
c    char1,2  - characteristic distances for +/- waves (for average)
c    c0       - characteristic distance (for lagrangean cell face)
c    dq, ql, qr, q6 - 1D field temporaries
c    
c  OUTPUT:
c    qla, qra - left and right state values (from char1,2)
c    ql0, qr0 - left and right state values (from c0)
c
c  EXTERNALS:
c
c  LOCALS:
c
c  PARAMETERS:
c    ft     - a constant used in eq. 1.124 (=2*2/3)
c
c-----------------------------------------------------------------------
c
      implicit NONE
c
c  argument declarations
c
      integer idim, i1, i2, ni, iflatten, isteep
      real c1(idim), c2(idim), c3(idim), c4(idim), c5(idim), c6(idim),
     &     char1(idim), char2(idim), c0(idim),
     &     qla(idim), qra(idim), ql0(idim), qr0(idim)
      real qslice(idim), steepen(idim), flatten(idim)
c
c  parameters
c
      real ft
      parameter(ft = 4.0/3.0)
c
c  local declarations (arrays passed as temps)
c
      integer i
      real qdiff(idim), temp1(idim), temp2(idim), temp3(idim), 
     $     temp4(idim), temp5(idim), temp22(idim), temp23(idim)
      real tempa(idim), tempb(idim), tempc(idim)
      real dq(idim), ql(idim), qr(idim), q6(idim)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////
c=======================================================================
c
c     Compute average linear slopes (eqn 1.7)
c      Monotonize (eqn 1.8)
c
c      do i = i1-2, i2+2
c         qplus = qslice(i+1)-qslice(i  )
c         qmnus = qslice(i  )-qslice(i-1)
c         dq(i) = c1(i)*qplus + c2(i)*qmnus
c         temp1 = min(abs(dq(i)), 2.0*abs(qmnus), 2.0*abs(qplus))
c         if (qplus*qmnus .gt. 0) then
c            dq(i) = temp1*sign(1.0, dq(i))
c         else
c            dq(i) = 0.0
c         endif
c      enddo
c
      ni = i2-i1+5
      MKL_SUB(ni, qslice(i1-1), qslice(i1-2), qdiff(i1-2))
      MKL_MUL(ni, c1(i1-2), qdiff(i1-1), temp1(i1-2))
      MKL_MUL(ni, c2(i1-2), qdiff(i1-2), temp2(i1-2))
      MKL_ADD(ni, temp1(i1-2), temp2(i1-2), dq(i1-2))
      MKL_MUL(ni, qdiff(i1-1), qdiff(i1-2), temp3(i1-2))
      do i = i1-2, i2+2
         if (temp3 .gt. 0) then
            dq(i) = temp1(i) * sign(1.0, dq(i))
         else
            dq(i) = 0.0
         endif
      enddo
c     
c     Construct left and right values (eqn 1.6)
c
c      do i = i1-1, i2+2
c         ql(i) = c3(i)*qslice(i-1) + c4(i)*qslice(i) +
c     &           c5(i)*    dq(i-1)   + c6(i)*dq(i)
c         qr(i-1) = ql(i)
c      enddo
c
      ni = i2-i1+4
      MKL_MUL(ni, c3(i1-1), qslice(i1-2), temp1)
      MKL_MUL(ni, c4(i1-1), qslice(i1-1), temp2)
      MKL_MUL(ni, c5(i1-1), dq(i1-2), temp3)
      MKL_MUL(ni, c6(i1-1), dq(i1-1), temp4)
      do i = i1-1, i2+2
         ql(i) = 0.0
      enddo
      MKL_ADD(ni, temp1,    temp2, ql(i1-1))
      MKL_ADD(ni, ql(i1-1), temp3, ql(i1-1))
      MKL_ADD(ni, ql(i1-1), temp4, ql(i1-1))
      do i = i1-1, i2+2
         qr(i-1) = ql(i)
      enddo
c
c     Steepen if asked for (use precomputed steepening parameter)
c
      if (isteep .ne. 0) then
         do i = i1-1, i2+1
            ql(i) = (1.0-steepen(i))*ql(i) + 
     &              steepen(i)*(qslice(i-1)+0.5*dq(i-1))
            qr(i) = (1.0-steepen(i))*qr(i) + 
     &              steepen(i)*(qslice(i+1)-0.5*dq(i+1))
         enddo
      endif
c
c     Monotonize again (eqn 1.10)
c
      ni = i2-i1+3
c      do i=i1-1,i2+1
c         temp1 = (qr(i)-qslice(i))*(qslice(i)-ql(i))

      MKL_SUB(ni, qr(i1-1), qslice(i1-1), tempa)
      MKL_SUB(ni, qslice(i1-1), ql(i1-1), tempb)
      MKL_MUL(ni, tempa, tempb, temp1(i1-1))

c         temp2 = qr(i)-ql(i)

      MKL_SUB(ni, qr(i1-1), ql(i1-1), temp2(i1-1))

c         temp3 = 6.0*(qslice(i)-0.5*(qr(i)+ql(i)))

      MKL_ADD(ni, qr(i1-1), ql(i1-1), tempa)
      do i=1,ni-1
         tempa(i) = 0.5 * tempa(i)
      enddo
      MKL_SUB(ni, qslice(i1-1), tempa, temp3(i1-1))
      do i=i1-1,i2+1
         temp3(i) = 6.0 * tempa(i)
      enddo

c         if (temp1 .le. 0.0) then
c            ql(i) = qslice(i)
c            qr(i) = qslice(i)
c         endif

      do i=i1-1,i2+1
         if (temp1(i) .le. 0.0) then
            ql(i) = qslice(i)
            qr(i) = qslice(i)
         endif
      enddo

c         temp22 = temp2**2
c         temp23 = temp2*temp3

      MKL_SQR(ni, temp2(i1-1), temp22(i1-1))
      MKL_MUL(ni, temp2(i1-1), temp3(i1-1), temp23(i1-1))

c         if (temp22 .lt. temp23)
c     &        ql(i) = 3.0*qslice(i) - 2.0*qr(i)
c         if (temp22 .lt. -temp23)
c     &        qr(i) = 3.0*qslice(i) - 2.0*ql(i)
c      enddo
c

      do i=i1-1,i2+1
         if (temp22(i) .lt. temp23(i))
     $        ql(i) = 3.0*qslice(i) - 2.0*qr(i)
         if (temp22(i) .lt. -temp23(i))
     $        qr(i) = 3.0*qslice(i) - 2.0*ql(i)
      enddo

c
c     If requested, flatten slopes with flatteners calculated in calcdiss (4.1)
c
      if (iflatten .ne. 0) then
         do i = i1-1, i2+1
            ql(i) = qslice(i)*flatten(i) + ql(i)*(1.0-flatten(i))
            qr(i) = qslice(i)*flatten(i) + qr(i)*(1.0-flatten(i))
         enddo
      endif
c
c    Now construct left and right interface values (eqn 1.12 and 3.3)
c
      ni = i2-i1+3
c      do i = i1-1, i2+1
c         q6(i) = 6.0*(qslice(i)-0.5*(ql(i)+qr(i)))

      MKL_ADD(ni, qr(i1-1), ql(i1-1), tempa)
      do i=1,ni-1
         tempa(i) = 0.5 * tempa(i)
      enddo
      MKL_SUB(ni, qslice(i1-1), tempa, q6(i1-1))
      do i=i1-1,i2+1
         q6(i) = 6.0 * q6(i)
      enddo

c         dq(i) = qr(i) - ql(i)
c      enddo

      MKL_SUB(ni, qr(i1-1), ql(i1-1), dq)
c
c
      ni = i2-i1+2
c      do i = i1, i2+1
c        qla(i)= qr(i-1)-char1(i-1)*(dq(i-1)-(1.0-ft*char1(i-1))*q6(i-1))
      
      do i = i1, i2+1
         temp1(i) = 1.0 - ft * char1(i-1)
      enddo
      MKL_SUB(ni, dq(i1-1), temp1(i1), temp2(i1))
      MKL_MUL(ni, char1(i1-1), temp2(i1), temp3(i1))
      MKL_MUL(ni, temp3(i1), q6(i1-1), temp4(i1))
      MKL_SUB(ni, qr(i1-1), temp4(i1), qla(i1))

c        qra(i)= ql(i  )+char2(i  )*(dq(i  )+(1.0-ft*char2(i  ))*q6(i  ))
c      enddo

      do i = i1, i2+1
         temp1(i) = 1.0 - ft * char2(i)
      enddo
      MKL_SUB(ni, dq(i1), temp1(i1), temp2(i1))
      MKL_MUL(ni, char2(i1), temp2(i1), temp3(i1))
      MKL_MUL(ni, temp3(i1), q6(i1), temp4(i1))
      MKL_ADD(ni, ql(i1), temp4(i1), qra(i1))

c
c      do i=i1,i2+1
c         ql0(i) = qr(i-1)-c0(i-1)*(dq(i-1)-(1.0-ft*c0(i-1))*q6(i-1))

      do i = i1, i2+1
         temp1(i) = 1.0 - ft * c0(i-1)
      enddo
      MKL_SUB(ni, dq(i1-1), temp1(i1), temp2(i1))
      MKL_MUL(ni, c0(i1-1), temp2(i1), temp3(i1))
      MKL_MUL(ni, temp3(i1), q6(i1-1), temp4(i1))
      MKL_SUB(ni, qr(i1-1), temp4(i1), ql0(i1))

c         qr0(i) = ql(i  )-c0(i  )*(dq(i  )+(1.0+ft*c0(i  ))*q6(i  ))
c      enddo

      do i = i1, i2+1
         temp1(i) = 1.0 + ft * c0(i)
      enddo
      MKL_SUB(ni, dq(i1), temp1(i1), temp2(i1))
      MKL_MUL(ni, c0(i1), temp2(i1), temp3(i1))
      MKL_MUL(ni, temp3(i1), q6(i1), temp4(i1))
      MKL_SUB(ni, ql(i1), temp4(i1), qr0(i1))

c
      return
      end

#endif /* USE_MKL */
