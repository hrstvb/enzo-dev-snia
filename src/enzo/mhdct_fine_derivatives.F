c#include "fortran.def"
c#include "error.def"

      subroutine mhdct_fine_derivatives(bx, by, bz, parentdim, rf,
     +                          childx, childy, childz,
     +                          cd, childstart, refinedim,
     +                          otherbx, otherby, otherbz, otherdim,
     +                          otherstart,
     +     DyBx, DzBx, DyzBx,
     +     DxBy, DzBy, DxzBy,
     +     DxBz, DyBz, DxyBz,
     +     DBxFlag, DByFlag, DBzFlag,
     +                          dx, dy, dz, method, step, cycle)


c     This routine (duh) interpolates the parent Face Centered Magnetic Field
c     to a child grid.
c     
c     Variables:
c     b[x,y,z](3)      xPartent Magnetic Field. 
c     parenddim(3)     Parent Dimension, in Parent Units.
c     parentstart(3)   Parent Start Index, in Parent Dims
c     parenttempdim(3) dimension of subgrid + 1 layer of boundary (Parent Units)
c     refine(3)        inter refinement factor (parent cell size)/(child cell size)
c     child[x,y,z]     Child Magnetic Field
c     cd(3)            Dimension of child grid.
c     childstart(3)    Position within child grid to start.  (Probably 0)
c     refinedim(3)     Dimension of region to be refined.
c     dx, dy, dz       Parent Grid Cell Size (Note: Constant over the whole subgrid.)
c                                            (If you want that changed, you have a lot
c                                            of work to do outside of this routine.)
c     method           Interpolation method.  
c     otherb[x,y,z]    For prolongation: Subgrids on the same level as child[x,y,z], 
c                          from a previous timestep.      
c                      See Balsara's AMR MHD paper for details.
c     otherdim(3)      Dimension of otherb.  As in all fortran routines, this is the
c                          cell centered dimension.
c     otherstart(3)    Index of the beginning of the overlap between the old subgrid
c                          and the new one.
c     
c     step             Portion of the interpolation routine used.  
c                      The code is broken into two parts:  The derivative taking and the reconstruction.
c                      Since this is a face centered method, the derivatives need to be globaly filled
c                      before the reconstruction happens, to ensure data consistancy between grids.
c                      Step=0 => Derivatives + Reconstruction
c                      Step=1 => derivatives only
c                      Step=2 => Reconstruction Only.
c                      If you want to see the flow in the code, look for 'goto' and 'continue'
c                      
c     Current Methods:  (dcc 01/23/06 this comment is garbage.)
c                   0: Divergence Free Quadratic interpolation by Balsara, simple sample.
c                   1:  Same routine, area weighted average.
c                   2:  undefined.
c     written David Collins, 8/4/04, 12:50 am.
c

      implicit none

c     Arguments

      integer parentdim(3)
      integer cd(3), childstart(3), refinedim(3), otherdim(3)
      integer rf(3), otherstart(3), method, step, cycle

      real  dx, dy, dz
c      REALSUB  dx, dy, dz

      real bx(parentdim(1)+1, parentdim(2), parentdim(3))
      real by(parentdim(1), parentdim(2)+1, parentdim(3))
      real bz(parentdim(1), parentdim(2), parentdim(3)+1)

      real childx(cd(1)+1, cd(2), cd(3) )
      real childy(cd(1), cd(2)+1, cd(3) )
      real childz(cd(1), cd(2), cd(3)+1 )

      real otherbx(otherdim(1)+1,otherdim(2), otherdim(3))
      real otherby(otherdim(1),otherdim(2)+1, otherdim(3))
      real otherbz(otherdim(1),otherdim(2), otherdim(3)+1)

c     local variables
c     Li, Lj, and Lk are 'local' index (which sub-cell in a parent cell)
c     Lx, Ly, and Lz are 'local' position (the center of the cell is 0, the edges
c                                         (are +- dx/2)
c     a(7,nx,ny,nz) 
c     c(7,nx,ny,nz)    coefficients for the reconstruction.
c     c(7,nx,ny,nz)    i.e. 
c     bx= a(1) + a(2)*x   + a(3)*y   + a(4)*z 
c              + a(5)*x^2 + a(6)*y*x + a(7)*z*x
c     by= b(1) + b(2)*x   + a(3)*y   + a(4)*z)
c              + b(5)*y^2 + b(6)*x*y + b(7)*z*y
c     bz= c(1) + c(2)*x)  + c(3)*y   +c(4)*z
c              + c(4)*z^2 + c(6)*x*z + c(7)*y*z
c     PLUS HIGHER TERMS.  NEEDS MORE COMMENTS.
c     DyBx (etc) Y derivative of Bx. (etc.)
c                Derivatives taken using the MinMod slope limiter,
c                for TVD.
c     DxBy, etc       Derivatives for reconstruction. 
c                     Note: the size of coefficients and derivatives
c                           is the size of the parent grid getting refined
c                           This is smaller than the portion of the parent grid
c                           that's actually passed in, because the refinement stencil
c	                    is more than a single cell.
c     Li, Lj, Lk      integer position of subcell within parent cell
c     Pi, Pj, Pk      index of parent position, global parent units 
c     Lx, Ly, Lz      Local position of this cell. (cell boundaries are [-dx/2, +dx/2]
c     oi, oj, ok      Position of 'other' subgrid
c     offset2(3)      Starting position of child grid w/ respect to the Parent grid.(subgrid units)
c     offset1(3)      Starting position of coefficients wrt ParentGrid.  (Actually, 1.)

      real DyBx(1+cd(1)/rf(1)+1, 1+cd(2)/rf(2), 1+cd(3)/rf(3))
      real DzBx(1+cd(1)/rf(1)+1, 1+cd(2)/rf(2), 1+cd(3)/rf(3))
      real DyzBx(1+cd(1)/rf(1)+1, 1+cd(2)/rf(2), 1+cd(3)/rf(3))
      integer DBxFlag(1+cd(1)/rf(1)+1, 1+cd(2)/rf(2), 1+cd(3)/rf(3))

      real DxBy(1+cd(1)/rf(1), 1+cd(2)/rf(2)+1, 1+cd(3)/rf(3))
      real DzBy(1+cd(1)/rf(1), 1+cd(2)/rf(2)+1, 1+cd(3)/rf(3))
      real DxzBy(1+cd(1)/rf(1), 1+cd(2)/rf(2)+1, 1+cd(3)/rf(3))
      integer DByFlag(1+cd(1)/rf(1), 1+cd(2)/rf(2)+1, 1+cd(3)/rf(3))

      real DxBz(1+cd(1)/rf(1), 1+cd(2)/rf(2), 1+cd(3)/rf(3)+1)
      real DyBz(1+cd(1)/rf(1), 1+cd(2)/rf(2), 1+cd(3)/rf(3)+1)
      real DxyBz(1+cd(1)/rf(1), 1+cd(2)/rf(2), 1+cd(3)/rf(3)+1)
      integer DBzFlag(1+cd(1)/rf(1), 1+cd(2)/rf(2), 1+cd(3)/rf(3)+1)
      
      real a(11,1+cd(1)/rf(1), 1+cd(2)/rf(2), 1+cd(3)/rf(3))
      real b(11,1+cd(1)/rf(1), 1+cd(2)/rf(2), 1+cd(3)/rf(3))
      real c(11,1+cd(1)/rf(1), 1+cd(2)/rf(2), 1+cd(3)/rf(3))

      integer i, j, k, l, i1, i2, j1, j2, k1, k2, m
      integer Li, Lj, Lk, Pi, Pj, Pk, oi, oj, ok
      integer offset1(3), offset2(3), t, pstart(3), s(3)
      integer prdim(3), pend(3)
      integer file0, file1, file2,file3, file4
      integer TestX, TestY, TestZ
c      character*15 filename
      real Lx, Ly, Lz, dxi, dyi, dzi
c      REALSUB Lx, Ly, Lz, dxi, dyi, dzi
      real minmod
      real tol
      integer NanLoo
c     this is just what it sounds like:  it's a dummy variable that holds things together while
c     I debug, before I rip this routine appart.
      integer methodkludge

      tol = 1e-9
      dxi = 1.0/dx
      dyi = 1.0/dy
      dzi = 1.0/dz
      
      do i=1, 3
         offset1(i) = 1

c     Offset2 is the number of child cells between the first child cell and a parent cell edge.
c     Since an odd number of ghost zones may be used, the physical location of the edge of the 
c     ghost zones of the subgrid may not line up with a parent grid edge.
c     offset2 reflects this space.

c     On the actual calculation:  the first line is the number of cells difference.  
c     The second line removes the boarder of parent cells.

         offset2(i) = (rf(i)*ParentDim(i)- cd(i) )/2
         offset2(i) = offset2(i) - rf(i)

c     pstart is the starting position of the coefficient arrays.
         pstart(i) = (childstart(i)+offset2(i))/rf(i)+1

c     s() is a shift, used to shift the index for certain prolongation faces.
         s(i) = 0

c     prdim is the dimension of the coeffecient grid to be filled.
c     the +1 is to take care of rounding errors.

         prdim(i) = max( (refinedim(i)+1)/rf(i),1)
         pend(i) = pstart(i) + prdim(i) -1
      enddo
      

c     A Negative method flag indicates prolongation to newly interpolated fields.
c     the variable s takes care of a necessary position offset.  See documentation maps
c     If david hasn't written the documentation map yet, throw eggs at him.

      if(method .eq. -1) s(1) = 1
      if(method .eq. -3) s(2) = 1
      if(method .eq. -5) s(3) = 1



c
c     Fill derivaives of Bx
c


c     Prolonging a left x face of the other grid
c     This is index hell.  Remember that DyzBx lives on the Parent Grid,
c     which is generated in Grid_InterpolateFieldValues.
c     Also remember that otherstart(i) comes from c++, which starts arrays at 0.
c     Also remember that Sonic Youth is important for you to listen to.

c      if( (method .eq. -1 .or. method .eq. -2) .and. 1.eq.0) then


      if( (method .eq. -1 .or. method .eq. -2)  ) then
c dccf
         if( method .eq. -1 ) then
            i = pstart(1) + 1
         endif
         if( method .eq. -2 ) then
            i = pstart(1) 
         endif

         oi = otherstart(1)+1

         j1 = pstart(2)
         j2 = pend(2)
         
         k1 = pstart(3)
         k2 = pend(3)
         
         do k=k1,k2
            do j=j1,j2
               
               ok = rf(3)*(k-k1) + otherstart(3) + 1
               oj = rf(2)*(j-j1) + otherstart(2) + 1

               DBxFlag(i,j,k) = DBxFlag(i,j,k) + 1

               DyzBx(i,j,k) = (
     +               otherbx(oi,oj+1,ok+1)
     +              - otherbx(oi,oj  ,ok+1)
     +              + otherbx(oi,oj  ,ok  )
     +              - otherbx(oi,oj+1,ok  )
     +              )
     +              *dzi*dyi*rf(2)*rf(3)

               DyBx(i,j,k) = (
     +               otherbx(oi,oj+1,ok+1)
     +              - otherbx(oi,oj  ,ok+1)
     +              + otherbx(oi,oj+1,ok  )
     +              - otherbx(oi,oj  ,ok  )
     +              )
     +              *dyi*rf(2)*0.5

               DzBx(i,j,k) = (
     +               otherbx(oi,oj+1,ok+1)
     +              - otherbx(oi,oj+1,ok  )
     +              + otherbx(oi,oj  ,ok+1)
     +              - otherbx(oi,oj  ,ok  )
     +              )
     +              *dzi*rf(3)*0.5


               
      enddo
      enddo


      endif      


      if( method .eq. -3 .or. method .eq. -4 ) then

         i1 = pstart(1)
         i2 = pend(1)

c         j1 = pstart(2)
c         j2 = j1 + prdim(2) -1
         
         k1 = pstart(3)
         k2 = pend(3)
         
         if( method .eq. -3 ) then
            j = pstart(2) + 1

         endif
         if( method .eq. -4 ) then
            j = pstart(2) 
         endif

            oj = otherstart(2) +1

c         oi = otherstart(1)
c         ok = otherstart(3)
         
         do k=k1,k2
c            do j=j1, j2
            do i=i1, i2

               oi = rf(1)*(i-i1) + otherstart(1) + 1
c               oj = rf(2)*(j-j1) + otherstart(2) + 1
               ok = rf(3)*(k-k1) + otherstart(3) + 1

               DByFlag(i,j,k) = DByFlag(i,j,k) + 1

               DxzBy(i,j,k) = (
     +              otherby(oi+1, oj, ok+1)
     +              -otherby(oi, oj, ok+1)
     +              +otherby(oi, oj, ok)
     +             - otherby(oi+1, oj, ok)
     +              )
     +              *dxi*dzi*rf(1)*rf(3)

               DxBy(i,j,k) = (
     +              otherby(oi+1, oj, ok+1)
     +              -otherby(oi, oj, ok+1)
     +              +otherby(oi+1, oj, ok)
     +              -otherby(oi, oj, ok)
     +              )
     +              *dxi*rf(1)*0.5

               DzBy(i,j,k) = (
     +              otherby(oi+1, oj, ok+1)
     +              -otherby(oi+1, oj, ok)
     +              +otherby(oi, oj, ok+1)
     +              -otherby(oi, oj, ok)
     +              )
     +              *dzi*rf(3)*0.5


      enddo
      enddo

      endif      

      if( method .eq. -5 .or. method .eq. -6 ) then

         i1 = pstart(1)
         i2 = pend(1)

         j1 = pstart(2)
         j2 = pend(2)
         
c         k1 = pstart(3)
c         k2 = k1 + prdim(3) -1
         
         if( method .eq. -5 ) then
            k = pstart(3) + 1
         endif
         if( method .eq. -6 ) then
            k = pstart(3) 
         endif
         ok = otherstart(3)+1

c         oi = otherstart(1)
c         oj = otherstart(2)

         
c         do k=k1,k2
            do j=j1, j2
            do i=i1, i2

               oi = rf(1)*(i-i1) + otherstart(1) + 1
               oj = rf(2)*(j-j1) + otherstart(2) + 1
c               ok = rf(3)*(k-k1) + otherstart(3) + 1

               DBzFlag(i,j,k) = DBzFlag(i,j,k) + 1

               DxyBz(i,j,k) = (  
     +               otherbz(oi+1,oj+1,ok)
     +              -otherbz(oi  ,oj+1,ok)
     +              +otherbz(oi  ,oj  ,ok)
     +              -otherbz(oi+1,oj  ,ok)
     +              )
     +              *dxi*dyi*rf(1)*rf(2)


               DxBz(i,j,k) = (  
     +               otherbz(oi+1,oj+1,ok)
     +              -otherbz(oi  ,oj+1,ok)
     +              +otherbz(oi+1,oj  ,ok)
     +              -otherbz(oi  ,oj  ,ok)
     +              )
     +              *dxi*rf(1)*0.5

               DyBz(i,j,k) = (  
     +               otherbz(oi+1,oj+1,ok)
     +              -otherbz(oi+1,oj  ,ok)
     +              +otherbz(oi  ,oj+1,ok)
     +              -otherbz(oi  ,oj  ,ok)
     +              )
     +              *dyi*rf(2)*0.5


c     See the Bx=<bx> comment.

               if( cycle .eq. -12 ) then
                  Bz(i+offset1(1),j+offset1(2),k+offset1(3)) = (  
     +                 otherbz(oi+1,oj+1,ok)+
     +                 otherbz(oi+1,oj  ,ok)+
     +                 otherbz(oi  ,oj+1,ok)+
     +                 otherbz(oi  ,oj  ,ok)
     +              )*0.25


               endif


      enddo
      enddo

c     endif method= -5 or -6
      endif      

 666  continue
      end  
      
