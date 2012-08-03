c     $Id: interpolate.f,v 1.1 2006/02/24 09:07:57 alessia Exp $
c----------------------------------------------------------------------

      subroutine lin_interpolate(x1,x2,y1,y2,x,y)
c----------------------------------------------------------------------
c
c   Linear interpolation between two points defined by (x1,y1) and
c   (x2,y2).  Queried point has x coordinate (x).
c
c----------------------------------------------------------------------
      real x1,x2,y1,y2,x,y
      
      y=y1+(x-x1)*(y2-y1)/(x2-x1)
      
      end subroutine lin_interpolate
c----------------------------------------------------------------------

      subroutine step_interpolate(x1,x2,y1,y2,dx,x,y,n,nx,ok)
c----------------------------------------------------------------------
c     linear interpolation in regular steps of dx between x1 and x2
c----------------------------------------------------------------------
      real x1,x2,y1,y2,dx,x,y
      integer n,nx
      logical ok

      dimension x(n), y(n)


c     find number of points in interpolated array
      nx=int((x2-x1)/dx+1)

c     if the output arrays are obviously too small, return ok=false
      if (nx .gt. n) then
         ok=.false.
         write(*,*) 'Output arrays too small in step_interpolate'
         return
      endif

c     first point in output array is (x1,y1)
      x(1)=x1
      y(1)=y1
c     interpolate in equal steps for intermediate points
      do i = 2, nx-1
        x(i)=x(i-1)+dx
        call lin_interpolate(x1,x2,y1,y2,x(i),y(i))
        write(*,*) y(i)
      enddo
c     last point in output array is (x2,y2)
      x(nx)=x2
      y(nx)=y2

      ok=.true.
      return
      end subroutine step_interpolate
c----------------------------------------------------------------------


      subroutine layer_interpolate(x1,x2,y1,y2,dx,x,y,h,y_av,n,nx,ok)
c----------------------------------------------------------------------
c     Split an linear interpolation interval into layers of given
c     thickness and constant properties within each layer.  Note: 
c     nx = number of interpolated points
c     nx-1 = number of layers
c----------------------------------------------------------------------
      real x1,x2,y1,y2,dx,x,y,h,y_av
      integer n,nx
      logical ok

      dimension x(n), y(n), h(n), y_av(n)

c     do step interpolation
      call step_interpolate(x1,x2,y1,y2,dx,x,y,n,nx,ok)
      if (.not.ok) then
        write(*,*) 'Error in call to step_interpolate'
        return
      endif
      
     
c     turn interpolated points into layers with average properties
      do i = 1, nx-1
        h(i) = x(i+1)-x(i)
        y_av(i)=0.5*(y(i+1)+y(i)) 
      enddo

      ok=.true.
      return

      end subroutine layer_interpolate
c----------------------------------------------------------------------
     
