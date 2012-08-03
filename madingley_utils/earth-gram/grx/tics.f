c
c*******************************************************************************
c
c    Subroutine tics
c
c    Author - Danny Harvey
c
c*******************************************************************************
c
      subroutine tics(x1,y1,x2,y2,stic,ntic,dtic,tlen,thick,idir)
c
c    routine tics draws tic marks
c
c    inputs  - x1     = x coordinate of starting point of axis in user units
c              y1     = y coordinate of starting point af axis in user units
c              x2     = x coordinate of ending point af axis in user units
c              y2     = y coordinate of ending point af axis in user units
c              stic   = position along axis from (x1,y1) for first tic mark
c			in user units
c              ntic   = number of tic marks
c              dtic   = increment along axis for tic marks in user units
c              tlen   = length of tic marks IN INCHES
c              thick  = thickness of tic marks IN INCHES
c              idir   = flag which defines orientation of tic marks
c			(directions given for (x2,y2) directly to the
c			right of (x1,y1))
c			> 0 - tic marks point up
c			= 0 - tic marks point on both sides
c			< 0 - tic marks point down
c
      common /spc/ xll,yll,xur,yur,ASPECT,xc,xcl,fl
c
      common /pdim/ xdim,ydim,xlow,ylow,rxdim,rydim,rxlow,rylow,
     1              xbl,ybl,xbh,ybh,xbm,ybm,itran,tangle,
     2              ca,sa,cellht,cellwd,ixtype,iytype
c
      common /pscl/ xmin,xmax,ymin,ymax,xscale,yscale,xrange,yrange
c
      common  /ocflag/  imflag,iocflg,iltp
      integer*4 ltp
c
      dimension x(3),y(3)
c
      imfl = imflag
      ltp = iltp
      imflag = 0
      if (imfl .eq. 1)  call Nopen
      if (imfl .eq. 1)  call Nlstyl(ltp)
      xr1 = xmap(x1)
      yr1 = ymap(y1)
      xr2 = xmap(x2)
      yr2 = ymap(y2)
      dx = x2 - x1
      dy = y2 - y1
      axis = sqrt(dx**2 + dy**2)
      udx = dx/axis
      udy = dy/axis
      rdx = xr2 - xr1
      rdy = yr2 - yr1
      rrdx = rdx/axis
      rrdy = rdy/axis
      raxis = sqrt(rdx**2 + rdy**2)
      rxdtic = xc*tlen*(-rdy)/raxis
      rydtic = xc*tlen*rdx/raxis
      if (idir .lt. 0) rxdtic = -rxdtic
      if (idir .lt. 0) rydtic = -rydtic
c
      do 100 i = 1,ntic
      rxtic = xmap(x1 + (stic+dtic*(i-1))*udx)
      rytic = ymap(y1 + (stic+dtic*(i-1))*udy)
      x(1) = rxtic + rxdtic
      y(1) = rytic + rydtic
      x(2) = rxtic
      y(2) = rytic
      x(3) = rxtic - rxdtic
      y(3) = rytic - rydtic
      jplot = 2
      if (idir .eq. 0) jplot = 3
      do 110  j = 1,jplot
      x(j) = (x(j) - rxlow)/xscale + xmin
      y(j) = (y(j) - rylow)/yscale + ymin
      if (ixtype .eq. 1) x(j) = 10.**x(j)
      if (iytype .eq. 1) y(j) = 10.**y(j)
  110 continue
  100 call nplot(jplot,x,y,0,1,thick,0,' ')
      if (imfl .eq. 1)  call Nclose
      imflag = imfl
c
      return
      end
