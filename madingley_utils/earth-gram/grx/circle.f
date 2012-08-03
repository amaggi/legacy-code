c
c*******************************************************************************
c
c    Subroutine circle
c
c    Author - Danny Harvey
c
c*******************************************************************************
c
      subroutine circle(xc,yc,r,narc,iclose,iclip,thick,ithick)
c
c    routine circle will draw either a open or closed circle (open meaning
c    a circle with a white center, closed a circle with a black center)
c
c    inputs  - xc     = x-coordinate of the center of the circle in user units
c	       yc     = y-coordinate of the center of the circle in user units
c              r      = radius of the circle in user units
c			note:  If the user units are not isotropic with respect
c			       to the physical plot map, then you will get
c			       open and closed ellipses (whoopee!)
c              narc   = number of straight line segments to approximate the
c			circle with
c              iclose = open-closed flag
c			= 0 - circle open
c			.ne. 0 - circle closed
c              iclip  = clip flag see nplot
c              thick  = line thickness see nplot
c			note: this is passed to nplot along with iclip.
c 			      This allows big open circles with thick
c			      circumferences. Probably should not be used
c			      for closed circles.
c              ithick = thickness flag see nplot
c
      dimension xbuf(100),ybuf(100)
c
      common /spc/ xll,yll,xur,yur,ASPECT,xcc,xcl,fl
c
      common /pscl/ xmin,xmax,ymin,ymax,xscale,yscale,xrange,yrange
c
      common /pdim/ xdim,ydim,xlow,ylow,rxdim,rydim,rxlow,rylow,
     1              xbl,ybl,xbh,ybh,xbm,ybm,iitran,tangle,
     2              ca,sa,cellht,cellwd,ixtype,iytype
c
      x = xmap(xc)
      y = ymap(yc)
      if (iclip .eq. 0) then
c       if (x+r*xscale .lt. xbl) return
c       if (x-r*xscale .gt. xbh) return
c       if (y+r*xscale .lt. ybl) return
c       if (y-r*xscale .gt. ybh) return
      end if
      xplt = x + xll
      yplt = y + yll
      yplt = yur - yplt
      if (xplt .lt. xll) xplt = xll
      if (yplt .lt. yll) yplt = yll
      if (xplt .gt. xur) xplt = xur
      if (yplt .gt. yur) yplt = yur
      ixplt = xplt + 0.5
      iyplt = yplt + 0.5
      ir = r*xscale + 0.5
      call hdcircle (ixplt, iyplt, ir, iclose, iclip)
      call hdcirclel (xc, yc, r, iclose, iclip)
      return
c
      end
