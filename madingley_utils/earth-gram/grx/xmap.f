c
c*******************************************************************************
c
c    Subroutines xmap,ymap
c
c    Author - Danny Harvey
c
c*******************************************************************************
c
c    The following function routine maps from plot data units to
c    raster units 
c
      function xmap(x)
c
c    horizontal plot map
c
      common /pdim/ xdim,ydim,xlow,ylow,rxdim,rydim,rxlow,rylow,
     1              xbl,ybl,xbh,ybh,xbm,ybm,itran,tangle,
     2              ca,sa,cellht,cellwd,ixtype,iytype
c
      common /pscl/ xmin,xmax,ymin,ymax,xscale,yscale,xrange,yrange
c
      if (ixtype .eq. 0) then
        xmap = rxlow + (x-xmin)*xscale
      else
        xx = x
        if (x .le. 0.) xx = 1.e-30
        xmap = rxlow + (alog10(xx)-xmin)*xscale
      end if
c
      return
      end
c
      function ymap(y)
c
c    vertical plot map
c
      common /pdim/ xdim,ydim,xlow,ylow,rxdim,rydim,rxlow,rylow,
     1              xbl,ybl,xbh,ybh,xbm,ybm,itran,tangle,
     2              ca,sa,cellht,cellwd,ixtype,iytype
c
      common /pscl/ xmin,xmax,ymin,ymax,xscale,yscale,xrange,yrange
c
      if (iytype .eq. 0) then
        ymap = rylow + (y-ymin)*yscale
      else
        yy = y
        if (y .le. 0.) yy = 1.e-30
        ymap = rylow + (alog10(yy)-ymin)*yscale
      end if
c
      return
      end
      subroutine getxmp (x, xm)
      xm = xmap(x)
      return
      end
      subroutine getymp (y, ym)
      ym = ymap(y)
      return
      end
      subroutine getxscale (xsc)
c
      common /pscl/ xmin,xmax,ymin,ymax,xscale,yscale,xrange,yrange
c
      xsc = xscale
      return
      end
      function rxmap(x)
c
      common /pdim/ xdim,ydim,xlow,ylow,rxdim,rydim,rxlow,rylow,
     1              xbl,ybl,xbh,ybh,xbm,ybm,itran,tangle,
     2              ca,sa,cellht,cellwd,ixtype,iytype
c
      common /pscl/ xmin,xmax,ymin,ymax,xscale,yscale,xrange,yrange
c
      if (ixtype .eq. 0) then
        rxmap = (x - rxlow)/xscale + xmin
      else
        rxmap = 10.0**((x - rxlow)/xscale+xmin)
      end if
c
      return
      end
c
      function rymap(y)
c
      common /pdim/ xdim,ydim,xlow,ylow,rxdim,rydim,rxlow,rylow,
     1              xbl,ybl,xbh,ybh,xbm,ybm,itran,tangle,
     2              ca,sa,cellht,cellwd,ixtype,iytype
c
      common /pscl/ xmin,xmax,ymin,ymax,xscale,yscale,xrange,yrange
c
      if (iytype .eq. 0) then
        rymap = (y - rylow)/yscale + ymin
      else
        rymap = 10.0**((y - rylow)/yscale+ymin)
      end if
c
      return
      end
