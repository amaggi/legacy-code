c
c*******************************************************************************
c
c    Subroutine setdim
c
c    Author - Danny Harvey
c
c*******************************************************************************
c
      subroutine setdim(xdim,ydim,xlow,ylow)
c
c    routine setdim sets the physical dimensions of the plot and must be
c    called prior to a call to setscl.
c
c    inputs  - xdim   = horizontal dimension in inches
c 	       ydim   = vertical dimension in inches
c	       xlow   = horizontal location of lower left hand corner of plot
c			in inches
c	       ylow   = vertical location of lower left hand corner of plot
c			in inches
c
      common /spc/ xll,yll,xur,yur,ASPECT,xc,xcl,fl
c
      common /pdim/ xxdim,yydim,xxlow,yylow,rxdim,rydim,rxlow,rylow,
     1              xbl,ybl,xbh,ybh,xbm,ybm,itran,tangle,
     2              ca,sa,cellht,cellwd,ixtype,iytype
c
c    xc = standard raster conversion factor in rasters per inch
c
      xxdim = xdim
      yydim = ydim
      xxlow = xlow
      yylow = ylow
c
      ixdim = xdim*xc + 0.5
      iydim = ydim*xc + 0.5
      ixlow = xlow*xc + 0.5
      iylow = ylow*xc + 0.5
      rxdim = ixdim
      rydim = iydim
      rxlow = ixlow
      rylow = iylow
c
      xbm = xur - xll
      ybm = yur - yll
      xbl = rxlow
      if (xbl .lt. 0.) xbl = 0.
      ybl = rylow
      if (ybl .lt. 0.) ybl = 0.
      xbh = rxlow + rxdim
      if (xbh .gt. xbm) xbh = xbm
      ybh = rylow + rydim
      if (ybh .gt. ybm) ybh = ybm
c
      return
      end
      subroutine getfrm (xmin, xmax, ymin, ymax)
c
      common /spc/ xll,yll,xur,yur,ASPECT,xc,xcl,fl
c
      xmin = xll
      xmax = xur
      ymin = yll
      ymax = yur
c
      return
      end
      subroutine getdim (xdim,ydim,xlow,ylow)
c
      common /pdim/ xxdim,yydim,xxlow,yylow,rxdim,rydim,rxlow,rylow,
     1              xbl,ybl,xbh,ybh,xbm,ybm,itran,tangle,
     2              ca,sa,cellht,cellwd,ixtype,iytype
c
      xdim = xxdim
      ydim = yydim
      xlow = xxlow
      ylow = yylow
c
      return
      end
