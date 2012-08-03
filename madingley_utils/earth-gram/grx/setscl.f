c
c*******************************************************************************
c
c    Subroutine setscl
c
c    Author - Danny Harvey
c
c*******************************************************************************
c
      subroutine setscl(xmin,xmax,ymin,ymax)
c
c    routine setscl sets the scale factors for the plot.  It must be called
c    after the call to setdim
c
c    inputs  - xmin   = value of plot data corresponding to the left hand edge
c			of the plot
c	       xmax   = value of plot data corresponding to the right hand edge
c			of the plot
c	       ymin   = value of plot data corresponding to the bottom edge
c			of the plot
c	       ymax   = value of plot data corresponding to the top edge
c			of the plot
c
      common /pdim/ xdim,ydim,xlow,ylow,rxdim,rydim,rxlow,rylow,
     1              xbl,ybl,xbh,ybh,xbm,ybm,itran,tangle,
     2              ca,sa,cellht,cellwd,ixtype,iytype
c
      common /pscl/ xxmin,xxmax,yymin,yymax,xscale,yscale,xrange,yrange
      common /pscl2/ xrmin, xrmax, yrmin, yrmax
c
      xrmin = xmin
      xrmax = xmax
      yrmin = ymin
      yrmax = ymax
      if (ixtype .eq. 0) then
        xxmin = xmin
        xxmax = xmax
        xrange = xmax - xmin
        if (xrange .eq. 0.)  go to 900
        xscale = rxdim/xrange
      else
        if (xmin .le. 0.)  go to 920
        if (xmax .le. 0.)  go to 920
        xxmin = alog10(xmin)
        xxmax = alog10(xmax)
        xrange = xxmax - xxmin
        if (xrange .eq. 0.)  go to 900
        xscale = rxdim/xrange
      end if
c
      if (iytype .eq. 0) then
        yymin = ymin
        yymax = ymax
        yrange = ymax - ymin
        if (yrange .eq. 0.)  go to 910
        yscale = rydim/yrange
      else
        if (ymin .le. 0.)  go to 930
        if (ymax .le. 0.)  go to 930
        yymin = alog10(ymin)
        yymax = alog10(ymax)
        yrange = yymax - yymin
        if (yrange .eq. 0.)  go to 910
        yscale = rydim/yrange
      end if
c
      return
c
c    error abort
c
  900 call nclose
      call nalpha
      WRITE (6, 901)
  901 FORMAT(' SETSCL: horizontal plot range set to zero - run abo',
     1   'rted')
      stop
c
  910 call nclose
      call nalpha
      WRITE (6, 911)
  911 FORMAT(' SETSCL: vertical plot range set to zero - run abort',
     1   'ed')
      stop
c
  920 call nclose
      call nalpha
      write (6, '(/1x,a,a/)')
     1  'SETSCL: horizontal plot limit less than or equal to zero ',
     2  'for log plot'
      stop
c
  930 call nclose
      call nalpha
      write (6, '(/1x,a,a/)')
     1  'SETSCL: vertical plot limit less than or equal to zero ',
     2  'for log plot'
      stop
c
      end
      subroutine getclp (xmin, xmax, ymin, ymax)
c
      common /pdim/ xdim,ydim,xlow,ylow,rxdim,rydim,rxlow,rylow,
     1              xbl,ybl,xbh,ybh,xbm,ybm,itran,tangle,
     2              ca,sa,cellht,cellwd,ixtype,iytype
c
      xmin = rxlow
      xmax = rxlow + rxdim
      ymin = rylow
      ymax = rylow + rydim
      return
      end
      subroutine getscl (xmin, xmax, ymin, ymax)
c
      common /pscl2/ xrmin, xrmax, yrmin, yrmax
c
      xmin = xrmin
      xmax = xrmax
      ymin = yrmin
      ymax = yrmax
      return
      end
