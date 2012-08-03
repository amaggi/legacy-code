      subroutine tdraw(x,y,ifirst,itran)
c
c    routine tdraw plots a single point
c
c    inputs  - x      = x coordinate in rasters
c	       y      = y coordinate in rasters
c	       ifirst = 1 - first point
c			.ne. 1 - not first point
c	       itran  = plot mode flag
c
      if (ifirst .eq. 1) call movev(x,y,itran)
      call drawv(x,y,itran)
      ifirst = 0
c
      return
      end
