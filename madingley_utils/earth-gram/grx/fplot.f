      subroutine fplot(n,x,y)
c
c    routine fplot is a fast plotting routine. This routine was written
c    to take advantage of extensions to the X server for fast vector
c    plotting.
c
c    inputs  - n      = number of points to plot
c              x      = plot data horizontal coordinate array
c              y      = plot data vertical coordinate array
c
      real*4 x(*),y(*)
c
      common /spc/ xll,yll,xur,yur,ASPECT,xc,xcl,fl
c
      common /pdim/ xdim,ydim,xlow,ylow,rxdim,rydim,rxlow,rylow,
     1              xbl,ybl,xbh,ybh,xbm,ybm,itran,tangle,
     2              ca,sa,cellht,cellwd,ixtype,iytype
c
      common /pscl/ xmin,xmax,ymin,ymax,xscale,yscale,xrange,yrange
c
      xl  =  rxlow*2.0/xbm - 1.0
      yl  =  rylow*2.0/ybm - 1.0
      xd  =  rxdim*2.0/xbm 
      yd  =  rydim*2.0/ybm
      call hdfplt (n,x,y,xl,yl,xd,yd,xmin,xmax,ymin,ymax)
      call hdstrk
      call hdstrkl
      return
      end
