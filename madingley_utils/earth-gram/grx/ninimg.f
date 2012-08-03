c
c*******************************************************************************
c
c    Subroutine ninimg
c
c    Author - Danny Harvey
c
c*******************************************************************************
c
      subroutine ninimg(image, nx, nxmax, ny, x, y, z, zmin, dz)
c
c    routine ninimg will interpolate a z=z(x,y) grid onto an image.
c
c
      common /pdim/ xdim,ydim,xlow,ylow,rxdim,rydim,rxlow,rylow,
     1              xbl,ybl,xbh,ybh,xbm,ybm,iitran,tangle,
     2              ca,sa,cellht,cellwd,ixtype,iytype
c
      common /pscl/ xxmin,xxmax,yymin,yymax,xscale,yscale,xrange,yrange
c
      common /spc/ xll,yll,xur,yur,ASPECT,xc,xcl,fl
c
      call hdinimg(image,nx,nxmax,ny,x,y,z,xxmin,xxmax,yymin,yymax,
     +             zmin,dz)
      return
      end
