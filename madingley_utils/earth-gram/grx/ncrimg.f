c
c*******************************************************************************
c
c    Subroutine ncrimg
c
c    Author - Danny Harvey
c
c*******************************************************************************
c
      subroutine ncrimg(nz, image)
c
c    routine ncrimg will create an image
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
      xmin = xmap(xxmin)
      xmax = xmap(xxmax)
      ymin = ymap(yymin)
      ymax = ymap(yymax)
      xpmin = xmin + xll
      ypmin = ymin + yll
      xpmax = xmax + xll
      ypmax = ymax + yll
      ypmin = yur - ypmin
      ypmax = yur - ypmax
      if (xpmin .lt. xll) xpmin = xll
      if (ypmin .lt. yll) ypmin = yll
      if (xpmin .gt. xur) xpmin = xur
      if (ypmin .gt. yur) ypmin = yur
      if (xpmax .lt. xll) xpmax = xll
      if (ypmax .lt. yll) ypmax = yll
      if (xpmax .gt. xur) xpmax = xur
      if (ypmax .gt. yur) ypmax = yur
      ixmin = xpmin + 0.5
      iymin = ypmin + 0.5
      ixmax = xpmax + 0.5
      iymax = ypmax + 0.5
      call hdcrimg(ixmax-ixmin+1,iymin-iymax+2,nz,image)
      return
      end
