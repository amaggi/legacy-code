c
c*******************************************************************************
c
c    Subroutine clrrgn
c
c    Author - Danny Harvey
c
c*******************************************************************************
c
      subroutine clrrgn(xmn,xmx,ymn,ymx)
c
c    routine clrrgn will clear out a region.
c
c
      common /pdim/ xdim,ydim,xlow,ylow,rxdim,rydim,rxlow,rylow,
     1              xbl,ybl,xbh,ybh,xbm,ybm,iitran,tangle,
     2              ca,sa,cellht,cellwd,ixtype,iytype
c
      common /spc/ xll,yll,xur,yur,ASPECT,xc,xcl,fl
c
      xmin = xmap(xmn)
      xmax = xmap(xmx)
      ymin = ymap(ymn)
      ymax = ymap(ymx)
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
      call hdclrg(ixmin,iymax,ixmax-ixmin,iymin-iymax+1)
      ypmin = ymap(ymn) + yll
      ypmax = ymap(ymx) + yll
      if (xpmin .lt. xll) xpmin = xll
      if (ypmin .lt. yll) ypmin = yll
      if (xpmin .gt. xur) xpmin = xur
      if (ypmin .gt. yur) ypmin = yur
      if (xpmax .lt. xll) xpmax = xll
      if (ypmax .lt. yll) ypmax = yll
      if (xpmax .gt. xur) xpmax = xur
      if (ypmax .gt. yur) ypmax = yur
      ixmin = xpmin*fl + 0.5
      iymin = ypmin*fl + 0.5
      ixmax = xpmax*fl + 0.5
      iymax = ypmax*fl + 0.5
      call hdclrgl(ixmin,iymin,ixmax-ixmin+1,iymax-iymin+1)
      return
      end
