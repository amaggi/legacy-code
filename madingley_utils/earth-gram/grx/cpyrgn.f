c
c*******************************************************************************
c
c    Subroutine cpyrgn
c
c    Author - Danny Harvey
c
c*******************************************************************************
c
      subroutine cpyrgn(xmn,xmx,ymn,ymx)
c
c    routine cpyrgn will copy a region from the pixmap to the window. This is
c    usually used in conjunction with plot batching.
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
      call hdcprg(ixmin,iymax,ixmax-ixmin,iymin-iymax+1)
      return
      end
