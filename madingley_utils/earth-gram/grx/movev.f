      subroutine movev(x,y,itran)
c
c    routine movev moves the cursor to position x,y and then initiates
c    the graph mode.
c
      common /spc/ xll,yll,xur,yur,ASPECT,xc,xcl,fl
c
      xplt = x + xll
      yplt = yur - y - yll
      if (xplt .lt. xll) xplt = xll
      if (yplt .lt. yll) yplt = yll
      if (xplt .gt. xur) xplt = xur
      if (yplt .gt. yur) yplt = yur
      call hdmove(xplt,yplt,fl)
      yplt = y + yll
      if (yplt .lt. yll) yplt = yll
      if (yplt .gt. yur) yplt = yur
      call hdmovel(xplt,yplt,fl)
      return
      end
