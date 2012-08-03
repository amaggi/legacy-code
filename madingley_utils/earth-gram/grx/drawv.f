      subroutine drawv(x,y,itran)
c
c    routine drawv will draw a straight line from the existing position
c    to position x,y.
c
      common /spc/ xll,yll,xur,yur,ASPECT,xc,xcl,fl
c
      common /drwclr/ idraw
c
      xplt = x + xll
      yplt = yur - y - yll
      if (xplt .lt. xll) xplt = xll
      if (yplt .lt. yll) yplt = yll
      if (xplt .gt. xur) xplt = xur
      if (yplt .gt. yur) yplt = yur
      idraw = 1
      call hddraw(xplt,yplt,fl)
      yplt = y + yll
      if (yplt .lt. yll) yplt = yll
      if (yplt .gt. yur) yplt = yur
      call hddrawl(xplt,yplt,fl)
      return
      end
      subroutine getfl (xxl,yyl,ffl)
      common /spc/ xll,yll,xur,yur,ASPECT,xc,xcl,fl
c
      xxl = xll
      yyl = yll
      ffl = fl
      return
      end
