      FUNCTION golden(ax,bx,cx,f,tol,xmin)
      REAL golden,ax,bx,cx,tol,xmin,f,R,C
      EXTERNAL f
      PARAMETER (R=.61803399,C=1.-R)
      REAL f1,f2,x0,x1,x2,x3
      x0=ax
      x3=cx
      if(abs(cx-bx).gt.abs(bx-ax))then
        x1=bx
        x2=bx+C*(cx-bx)
      else
        x2=bx
        x1=bx-C*(bx-ax)
      endif
      f1=f(x1)
      f2=f(x2)
1     if(abs(x3-x0).gt.tol*(abs(x1)+abs(x2)))then
        if(f2.lt.f1)then
          x0=x1
          x1=x2
          x2=R*x1+C*x3
          f1=f2
          f2=f(x2)
        else
          x3=x2
          x2=x1
          x1=R*x2+C*x0
          f2=f1
          f1=f(x1)
        endif
      goto 1
      endif
      if(f1.lt.f2)then
        golden=f1
        xmin=x1
      else
        golden=f2
        xmin=x2
      endif
      return
      END
