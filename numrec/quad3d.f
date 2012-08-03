      SUBROUTINE quad3d(x1,x2,ss)
      REAL ss,x1,x2,h
      EXTERNAL h
CU    USES h,qgausx
      call qgausx(h,x1,x2,ss)
      return
      END

      FUNCTION f(zz)
      REAL f,zz,func,x,y,z
      COMMON /xyz/ x,y,z
CU    USES func
      z=zz
      f=func(x,y,z)
      return
      END

      FUNCTION g(yy)
      REAL g,yy,f,z1,z2,x,y,z
      EXTERNAL f
      COMMON /xyz/ x,y,z
CU    USES f,qgausz,z1,z2
      REAL ss
      y=yy
      call qgausz(f,z1(x,y),z2(x,y),ss)
      g=ss
      return
      END

      FUNCTION h(xx)
      REAL h,xx,g,y1,y2,x,y,z
      EXTERNAL g
      COMMON /xyz/ x,y,z
CU    USES g,qgausy,y1,y2
      REAL ss
      x=xx
      call qgausy(g,y1(x),y2(x),ss)
      h=ss
      return
      END
