      SUBROUTINE kermom(w,y,m)
      INTEGER m
      DOUBLE PRECISION w(m),y,x,d,df,clog,x2,x3,x4
      COMMON /momcom/ x
      if (y.ge.x) then
        d=y-x
        df=2.d0*sqrt(d)*d
        w(1)=df/3.d0
        w(2)=df*(x/3.d0+d/5.d0)
        w(3)=df*((x/3.d0 + 0.4d0*d)*x + d**2/7.d0)
        w(4)=df*(((x/3.d0 + 0.6d0*d)*x + 3.d0*d**2/7.d0)*x+ d**3/9.d0)
      else
        x2=x**2
        x3=x2*x
        x4=x2*x2
        d=x-y
        clog=log(d)
        w(1)=d*(clog-1.d0)
        w(2)=-0.25d0*(3.d0*x+y-2.d0*clog*(x+y))*d
        w(3)=(-11.d0*x3+y*(6.d0*x2+y*(3.d0*x+2.d0*y))+6.d0*clog*(x3-y**
     *3))/18.d0
        w(4)=(-25.d0*x4+y*(12.d0*x3+y*(6.d0*x2+y*(4.d0*x+3.d0*y)))+
     *12.d0*clog*(x4-y**4))/48.d0
      endif
      return
      END
