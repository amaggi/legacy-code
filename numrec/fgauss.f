      SUBROUTINE fgauss(x,a,y,dyda,na)
      INTEGER na
      REAL x,y,a(na),dyda(na)
      INTEGER i
      REAL arg,ex,fac
      y=0.
      do 11 i=1,na-1,3
        arg=(x-a(i+1))/a(i+2)
        ex=exp(-arg**2)
        fac=a(i)*ex*2.*arg
        y=y+a(i)*ex
        dyda(i)=ex
        dyda(i+1)=fac/a(i+2)
        dyda(i+2)=fac*arg/a(i+2)
11    continue
      return
      END
