      SUBROUTINE chebft(a,b,c,n,func)
      INTEGER n,NMAX
      REAL a,b,c(n),func,PI
      EXTERNAL func
      PARAMETER (NMAX=50, PI=3.141592653589793d0)
      INTEGER j,k
      REAL bma,bpa,fac,y,f(NMAX)
      DOUBLE PRECISION sum
      bma=0.5*(b-a)
      bpa=0.5*(b+a)
      do 11 k=1,n
        y=cos(PI*(k-0.5)/n)
        f(k)=func(y*bma+bpa)
11    continue
      fac=2./n
      do 13 j=1,n
        sum=0.d0
        do 12 k=1,n
          sum=sum+f(k)*cos((PI*(j-1))*((k-0.5d0)/n))
12      continue
        c(j)=fac*sum
13    continue
      return
      END
