      FUNCTION dawson(x)
      INTEGER NMAX
      REAL dawson,x,H,A1,A2,A3
      PARAMETER (NMAX=6,H=0.4,A1=2./3.,A2=0.4,A3=2./7.)
      INTEGER i,init,n0
      REAL d1,d2,e1,e2,sum,x2,xp,xx,c(NMAX)
      SAVE init,c
      DATA init/0/
      if(init.eq.0)then
        init=1
        do 11 i=1,NMAX
          c(i)=exp(-((2.*float(i)-1.)*H)**2)
11      continue
      endif
      if(abs(x).lt.0.2)then
        x2=x**2
        dawson=x*(1.-A1*x2*(1.-A2*x2*(1.-A3*x2)))
      else
        xx=abs(x)
        n0=2*nint(0.5*xx/H)
        xp=xx-float(n0)*H
        e1=exp(2.*xp*H)
        e2=e1**2
        d1=float(n0+1)
        d2=d1-2.
        sum=0.
        do 12 i=1,NMAX
          sum=sum+c(i)*(e1/d1+1./(d2*e1))
          d1=d1+2.
          d2=d2-2.
          e1=e2*e1
12      continue
        dawson=0.5641895835*sign(exp(-xp**2),x)*sum
      endif
      return
      END
