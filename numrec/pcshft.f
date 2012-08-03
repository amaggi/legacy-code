      SUBROUTINE pcshft(a,b,d,n)
      INTEGER n
      REAL a,b,d(n)
      INTEGER j,k
      REAL const,fac
      const=2./(b-a)
      fac=const
      do 11 j=2,n
        d(j)=d(j)*fac
        fac=fac*const
11    continue
      const=0.5*(a+b)
      do 13 j=1,n-1
        do 12 k=n-1,j,-1
          d(k)=d(k)-const*d(k+1)
12      continue
13    continue
      return
      END
