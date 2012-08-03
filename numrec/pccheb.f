      SUBROUTINE pccheb(d,c,n)
      INTEGER n
      REAL c(n),d(n)
      INTEGER j,jm,jp,k
      REAL fac,pow
      pow=1.
      c(1)=2.*d(1)
      do 12 k=2,n
        c(k)=0.
        fac=d(k)/pow
        jm=k-1
        jp=1
        do 11 j=k,1,-2
          c(j)=c(j)+fac
          fac=fac*float(jm)/float(jp)
          jm=jm-1
          jp=jp+1
11      continue
        pow=2.*pow
12    continue
      return
      END
