      SUBROUTINE chebpc(c,d,n)
      INTEGER n,NMAX
      REAL c(n),d(n)
      PARAMETER (NMAX=50)
      INTEGER j,k
      REAL sv,dd(NMAX)
      do 11 j=1,n
        d(j)=0.
        dd(j)=0.
11    continue
      d(1)=c(n)
      do 13 j=n-1,2,-1
        do 12 k=n-j+1,2,-1
          sv=d(k)
          d(k)=2.*d(k-1)-dd(k)
          dd(k)=sv
12      continue
        sv=d(1)
        d(1)=-dd(1)+c(j)
        dd(1)=sv
13    continue
      do 14 j=n,2,-1
        d(j)=d(j-1)-dd(j)
14    continue
      d(1)=-dd(1)+0.5*c(1)
      return
      END
