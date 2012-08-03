      SUBROUTINE fill0(u,n)
      INTEGER n
      DOUBLE PRECISION u(n,n)
      INTEGER i,j
      do 12 j=1,n
        do 11 i=1,n
          u(i,j)=0.d0
11      continue
12    continue
      return
      END
