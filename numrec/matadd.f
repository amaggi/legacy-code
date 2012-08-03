      SUBROUTINE matadd(a,b,c,n)
      INTEGER n
      DOUBLE PRECISION a(n,n),b(n,n),c(n,n)
      INTEGER i,j
      do 12 j=1,n
        do 11 i=1,n
          c(i,j)=a(i,j)+b(i,j)
11      continue
12    continue
      return
      END
