      DOUBLE PRECISION FUNCTION anorm2(a,n)
      INTEGER n
      DOUBLE PRECISION a(n,n)
      INTEGER i,j
      DOUBLE PRECISION sum
      sum=0.d0
      do 12 j=1,n
        do 11 i=1,n
          sum=sum+a(i,j)**2
11      continue
12    continue
      anorm2=sqrt(sum)/n
      return
      END
