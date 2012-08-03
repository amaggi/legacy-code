      SUBROUTINE rsolv(a,n,np,d,b)
      INTEGER n,np
      REAL a(np,np),b(n),d(n)
      INTEGER i,j
      REAL sum
      b(n)=b(n)/d(n)
      do 12 i=n-1,1,-1
        sum=0.
        do 11 j=i+1,n
          sum=sum+a(i,j)*b(j)
11      continue
        b(i)=(b(i)-sum)/d(i)
12    continue
      return
      END
