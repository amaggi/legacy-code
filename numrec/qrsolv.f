      SUBROUTINE qrsolv(a,n,np,c,d,b)
      INTEGER n,np
      REAL a(np,np),b(n),c(n),d(n)
CU    USES rsolv
      INTEGER i,j
      REAL sum,tau
      do 13 j=1,n-1
        sum=0.
        do 11 i=j,n
          sum=sum+a(i,j)*b(i)
11      continue
        tau=sum/c(j)
        do 12 i=j,n
          b(i)=b(i)-tau*a(i,j)
12      continue
13    continue
      call rsolv(a,n,np,d,b)
      return
      END
