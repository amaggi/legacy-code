      SUBROUTINE polcoe(x,y,n,cof)
      INTEGER n,NMAX
      REAL cof(n),x(n),y(n)
      PARAMETER (NMAX=15)
      INTEGER i,j,k
      REAL b,ff,phi,s(NMAX)
      do 11 i=1,n
        s(i)=0.
        cof(i)=0.
11    continue
      s(n)=-x(1)
      do 13 i=2,n
        do 12 j=n+1-i,n-1
          s(j)=s(j)-x(i)*s(j+1)
12      continue
        s(n)=s(n)-x(i)
13    continue
      do 16 j=1,n
        phi=n
        do 14 k=n-1,1,-1
          phi=k*s(k+1)+x(j)*phi
14      continue
        ff=y(j)/phi
        b=1.
        do 15 k=n,1,-1
          cof(k)=cof(k)+b*ff
          b=s(k)+x(j)*b
15      continue
16    continue
      return
      END
