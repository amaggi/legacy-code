      SUBROUTINE polcof(xa,ya,n,cof)
      INTEGER n,NMAX
      REAL cof(n),xa(n),ya(n)
      PARAMETER (NMAX=15)
CU    USES polint
      INTEGER i,j,k
      REAL dy,xmin,x(NMAX),y(NMAX)
      do 11 j=1,n
        x(j)=xa(j)
        y(j)=ya(j)
11    continue
      do 14 j=1,n
        call polint(x,y,n+1-j,0.,cof(j),dy)
        xmin=1.e38
        k=0
        do 12 i=1,n+1-j
          if (abs(x(i)).lt.xmin)then
            xmin=abs(x(i))
            k=i
          endif
          if(x(i).ne.0.)y(i)=(y(i)-cof(j))/x(i)
12      continue
        do 13 i=k+1,n+1-j
          y(i-1)=y(i)
          x(i-1)=x(i)
13      continue
14    continue
      return
      END
