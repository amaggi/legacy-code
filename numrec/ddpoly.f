      SUBROUTINE ddpoly(c,nc,x,pd,nd)
      INTEGER nc,nd
      REAL x,c(nc),pd(nd)
      INTEGER i,j,nnd
      REAL const
      pd(1)=c(nc)
      do 11 j=2,nd
        pd(j)=0.
11    continue
      do 13 i=nc-1,1,-1
        nnd=min(nd,nc+1-i)
        do 12 j=nnd,2,-1
          pd(j)=pd(j)*x+pd(j-1)
12      continue
        pd(1)=pd(1)*x+c(i)
13    continue
      const=2.
      do 14 i=3,nd
        pd(i)=const*pd(i)
        const=const*i
14    continue
      return
      END
