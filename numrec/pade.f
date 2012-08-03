      SUBROUTINE pade(cof,n,resid)
      INTEGER n,NMAX
      REAL resid,BIG
      DOUBLE PRECISION cof(2*n+1)
      PARAMETER (NMAX=20,BIG=1.E30)
CU    USES lubksb,ludcmp,mprove
      INTEGER j,k,indx(NMAX)
      REAL d,rr,rrold,sum,q(NMAX,NMAX),qlu(NMAX,NMAX),x(NMAX),y(NMAX),
     *z(NMAX)
      do 12 j=1,n
        x(j)=cof(n+j+1)
        y(j)=x(j)
        do 11 k=1,n
          q(j,k)=cof(j-k+n+1)
          qlu(j,k)=q(j,k)
11      continue
12    continue
      call ludcmp(qlu,n,NMAX,indx,d)
      call lubksb(qlu,n,NMAX,indx,x)
      rr=BIG
1     continue
        rrold=rr
        do 13 j=1,n
          z(j)=x(j)
13      continue
        call mprove(q,qlu,n,NMAX,indx,y,x)
        rr=0.
        do 14 j=1,n
          rr=rr+(z(j)-x(j))**2
14      continue
      if(rr.lt.rrold)goto 1
      resid=sqrt(rr)
      do 16 k=1,n
        sum=cof(k+1)
        do 15 j=1,k
          sum=sum-x(j)*cof(k-j+1)
15      continue
        y(k)=sum
16    continue
      do 17 j=1,n
        cof(j+1)=y(j)
        cof(j+n+1)=-x(j)
17    continue
      return
      END
