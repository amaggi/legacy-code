      SUBROUTINE mprove(a,alud,n,np,indx,b,x)
      INTEGER n,np,indx(n),NMAX
      REAL a(np,np),alud(np,np),b(n),x(n)
      PARAMETER (NMAX=500)
CU    USES lubksb
      INTEGER i,j
      REAL r(NMAX)
      DOUBLE PRECISION sdp
      do 12 i=1,n
        sdp=-b(i)
        do 11 j=1,n
          sdp=sdp+dble(a(i,j))*dble(x(j))
11      continue
        r(i)=sdp
12    continue
      call lubksb(alud,n,np,indx,r)
      do 13 i=1,n
        x(i)=x(i)-r(i)
13    continue
      return
      END
