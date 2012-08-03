      SUBROUTINE fred2(n,a,b,t,f,w,g,ak)
      INTEGER n,NMAX
      REAL a,b,f(n),t(n),w(n),g,ak
      EXTERNAL ak,g
      PARAMETER (NMAX=200)
CU    USES ak,g,gauleg,lubksb,ludcmp
      INTEGER i,j,indx(NMAX)
      REAL d,omk(NMAX,NMAX)
      if(n.gt.NMAX) pause 'increase NMAX in fred2'
      call gauleg(a,b,t,w,n)
      do 12 i=1,n
        do 11 j=1,n
          if(i.eq.j)then
            omk(i,j)=1.
          else
            omk(i,j)=0.
          endif
          omk(i,j)=omk(i,j)-ak(t(i),t(j))*w(j)
11      continue
        f(i)=g(t(i))
12    continue
      call ludcmp(omk,n,NMAX,indx,d)
      call lubksb(omk,n,NMAX,indx,f)
      return
      END
