      PROGRAM fredex
      INTEGER NMAX
      REAL PI
      PARAMETER (NMAX=100,PI=3.14159265)
      INTEGER indx(NMAX),j,n
      REAL a(NMAX,NMAX),g(NMAX),x,d
CU    USES quadmx,ludcmp,lubksb
      n=40
      call quadmx(a,n,NMAX)
      call ludcmp(a,n,NMAX,indx,d)
      do 11 j=1,n
        x=(j-1)*PI/(n-1)
        g(j)=sin(x)
11    continue
      call lubksb(a,n,NMAX,indx,g)
      do 12 j=1,n
        x=(j-1)*PI/(n-1)
        write (*,*) j,x,g(j)
12    continue
      write (*,*) 'normal completion'
      END
