      SUBROUTINE atimes(n,x,r,itrnsp)
      INTEGER n,itrnsp,ija,NMAX
      DOUBLE PRECISION x(n),r(n),sa
      PARAMETER (NMAX=1000)
      COMMON /mat/ sa(NMAX),ija(NMAX)
CU    USES dsprsax,dsprstx
      if (itrnsp.eq.0) then
        call dsprsax(sa,ija,x,r,n)
      else
        call dsprstx(sa,ija,x,r,n)
      endif
      return
      END
