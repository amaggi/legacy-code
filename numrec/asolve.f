      SUBROUTINE asolve(n,b,x,itrnsp)
      INTEGER n,itrnsp,ija,NMAX,i
      DOUBLE PRECISION x(n),b(n),sa
      PARAMETER (NMAX=1000)
      COMMON /mat/ sa(NMAX),ija(NMAX)
      do 11 i=1,n
        x(i)=b(i)/sa(i)
11    continue
      return
      END
