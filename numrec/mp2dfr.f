      SUBROUTINE mp2dfr(a,s,n,m)
      INTEGER m,n,IAZ
      CHARACTER*1 a(*),s(*)
      PARAMETER (IAZ=48)
CU    USES mplsh,mpsmu
      INTEGER j
        m=2.408*n
        do 11 j=1,m
          call mpsmu(a,a,n,10)
          s(j)=char(ichar(a(1))+IAZ)
          call mplsh(a,n)
11      continue
      return
      END
