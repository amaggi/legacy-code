      SUBROUTINE quadmx(a,n,np)
      INTEGER n,np,NMAX
      REAL a(np,np),PI
      DOUBLE PRECISION xx
      PARAMETER (PI=3.14159265,NMAX=257)
      COMMON /momcom/ xx
      EXTERNAL kermom
CU    USES wwghts,kermom
      INTEGER j,k
      REAL h,wt(NMAX),x,cx,y
      h=PI/(n-1)
      do 12 j=1,n
        x=(j-1)*h
        xx=x
        call wwghts(wt,n,h,kermom)
        cx=cos(x)
        do 11 k=1,n
          y=(k-1)*h
          a(j,k)=wt(k)*cx*cos(y)
11      continue
        a(j,j)=a(j,j)+1.
12    continue
      return
      END
