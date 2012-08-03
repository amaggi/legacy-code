      SUBROUTINE gaucof(n,a,b,amu0,x,w)
      INTEGER n,NMAX
      REAL amu0,a(n),b(n),w(n),x(n)
      PARAMETER (NMAX=64)
CU    USES eigsrt,tqli
      INTEGER i,j
      REAL z(NMAX,NMAX)
      do 12 i=1,n
        if(i.ne.1)b(i)=sqrt(b(i))
        do 11 j=1,n
          if(i.eq.j)then
            z(i,j)=1.
          else
            z(i,j)=0.
          endif
11      continue
12    continue
      call tqli(a,b,n,NMAX,z)
      call eigsrt(a,z,n,NMAX)
      do 13 i=1,n
        x(i)=a(i)
        w(i)=amu0*z(1,i)**2
13    continue
      return
      END
