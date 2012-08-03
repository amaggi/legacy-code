      SUBROUTINE fixrts(d,m)
      INTEGER m,MMAX
      REAL d(m)
      PARAMETER (MMAX=100)
CU    USES zroots
      INTEGER i,j
      LOGICAL polish
      COMPLEX a(MMAX),roots(MMAX)
      a(m+1)=cmplx(1.,0.)
      do 11 j=m,1,-1
        a(j)=cmplx(-d(m+1-j),0.)
11    continue
      polish=.true.
      call zroots(a,m,roots,polish)
      do 12 j=1,m
        if(abs(roots(j)).gt.1.)then
          roots(j)=1./conjg(roots(j))
        endif
12    continue
      a(1)=-roots(1)
      a(2)=cmplx(1.,0.)
      do 14 j=2,m
        a(j+1)=cmplx(1.,0.)
        do 13 i=j,2,-1
          a(i)=a(i-1)-roots(j)*a(i)
13      continue
        a(1)=-roots(j)*a(1)
14    continue
      do 15 j=1,m
        d(m+1-j)=-real(a(j))
15    continue
      return
      END
