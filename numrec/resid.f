      SUBROUTINE resid(res,u,rhs,n)
      INTEGER n
      DOUBLE PRECISION res(n,n),rhs(n,n),u(n,n)
      INTEGER i,j
      DOUBLE PRECISION h,h2i
      h=1.d0/(n-1)
      h2i=1.d0/(h*h)
      do 12 j=2,n-1
        do 11 i=2,n-1
          res(i,j)=-h2i*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-4.d0*u(i,
     *j))+rhs(i,j)
11      continue
12    continue
      do 13 i=1,n
        res(i,1)=0.d0
        res(i,n)=0.d0
        res(1,i)=0.d0
        res(n,i)=0.d0
13    continue
      return
      END
