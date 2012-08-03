      SUBROUTINE lop(out,u,n)
      INTEGER n
      DOUBLE PRECISION out(n,n),u(n,n)
      INTEGER i,j
      DOUBLE PRECISION h,h2i
      h=1.d0/(n-1)
      h2i=1.d0/(h*h)
      do 12 j=2,n-1
        do 11 i=2,n-1
          out(i,j)=h2i*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-4.d0*u(i,
     *j))+u(i,j)**2
11      continue
12    continue
      do 13 i=1,n
        out(i,1)=0.d0
        out(i,n)=0.d0
        out(1,i)=0.d0
        out(n,i)=0.d0
13    continue
      return
      END
