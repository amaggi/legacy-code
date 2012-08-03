      SUBROUTINE relax2(u,rhs,n)
      INTEGER n
      DOUBLE PRECISION rhs(n,n),u(n,n)
      INTEGER i,ipass,isw,j,jsw
      DOUBLE PRECISION foh2,h,h2i,res
      h=1.d0/(n-1)
      h2i=1.d0/(h*h)
      foh2=-4.d0*h2i
      jsw=1
      do 13 ipass=1,2
        isw=jsw
        do 12 j=2,n-1
          do 11 i=isw+1,n-1,2
            res=h2i*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-4.d0*u(i,j))+
     *u(i,j)**2-rhs(i,j)
            u(i,j)=u(i,j)-res/(foh2+2.d0*u(i,j))
11        continue
          isw=3-isw
12      continue
        jsw=3-jsw
13    continue
      return
      END
