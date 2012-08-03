      SUBROUTINE relax(u,rhs,n)
      INTEGER n
      DOUBLE PRECISION rhs(n,n),u(n,n)
      INTEGER i,ipass,isw,j,jsw
      DOUBLE PRECISION h,h2
      h=1.d0/(n-1)
      h2=h*h
      jsw=1
      do 13 ipass=1,2
        isw=jsw
        do 12 j=2,n-1
          do 11 i=isw+1,n-1,2
            u(i,j)=0.25d0*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-h2*rhs(i,
     *j))
11        continue
          isw=3-isw
12      continue
        jsw=3-jsw
13    continue
      return
      END
