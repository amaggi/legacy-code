      SUBROUTINE rzextr(iest,xest,yest,yz,dy,nv)
      INTEGER iest,nv,IMAX,NMAX
      REAL xest,dy(nv),yest(nv),yz(nv)
      PARAMETER (IMAX=13,NMAX=50)
      INTEGER j,k
      REAL b,b1,c,ddy,v,yy,d(NMAX,IMAX),fx(IMAX),x(IMAX)
      SAVE d,x
      x(iest)=xest
      if(iest.eq.1) then
        do 11 j=1,nv
          yz(j)=yest(j)
          d(j,1)=yest(j)
          dy(j)=yest(j)
11      continue
      else
        do 12 k=1,iest-1
          fx(k+1)=x(iest-k)/xest
12      continue
        do 14 j=1,nv
          yy=yest(j)
          v=d(j,1)
          c=yy
          d(j,1)=yy
          do 13 k=2,iest
            b1=fx(k)*v
            b=b1-c
            if(b.ne.0.) then
              b=(c-v)/b
              ddy=c*b
              c=b1*b
            else
              ddy=v
            endif
            if (k.ne.iest) v=d(j,k)
            d(j,k)=ddy
            yy=yy+ddy
13        continue
          dy(j)=ddy
          yz(j)=yy
14      continue
      endif
      return
      END
