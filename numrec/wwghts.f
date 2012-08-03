      SUBROUTINE wwghts(wghts,n,h,kermom)
      INTEGER n
      REAL wghts(n),h
      EXTERNAL kermom
CU    USES kermom
      INTEGER j,k
      DOUBLE PRECISION wold(4),wnew(4),w(4),hh,hi,c,fac,a,b
      hh=h
      hi=1.d0/hh
      do 11 j=1,n
        wghts(j)=0.
11    continue
      call kermom(wold,0.d0,4)
      if (n.ge.4) then
        b=0.d0
        do 14 j=1,n-3
          c=j-1
          a=b
          b=a+hh
          if (j.eq.n-3) b=(n-1)*hh
          call kermom(wnew,b,4)
          fac=1.d0
          do 12 k=1,4
            w(k)=(wnew(k)-wold(k))*fac
            fac=fac*hi
12        continue
          wghts(j)=wghts(j)+((c+1.d0)*(c+2.d0)*(c+3.d0)*w(1)-(11.d0+c*
     *(12.d0+c*3.d0))*w(2)+3.d0*(c+2.d0)*w(3)-w(4))/6.d0
          wghts(j+1)=wghts(j+1)+(-c*(c+2.d0)*(c+3.d0)*w(1)+(6.d0+c*
     *(10.d0+c*3.d0))*w(2)-(3.d0*c+5.d0)*w(3)+w(4))*.5d0
          wghts(j+2)=wghts(j+2)+(c*(c+1.d0)*(c+3.d0)*w(1)-(3.d0+c*(8.d0+
     *c*3.d0))*w(2)+(3.d0*c+4.d0)*w(3)-w(4))*.5d0
          wghts(j+3)=wghts(j+3)+(-c*(c+1.d0)*(c+2.d0)*w(1)+(2.d0+c*
     *(6.d0+c*3.d0))*w(2)-3.d0*(c+1.d0)*w(3)+w(4))/6.d0
          do 13 k=1,4
            wold(k)=wnew(k)
13        continue
14      continue
      else if (n.eq.3) then
        call kermom(wnew,hh+hh,3)
        w(1)=wnew(1)-wold(1)
        w(2)=hi*(wnew(2)-wold(2))
        w(3)=hi**2*(wnew(3)-wold(3))
        wghts(1)=w(1)-1.5d0*w(2)+0.5d0*w(3)
        wghts(2)=2.d0*w(2)-w(3)
        wghts(3)=0.5d0*(w(3)-w(2))
      else if (n.eq.2) then
        call kermom(wnew,hh,2)
        wghts(2)=hi*(wnew(2)-wold(2))
        wghts(1)=wnew(1)-wold(1)-wghts(2)
      endif
      END
