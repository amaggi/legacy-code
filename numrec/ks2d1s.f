      SUBROUTINE ks2d1s(x1,y1,n1,quadvl,d1,prob)
      INTEGER n1
      REAL d1,prob,x1(n1),y1(n1)
      EXTERNAL quadvl
CU    USES pearsn,probks,quadct,quadvl
      INTEGER j
      REAL dum,dumm,fa,fb,fc,fd,ga,gb,gc,gd,r1,rr,sqen,probks
      d1=0.0
      do 11 j=1,n1
        call quadct(x1(j),y1(j),x1,y1,n1,fa,fb,fc,fd)
        call quadvl(x1(j),y1(j),ga,gb,gc,gd)
        d1=max(d1,abs(fa-ga),abs(fb-gb),abs(fc-gc),abs(fd-gd))
11    continue
      call pearsn(x1,y1,n1,r1,dum,dumm)
      sqen=sqrt(float(n1))
      rr=sqrt(1.0-r1**2)
      prob=probks(d1*sqen/(1.0+rr*(0.25-0.75/sqen)))
      return
      END
