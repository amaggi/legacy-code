      SUBROUTINE ks2d2s(x1,y1,n1,x2,y2,n2,d,prob)
      INTEGER n1,n2
      REAL d,prob,x1(n1),x2(n2),y1(n1),y2(n2)
CU    USES pearsn,probks,quadct
      INTEGER j
      REAL d1,d2,dum,dumm,fa,fb,fc,fd,ga,gb,gc,gd,r1,r2,rr,sqen,probks
      d1=0.0
      do 11 j=1,n1
        call quadct(x1(j),y1(j),x1,y1,n1,fa,fb,fc,fd)
        call quadct(x1(j),y1(j),x2,y2,n2,ga,gb,gc,gd)
        d1=max(d1,abs(fa-ga),abs(fb-gb),abs(fc-gc),abs(fd-gd))
11    continue
      d2=0.0
      do 12 j=1,n2
        call quadct(x2(j),y2(j),x1,y1,n1,fa,fb,fc,fd)
        call quadct(x2(j),y2(j),x2,y2,n2,ga,gb,gc,gd)
        d2=max(d2,abs(fa-ga),abs(fb-gb),abs(fc-gc),abs(fd-gd))
12    continue
      d=0.5*(d1+d2)
      sqen=sqrt(float(n1)*float(n2)/float(n1+n2))
      call pearsn(x1,y1,n1,r1,dum,dumm)
      call pearsn(x2,y2,n2,r2,dum,dumm)
      rr=sqrt(1.0-0.5*(r1**2+r2**2))
      prob=probks(d*sqen/(1.0+rr*(0.25-0.75/sqen)))
      return
      END
