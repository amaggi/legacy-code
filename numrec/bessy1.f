      FUNCTION bessy1(x)
      REAL bessy1,x
CU    USES bessj1
      REAL xx,z,bessj1
      DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,
     *s1,s2,s3,s4,s5,s6,s7,y
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,
     *s5,s6,s7
      DATA p1,p2,p3,p4,p5/1.d0,.183105d-2,-.3516396496d-4,
     *.2457520174d-5,-.240337019d-6/, q1,q2,q3,q4,q5/.04687499995d0,
     *-.2002690873d-3,.8449199096d-5,-.88228987d-6,.105787412d-6/
      DATA r1,r2,r3,r4,r5,r6/-.4900604943d13,.1275274390d13,
     *-.5153438139d11,.7349264551d9,-.4237922726d7,.8511937935d4/,s1,s2,
     *s3,s4,s5,s6,s7/.2499580570d14,.4244419664d12,.3733650367d10,
     *.2245904002d8,.1020426050d6,.3549632885d3,1.d0/
      if(x.lt.8.)then
        y=x**2
        bessy1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+
     *y*(s4+y*(s5+y*(s6+y*s7))))))+.636619772*(bessj1(x)*log(x)-1./x)
      else
        z=8./x
        y=z**2
        xx=x-2.356194491
        bessy1=sqrt(.636619772/x)*(sin(xx)*(p1+y*(p2+y*(p3+y*(p4+y*
     *p5))))+z*cos(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
      endif
      return
      END
