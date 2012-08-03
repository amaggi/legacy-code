      SUBROUTINE dftcor(w,delta,a,b,endpts,corre,corim,corfac)
      REAL a,b,corfac,corim,corre,delta,w,endpts(8)
      REAL a0i,a0r,a1i,a1r,a2i,a2r,a3i,a3r,arg,c,cl,cr,s,sl,sr,t,t2,t4,
     *t6
      DOUBLE PRECISION cth,ctth,spth2,sth,sth4i,stth,th,th2,th4,tmth2,
     *tth4i
      th=w*delta
      if (a.ge.b.or.th.lt.0.d0.or.th.gt.3.1416d0)pause
     *'bad arguments to dftcor'
      if(abs(th).lt.5.d-2)then
        t=th
        t2=t*t
        t4=t2*t2
        t6=t4*t2
        corfac=1.-(11./720.)*t4+(23./15120.)*t6
        a0r=(-2./3.)+t2/45.+(103./15120.)*t4-(169./226800.)*t6
        a1r=(7./24.)-(7./180.)*t2+(5./3456.)*t4-(7./259200.)*t6
        a2r=(-1./6.)+t2/45.-(5./6048.)*t4+t6/64800.
        a3r=(1./24.)-t2/180.+(5./24192.)*t4-t6/259200.
        a0i=t*(2./45.+(2./105.)*t2-(8./2835.)*t4+(86./467775.)*t6)
        a1i=t*(7./72.-t2/168.+(11./72576.)*t4-(13./5987520.)*t6)
        a2i=t*(-7./90.+t2/210.-(11./90720.)*t4+(13./7484400.)*t6)
        a3i=t*(7./360.-t2/840.+(11./362880.)*t4-(13./29937600.)*t6)
      else
        cth=cos(th)
        sth=sin(th)
        ctth=cth**2-sth**2
        stth=2.d0*sth*cth
        th2=th*th
        th4=th2*th2
        tmth2=3.d0-th2
        spth2=6.d0+th2
        sth4i=1./(6.d0*th4)
        tth4i=2.d0*sth4i
        corfac=tth4i*spth2*(3.d0-4.d0*cth+ctth)
        a0r=sth4i*(-42.d0+5.d0*th2+spth2*(8.d0*cth-ctth))
        a0i=sth4i*(th*(-12.d0+6.d0*th2)+spth2*stth)
        a1r=sth4i*(14.d0*tmth2-7.d0*spth2*cth)
        a1i=sth4i*(30.d0*th-5.d0*spth2*sth)
        a2r=tth4i*(-4.d0*tmth2+2.d0*spth2*cth)
        a2i=tth4i*(-12.d0*th+2.d0*spth2*sth)
        a3r=sth4i*(2.d0*tmth2-spth2*cth)
        a3i=sth4i*(6.d0*th-spth2*sth)
      endif
      cl=a0r*endpts(1)+a1r*endpts(2)+a2r*endpts(3)+a3r*endpts(4)
      sl=a0i*endpts(1)+a1i*endpts(2)+a2i*endpts(3)+a3i*endpts(4)
      cr=a0r*endpts(8)+a1r*endpts(7)+a2r*endpts(6)+a3r*endpts(5)
      sr=-a0i*endpts(8)-a1i*endpts(7)-a2i*endpts(6)-a3i*endpts(5)
      arg=w*(b-a)
      c=cos(arg)
      s=sin(arg)
      corre=cl+c*cr-s*sr
      corim=sl+s*cr+c*sr
      return
      END
