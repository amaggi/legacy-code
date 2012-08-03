      SUBROUTINE stiff(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
      INTEGER n,NMAX,MAXTRY
      REAL eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n),SAFETY,GROW,
     *PGROW,SHRNK,PSHRNK,ERRCON,GAM,A21,A31,A32,A2X,A3X,C21,C31,C32,C41,
     *C42,C43,B1,B2,B3,B4,E1,E2,E3,E4,C1X,C2X,C3X,C4X
      EXTERNAL derivs
      PARAMETER (NMAX=50,SAFETY=0.9,GROW=1.5,PGROW=-.25,SHRNK=0.5,
     *PSHRNK=-1./3.,ERRCON=.1296,MAXTRY=40)
      PARAMETER (GAM=1./2.,A21=2.,A31=48./25.,A32=6./25.,C21=-8.,
     *C31=372./25.,C32=12./5.,C41=-112./125.,C42=-54./125.,C43=-2./5.,
     *B1=19./9.,B2=1./2.,B3=25./108.,B4=125./108.,E1=17./54.,E2=7./36.,
     *E3=0.,E4=125./108.,C1X=1./2.,C2X=-3./2.,C3X=121./50.,C4X=29./250.,
     *A2X=1.,A3X=3./5.)
CU    USES derivs,jacobn,lubksb,ludcmp
      INTEGER i,j,jtry,indx(NMAX)
      REAL d,errmax,h,xsav,a(NMAX,NMAX),dfdx(NMAX),dfdy(NMAX,NMAX),
     *dysav(NMAX),err(NMAX),g1(NMAX),g2(NMAX),g3(NMAX),g4(NMAX),
     *ysav(NMAX)
      xsav=x
      do 11 i=1,n
        ysav(i)=y(i)
        dysav(i)=dydx(i)
11    continue
      call jacobn(xsav,ysav,dfdx,dfdy,n,NMAX)
      h=htry
      do 23 jtry=1,MAXTRY
        do 13 i=1,n
          do 12 j=1,n
            a(i,j)=-dfdy(i,j)
12        continue
          a(i,i)=1./(GAM*h)+a(i,i)
13      continue
        call ludcmp(a,n,NMAX,indx,d)
        do 14 i=1,n
          g1(i)=dysav(i)+h*C1X*dfdx(i)
14      continue
        call lubksb(a,n,NMAX,indx,g1)
        do 15 i=1,n
          y(i)=ysav(i)+A21*g1(i)
15      continue
        x=xsav+A2X*h
        call derivs(x,y,dydx)
        do 16 i=1,n
          g2(i)=dydx(i)+h*C2X*dfdx(i)+C21*g1(i)/h
16      continue
        call lubksb(a,n,NMAX,indx,g2)
        do 17 i=1,n
          y(i)=ysav(i)+A31*g1(i)+A32*g2(i)
17      continue
        x=xsav+A3X*h
        call derivs(x,y,dydx)
        do 18 i=1,n
          g3(i)=dydx(i)+h*C3X*dfdx(i)+(C31*g1(i)+C32*g2(i))/h
18      continue
        call lubksb(a,n,NMAX,indx,g3)
        do 19 i=1,n
          g4(i)=dydx(i)+h*C4X*dfdx(i)+(C41*g1(i)+C42*g2(i)+C43*g3(i))/h
19      continue
        call lubksb(a,n,NMAX,indx,g4)
        do 21 i=1,n
          y(i)=ysav(i)+B1*g1(i)+B2*g2(i)+B3*g3(i)+B4*g4(i)
          err(i)=E1*g1(i)+E2*g2(i)+E3*g3(i)+E4*g4(i)
21      continue
        x=xsav+h
        if(x.eq.xsav)pause 'stepsize not significant in stiff'
        errmax=0.
        do 22 i=1,n
          errmax=max(errmax,abs(err(i)/yscal(i)))
22      continue
        errmax=errmax/eps
        if(errmax.le.1.)then
          hdid=h
          if(errmax.gt.ERRCON)then
            hnext=SAFETY*h*errmax**PGROW
          else
            hnext=GROW*h
          endif
          return
        else
          hnext=SAFETY*h*errmax**PSHRNK
          if(hnext.lt.SHRNK*h)hnext=SHRNK*h
          h=hnext
        endif
23    continue
      pause 'exceeded MAXTRY in stiff'
      END
