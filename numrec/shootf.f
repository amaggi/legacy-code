CU    SUBROUTINE shootf(n,v,f) is named "funcv" for use with "newt"
      SUBROUTINE funcv(n,v,f)
      INTEGER n,nvar,nn2,kmax,kount,KMAXX,NMAX
      REAL f(n),v(n),x1,x2,xf,dxsav,xp,yp,EPS
      PARAMETER (NMAX=50,KMAXX=200,EPS=1.e-6)
      COMMON /caller/ x1,x2,xf,nvar,nn2
      COMMON /path/ kmax,kount,dxsav,xp(KMAXX),yp(NMAX,KMAXX)
CU    USES derivs,load1,load2,odeint,rkqs,score
      INTEGER i,nbad,nok
      REAL h1,hmin,f1(NMAX),f2(NMAX),y(NMAX)
      EXTERNAL derivs,rkqs
      kmax=0
      h1=(x2-x1)/100.
      hmin=0.
      call load1(x1,v,y)
      call odeint(y,nvar,x1,xf,EPS,h1,hmin,nok,nbad,derivs,rkqs)
      call score(xf,y,f1)
      call load2(x2,v(nn2+1),y)
      call odeint(y,nvar,x2,xf,EPS,h1,hmin,nok,nbad,derivs,rkqs)
      call score(xf,y,f2)
      do 11 i=1,n
        f(i)=f1(i)-f2(i)
11    continue
      return
      END
