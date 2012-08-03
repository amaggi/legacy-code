CU    SUBROUTINE shoot(n2,v,f) is named "funcv" for use with "newt"
      SUBROUTINE funcv(n2,v,f)
      INTEGER n2,nvar,kmax,kount,KMAXX,NMAX
      REAL f(n2),v(n2),x1,x2,dxsav,xp,yp,EPS
      PARAMETER (NMAX=50,KMAXX=200,EPS=1.e-6)
      COMMON /caller/ x1,x2,nvar
      COMMON /path/ kmax,kount,dxsav,xp(KMAXX),yp(NMAX,KMAXX)
CU    USES derivs,load,odeint,rkqs,score
      INTEGER nbad,nok
      REAL h1,hmin,y(NMAX)
      EXTERNAL derivs,rkqs
      kmax=0
      h1=(x2-x1)/100.
      hmin=0.
      call load(x1,v,y)
      call odeint(y,nvar,x1,x2,EPS,h1,hmin,nok,nbad,derivs,rkqs)
      call score(x2,y,f)
      return
      END
