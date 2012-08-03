      SUBROUTINE stoerm(y,d2y,nv,xs,htot,nstep,yout,derivs)
      INTEGER nstep,nv,NMAX
      REAL htot,xs,d2y(nv),y(nv),yout(nv)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
CU    USES derivs
      INTEGER i,n,neqns,nn
      REAL h,h2,halfh,x,ytemp(NMAX)
      h=htot/nstep
      halfh=0.5*h
      neqns=nv/2
      do 11 i=1,neqns
        n=neqns+i
        ytemp(n)=h*(y(n)+halfh*d2y(i))
        ytemp(i)=y(i)+ytemp(n)
11    continue
      x=xs+h
      call derivs(x,ytemp,yout)
      h2=h*h
      do 13 nn=2,nstep
        do 12 i=1,neqns
          n=neqns+i
          ytemp(n)=ytemp(n)+h2*yout(i)
          ytemp(i)=ytemp(i)+ytemp(n)
12      continue
        x=x+h
        call derivs(x,ytemp,yout)
13    continue
      do 14 i=1,neqns
        n=neqns+i
        yout(n)=ytemp(n)/h+halfh*yout(i)
        yout(i)=ytemp(i)
14    continue
      return
      END
