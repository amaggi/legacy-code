      SUBROUTINE mmid(y,dydx,nvar,xs,htot,nstep,yout,derivs)
      INTEGER nstep,nvar,NMAX
      REAL htot,xs,dydx(nvar),y(nvar),yout(nvar)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
      INTEGER i,n
      REAL h,h2,swap,x,ym(NMAX),yn(NMAX)
      h=htot/nstep
      do 11 i=1,nvar
        ym(i)=y(i)
        yn(i)=y(i)+h*dydx(i)
11    continue
      x=xs+h
      call derivs(x,yn,yout)
      h2=2.*h
      do 13 n=2,nstep
        do 12 i=1,nvar
          swap=ym(i)+h2*yout(i)
          ym(i)=yn(i)
          yn(i)=swap
12      continue
        x=x+h
        call derivs(x,yn,yout)
13    continue
      do 14 i=1,nvar
        yout(i)=0.5*(ym(i)+yn(i)+h*yout(i))
14    continue
      return
      END
