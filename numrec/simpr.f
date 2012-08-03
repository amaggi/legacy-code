      SUBROUTINE simpr(y,dydx,dfdx,dfdy,nmax,n,xs,htot,nstep,yout,
     *derivs)
      INTEGER n,nmax,nstep,NMAXX
      REAL htot,xs,dfdx(n),dfdy(nmax,nmax),dydx(n),y(n),yout(n)
      EXTERNAL derivs
      PARAMETER (NMAXX=50)
CU    USES derivs,lubksb,ludcmp
      INTEGER i,j,nn,indx(NMAXX)
      REAL d,h,x,a(NMAXX,NMAXX),del(NMAXX),ytemp(NMAXX)
      h=htot/nstep
      do 12 i=1,n
        do 11 j=1,n
          a(i,j)=-h*dfdy(i,j)
11      continue
        a(i,i)=a(i,i)+1.
12    continue
      call ludcmp(a,n,NMAXX,indx,d)
      do 13 i=1,n
        yout(i)=h*(dydx(i)+h*dfdx(i))
13    continue
      call lubksb(a,n,NMAXX,indx,yout)
      do 14 i=1,n
        del(i)=yout(i)
        ytemp(i)=y(i)+del(i)
14    continue
      x=xs+h
      call derivs(x,ytemp,yout)
      do 17 nn=2,nstep
        do 15 i=1,n
          yout(i)=h*yout(i)-del(i)
15      continue
        call lubksb(a,n,NMAXX,indx,yout)
        do 16 i=1,n
          del(i)=del(i)+2.*yout(i)
          ytemp(i)=ytemp(i)+del(i)
16      continue
        x=x+h
        call derivs(x,ytemp,yout)
17    continue
      do 18 i=1,n
        yout(i)=h*yout(i)-del(i)
18    continue
      call lubksb(a,n,NMAXX,indx,yout)
      do 19 i=1,n
        yout(i)=ytemp(i)+yout(i)
19    continue
      return
      END
