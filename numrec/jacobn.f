      SUBROUTINE jacobn(x,y,dfdx,dfdy,n,nmax)
      INTEGER n,nmax,i
      REAL x,y(*),dfdx(*),dfdy(nmax,nmax)
      do 11 i=1,3
        dfdx(i)=0.
11    continue
      dfdy(1,1)=-.013-1000.*y(3)
      dfdy(1,2)=0.
      dfdy(1,3)=-1000.*y(1)
      dfdy(2,1)=0.
      dfdy(2,2)=-2500.*y(3)
      dfdy(2,3)=-2500.*y(2)
      dfdy(3,1)=-.013-1000.*y(3)
      dfdy(3,2)=-2500.*y(3)
      dfdy(3,3)=-1000.*y(1)-2500.*y(2)
      return
      END

      SUBROUTINE derivs(x,y,dydx)
      REAL x,y(*),dydx(*)
      dydx(1)=-.013*y(1)-1000.*y(1)*y(3)
      dydx(2)=-2500.*y(2)*y(3)
      dydx(3)=-.013*y(1)-1000.*y(1)*y(3)-2500.*y(2)*y(3)
      return
      END

