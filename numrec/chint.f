      SUBROUTINE chint(a,b,c,cint,n)
      INTEGER n
      REAL a,b,c(n),cint(n)
      INTEGER j
      REAL con,fac,sum
      con=0.25*(b-a)
      sum=0.
      fac=1.
      do 11 j=2,n-1
        cint(j)=con*(c(j-1)-c(j+1))/(j-1)
        sum=sum+fac*cint(j)
        fac=-fac
11    continue
      cint(n)=con*c(n-1)/(n-1)
      sum=sum+fac*cint(n)
      cint(1)=2.*sum
      return
      END
