      SUBROUTINE midsql(funk,aa,bb,s,n)
      INTEGER n
      REAL aa,bb,s,funk
      EXTERNAL funk
      INTEGER it,j
      REAL ddel,del,sum,tnm,x,func,a,b
      func(x)=2.*x*funk(aa+x**2)
      b=sqrt(bb-aa)
      a=0.
      if (n.eq.1) then
        s=(b-a)*func(0.5*(a+b))
      else
        it=3**(n-2)
        tnm=it
        del=(b-a)/(3.*tnm)
        ddel=del+del
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+ddel
          sum=sum+func(x)
          x=x+del
11      continue
        s=(s+(b-a)*sum/tnm)/3.
      endif
      return
      END
