      FUNCTION factrl(n)
      INTEGER n
      REAL factrl
CU    USES gammln
      INTEGER j,ntop
      REAL a(33),gammln
      SAVE ntop,a
      DATA ntop,a(1)/0,1./
      if (n.lt.0) then
        pause 'negative factorial in factrl'
      else if (n.le.ntop) then
        factrl=a(n+1)
      else if (n.le.32) then
        do 11 j=ntop+1,n
          a(j+1)=j*a(j)
11      continue
        ntop=n
        factrl=a(n+1)
      else
        factrl=exp(gammln(n+1.))
      endif
      return
      END
