      FUNCTION factln(n)
      INTEGER n
      REAL factln
CU    USES gammln
      REAL a(100),gammln
      SAVE a
      DATA a/100*-1./
      if (n.lt.0) pause 'negative factorial in factln'
      if (n.le.99) then
        if (a(n+1).lt.0.) a(n+1)=gammln(n+1.)
        factln=a(n+1)
      else
        factln=gammln(n+1.)
      endif
      return
      END
