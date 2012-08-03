      SUBROUTINE fpoly(x,p,np)
      INTEGER np
      REAL x,p(np)
      INTEGER j
      p(1)=1.
      do 11 j=2,np
        p(j)=p(j-1)*x
11    continue
      return
      END
