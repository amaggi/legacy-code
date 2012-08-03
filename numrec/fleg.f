      SUBROUTINE fleg(x,pl,nl)
      INTEGER nl
      REAL x,pl(nl)
      INTEGER j
      REAL d,f1,f2,twox
      pl(1)=1.
      pl(2)=x
      if(nl.gt.2) then
        twox=2.*x
        f2=x
        d=1.
        do 11 j=3,nl
          f1=d
          f2=f2+twox
          d=d+1.
          pl(j)=(f2*pl(j-1)-f1*pl(j-2))/d
11      continue
      endif
      return
      END
