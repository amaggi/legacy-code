      FUNCTION erfc(x)
      REAL erfc,x
CU    USES gammp,gammq
      REAL gammp,gammq
      if(x.lt.0.)then
        erfc=1.+gammp(.5,x**2)
      else
        erfc=gammq(.5,x**2)
      endif
      return
      END
