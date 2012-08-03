      SUBROUTINE slvsm2(u,rhs)
      DOUBLE PRECISION rhs(3,3),u(3,3)
CU    USES fill0
      DOUBLE PRECISION disc,fact,h
      call fill0(u,3)
      h=.5d0
      fact=2./h**2
      disc=sqrt(fact**2+rhs(2,2))
      u(2,2)=-rhs(2,2)/(fact+disc)
      return
      END
