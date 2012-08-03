      SUBROUTINE slvsml(u,rhs)
      DOUBLE PRECISION rhs(3,3),u(3,3)
CU    USES fill0
      DOUBLE PRECISION h
      call fill0(u,3)
      h=.5d0
      u(2,2)=-h*h*rhs(2,2)/4.d0
      return
      END
