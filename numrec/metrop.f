      SUBROUTINE metrop(de,t,ans)
      REAL de,t
      LOGICAL ans
CU    USES ran3
      INTEGER jdum
      REAL ran3
      SAVE jdum
      DATA jdum /1/
      ans=(de.lt.0.0).or.(ran3(jdum).lt.exp(-de/t))
      return
      END
