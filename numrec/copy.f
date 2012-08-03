      SUBROUTINE copy(aout,ain,n)
      INTEGER n
      DOUBLE PRECISION ain(n,n),aout(n,n)
      INTEGER i,j
      do 12 i=1,n
        do 11 j=1,n
          aout(j,i)=ain(j,i)
11      continue
12    continue
      return
      END
