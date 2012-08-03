      FUNCTION fmin(x)
      INTEGER n,NP
      REAL fmin,x(*),fvec
      PARAMETER (NP=40)
      COMMON /newtv/ fvec(NP),n
      SAVE /newtv/
CU    USES funcv
      INTEGER i
      REAL sum
      call funcv(n,x,fvec)
      sum=0.
      do 11 i=1,n
        sum=sum+fvec(i)**2
11    continue
      fmin=0.5*sum
      return
      END
