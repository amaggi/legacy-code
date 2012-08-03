      FUNCTION df1dim(x)
      INTEGER NMAX
      REAL df1dim,x
      PARAMETER (NMAX=50)
CU    USES dfunc
      INTEGER j,ncom
      REAL df(NMAX),pcom(NMAX),xicom(NMAX),xt(NMAX)
      COMMON /f1com/ pcom,xicom,ncom
      do 11 j=1,ncom
        xt(j)=pcom(j)+x*xicom(j)
11    continue
      call dfunc(xt,df)
      df1dim=0.
      do 12 j=1,ncom
        df1dim=df1dim+df(j)*xicom(j)
12    continue
      return
      END
