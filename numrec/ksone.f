      SUBROUTINE ksone(data,n,func,d,prob)
      INTEGER n
      REAL d,data(n),func,prob
      EXTERNAL func
CU    USES probks,sort
      INTEGER j
      REAL dt,en,ff,fn,fo,probks
      call sort(n,data)
      en=n
      d=0.
      fo=0.
      do 11 j=1,n
        fn=j/en
        ff=func(data(j))
        dt=max(abs(fo-ff),abs(fn-ff))
        if(dt.gt.d)d=dt
        fo=fn
11    continue
      en=sqrt(en)
      prob=probks((en+0.12+0.11/en)*d)
      return
      END
