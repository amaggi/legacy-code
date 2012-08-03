      FUNCTION rtnewt(funcd,x1,x2,xacc)
      INTEGER JMAX
      REAL rtnewt,x1,x2,xacc
      EXTERNAL funcd
      PARAMETER (JMAX=20)
      INTEGER j
      REAL df,dx,f
      rtnewt=.5*(x1+x2)
      do 11 j=1,JMAX
        call funcd(rtnewt,f,df)
        dx=f/df
        rtnewt=rtnewt-dx
        if((x1-rtnewt)*(rtnewt-x2).lt.0.)pause
     *'rtnewt jumped out of brackets'
        if(abs(dx).lt.xacc) return
11    continue
      pause 'rtnewt exceeded maximum iterations'
      END
