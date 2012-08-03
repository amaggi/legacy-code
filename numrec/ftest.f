      SUBROUTINE ftest(data1,n1,data2,n2,f,prob)
      INTEGER n1,n2
      REAL f,prob,data1(n1),data2(n2)
CU    USES avevar,betai
      REAL ave1,ave2,df1,df2,var1,var2,betai
      call avevar(data1,n1,ave1,var1)
      call avevar(data2,n2,ave2,var2)
      if(var1.gt.var2)then
        f=var1/var2
        df1=n1-1
        df2=n2-1
      else
        f=var2/var1
        df1=n2-1
        df2=n1-1
      endif
      prob=2.*betai(0.5*df2,0.5*df1,df2/(df2+df1*f))
      if(prob.gt.1.)prob=2.-prob
      return
      END
