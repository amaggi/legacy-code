      SUBROUTINE tutest(data1,n1,data2,n2,t,prob)
      INTEGER n1,n2
      REAL prob,t,data1(n1),data2(n2)
CU    USES avevar,betai
      REAL ave1,ave2,df,var1,var2,betai
      call avevar(data1,n1,ave1,var1)
      call avevar(data2,n2,ave2,var2)
      t=(ave1-ave2)/sqrt(var1/n1+var2/n2)
      df=(var1/n1+var2/n2)**2/((var1/n1)**2/(n1-1)+(var2/n2)**2/(n2-1))
      prob=betai(0.5*df,0.5,df/(df+t**2))
      return
      END
