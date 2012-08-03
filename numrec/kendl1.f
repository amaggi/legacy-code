      SUBROUTINE kendl1(data1,data2,n,tau,z,prob)
      INTEGER n
      REAL prob,tau,z,data1(n),data2(n)
CU    USES erfcc
      INTEGER is,j,k,n1,n2
      REAL a1,a2,aa,var,erfcc
      n1=0
      n2=0
      is=0
      do 12 j=1,n-1
        do 11 k=j+1,n
          a1=data1(j)-data1(k)
          a2=data2(j)-data2(k)
          aa=a1*a2
          if(aa.ne.0.)then
            n1=n1+1
            n2=n2+1
            if(aa.gt.0.)then
              is=is+1
            else
              is=is-1
            endif
          else
            if(a1.ne.0.)n1=n1+1
            if(a2.ne.0.)n2=n2+1
          endif
11      continue
12    continue
      tau=float(is)/sqrt(float(n1)*float(n2))
      var=(4.*n+10.)/(9.*n*(n-1.))
      z=tau/sqrt(var)
      prob=erfcc(abs(z)/1.4142136)
      return
      END
