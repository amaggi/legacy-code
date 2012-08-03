      SUBROUTINE predic(data,ndata,d,m,future,nfut)
      INTEGER ndata,nfut,m,MMAX
      REAL d(m),data(ndata),future(nfut)
      PARAMETER (MMAX=100)
      INTEGER j,k
      REAL discrp,sum,reg(MMAX)
      do 11 j=1,m
        reg(j)=data(ndata+1-j)
11    continue
      do 14 j=1,nfut
        discrp=0.
        sum=discrp
        do 12 k=1,m
          sum=sum+d(k)*reg(k)
12      continue
        do 13 k=m,2,-1
          reg(k)=reg(k-1)
13      continue
        reg(1)=sum
        future(j)=sum
14    continue
      return
      END
