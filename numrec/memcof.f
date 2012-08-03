      SUBROUTINE memcof(data,n,m,xms,d)
      INTEGER m,n,MMAX,NMAX
      REAL xms,d(m),data(n)
      PARAMETER (MMAX=60,NMAX=2000)
      INTEGER i,j,k
      REAL denom,p,pneum,wk1(NMAX),wk2(NMAX),wkm(MMAX)
      if (m.gt.MMAX.or.n.gt.NMAX) pause 'workspace too small in memcof'
      p=0.
      do 11 j=1,n
        p=p+data(j)**2
11    continue
      xms=p/n
      wk1(1)=data(1)
      wk2(n-1)=data(n)
      do 12 j=2,n-1
        wk1(j)=data(j)
        wk2(j-1)=data(j)
12    continue
      do 17 k=1,m
        pneum=0.
        denom=0.
        do 13 j=1,n-k
          pneum=pneum+wk1(j)*wk2(j)
          denom=denom+wk1(j)**2+wk2(j)**2
13      continue
        d(k)=2.*pneum/denom
        xms=xms*(1.-d(k)**2)
        do 14 i=1,k-1
          d(i)=wkm(i)-d(k)*wkm(k-i)
14      continue
        if(k.eq.m)return
        do 15 i=1,k
          wkm(i)=d(i)
15      continue
        do 16 j=1,n-k-1
          wk1(j)=wk1(j)-wkm(k)*wk2(j)
          wk2(j)=wk2(j+1)-wkm(k)*wk1(j+1)
16      continue
17    continue
      pause 'never get here in memcof'
      END
