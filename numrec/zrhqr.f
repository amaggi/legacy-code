      SUBROUTINE zrhqr(a,m,rtr,rti)
      INTEGER m,MAXM
      REAL a(m+1),rtr(m),rti(m)
      PARAMETER (MAXM=50)
CU    USES balanc,hqr
      INTEGER j,k
      REAL hess(MAXM,MAXM),xr,xi
      if (m.gt.MAXM.or.a(m+1).eq.0.) pause 'bad args in zrhqr'
      do 12 k=1,m
        hess(1,k)=-a(m+1-k)/a(m+1)
        do 11 j=2,m
          hess(j,k)=0.
11      continue
        if (k.ne.m) hess(k+1,k)=1.
12    continue
      call balanc(hess,m,MAXM)
      call hqr(hess,m,MAXM,rtr,rti)
      do 14 j=2,m
        xr=rtr(j)
        xi=rti(j)
        do 13 k=j-1,1,-1
          if(rtr(k).le.xr)goto 1
          rtr(k+1)=rtr(k)
          rti(k+1)=rti(k)
13      continue
        k=0
1       rtr(k+1)=xr
        rti(k+1)=xi
14    continue
      return
      END
