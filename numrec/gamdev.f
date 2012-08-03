      FUNCTION gamdev(ia,idum)
      INTEGER ia,idum
      REAL gamdev
CU    USES ran1
      INTEGER j
      REAL am,e,s,v1,v2,x,y,ran1
      if(ia.lt.1)pause 'bad argument in gamdev'
      if(ia.lt.6)then
        x=1.
        do 11 j=1,ia
          x=x*ran1(idum)
11      continue
        x=-log(x)
      else
1         v1=2.*ran1(idum)-1.
          v2=2.*ran1(idum)-1.
        if(v1**2+v2**2.gt.1.)goto 1
          y=v2/v1
          am=ia-1
          s=sqrt(2.*am+1.)
          x=s*y+am
        if(x.le.0.)goto 1
          e=(1.+y**2)*exp(am*log(x/am)-s*y)
        if(ran1(idum).gt.e)goto 1
      endif
      gamdev=x
      return
      END
