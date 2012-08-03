      FUNCTION bessk(n,x)
      INTEGER n
      REAL bessk,x
CU    USES bessk0,bessk1
      INTEGER j
      REAL bk,bkm,bkp,tox,bessk0,bessk1
      if (n.lt.2) pause 'bad argument n in bessk'
      tox=2.0/x
      bkm=bessk0(x)
      bk=bessk1(x)
      do 11 j=1,n-1
        bkp=bkm+j*tox*bk
        bkm=bk
        bk=bkp
11    continue
      bessk=bk
      return
      END
