      FUNCTION bessi(n,x)
      INTEGER n,IACC
      REAL bessi,x,BIGNO,BIGNI
      PARAMETER (IACC=40,BIGNO=1.0e10,BIGNI=1.0e-10)
CU    USES bessi0
      INTEGER j,m
      REAL bi,bim,bip,tox,bessi0
      if (n.lt.2) pause 'bad argument n in bessi'
      if (x.eq.0.) then
        bessi=0.
      else
        tox=2.0/abs(x)
        bip=0.0
        bi=1.0
        bessi=0.
        m=2*((n+int(sqrt(float(IACC*n)))))
        do 11 j=m,1,-1
          bim=bip+float(j)*tox*bi
          bip=bi
          bi=bim
          if (abs(bi).gt.BIGNO) then
            bessi=bessi*BIGNI
            bi=bi*BIGNI
            bip=bip*BIGNI
          endif
          if (j.eq.n) bessi=bip
11      continue
        bessi=bessi*bessi0(x)/bi
        if (x.lt.0..and.mod(n,2).eq.1) bessi=-bessi
      endif
      return
      END
