      FUNCTION ellpi(phi,en,ak)
      REAL ellpi,ak,en,phi
CU    USES rf,rj
      REAL cc,enss,q,s,rf,rj
      s=sin(phi)
      enss=en*s*s
      cc=cos(phi)**2
      q=(1.-s*ak)*(1.+s*ak)
      ellpi=s*(rf(cc,q,1.)-enss*rj(cc,q,1.,1.+enss)/3.)
      return
      END
