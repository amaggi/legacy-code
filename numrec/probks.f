      FUNCTION probks(alam)
      REAL probks,alam,EPS1,EPS2
      PARAMETER (EPS1=0.001, EPS2=1.e-8)
      INTEGER j
      REAL a2,fac,term,termbf
      a2=-2.*alam**2
      fac=2.
      probks=0.
      termbf=0.
      do 11 j=1,100
        term=fac*exp(a2*j**2)
        probks=probks+term
        if(abs(term).le.EPS1*termbf.or.abs(term).le.EPS2*probks)return
        fac=-fac
        termbf=abs(term)
11    continue
      probks=1.
      return
      END
