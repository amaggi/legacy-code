      SUBROUTINE quadvl(x,y,fa,fb,fc,fd)
      REAL fa,fb,fc,fd,x,y
      REAL qa,qb,qc,qd
      qa=min(2.,max(0.,1.-x))
      qb=min(2.,max(0.,1.-y))
      qc=min(2.,max(0.,x+1.))
      qd=min(2.,max(0.,y+1.))
      fa=0.25*qa*qb
      fb=0.25*qb*qc
      fc=0.25*qc*qd
      fd=0.25*qd*qa
      return
      END
