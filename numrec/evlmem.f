      FUNCTION evlmem(fdt,d,m,xms)
      INTEGER m
      REAL evlmem,fdt,xms,d(m)
      INTEGER i
      REAL sumi,sumr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      theta=6.28318530717959d0*fdt
      wpr=cos(theta)
      wpi=sin(theta)
      wr=1.d0
      wi=0.d0
      sumr=1.
      sumi=0.
      do 11 i=1,m
        wtemp=wr
        wr=wr*wpr-wi*wpi
        wi=wi*wpr+wtemp*wpi
        sumr=sumr-d(i)*sngl(wr)
        sumi=sumi-d(i)*sngl(wi)
11    continue
      evlmem=xms/(sumr**2+sumi**2)
      return
      END
