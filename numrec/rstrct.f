      SUBROUTINE rstrct(uc,uf,nc)
      INTEGER nc
      DOUBLE PRECISION uc(nc,nc),uf(2*nc-1,2*nc-1)
      INTEGER ic,if,jc,jf
      do 12 jc=2,nc-1
        jf=2*jc-1
        do 11 ic=2,nc-1
          if=2*ic-1
          uc(ic,jc)=.5d0*uf(if,jf)+.125d0*(uf(if+1,jf)+uf(if-1,jf)+
     *uf(if,jf+1)+uf(if,jf-1))
11      continue
12    continue
      do 13 ic=1,nc
        uc(ic,1)=uf(2*ic-1,1)
        uc(ic,nc)=uf(2*ic-1,2*nc-1)
13    continue
      do 14 jc=1,nc
        uc(1,jc)=uf(1,2*jc-1)
        uc(nc,jc)=uf(2*nc-1,2*jc-1)
14    continue
      return
      END
