      SUBROUTINE interp(uf,uc,nf)
      INTEGER nf
      DOUBLE PRECISION uc(nf/2+1,nf/2+1),uf(nf,nf)
      INTEGER ic,if,jc,jf,nc
      nc=nf/2+1
      do 12 jc=1,nc
        jf=2*jc-1
        do 11 ic=1,nc
          uf(2*ic-1,jf)=uc(ic,jc)
11      continue
12    continue
      do 14 jf=1,nf,2
        do 13 if=2,nf-1,2
          uf(if,jf)=.5d0*(uf(if+1,jf)+uf(if-1,jf))
13      continue
14    continue
      do 16 jf=2,nf-1,2
        do 15 if=1,nf
          uf(if,jf)=.5d0*(uf(if,jf+1)+uf(if,jf-1))
15      continue
16    continue
      return
      END
