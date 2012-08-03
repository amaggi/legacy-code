      SUBROUTINE addint(uf,uc,res,nf)
      INTEGER nf
      DOUBLE PRECISION res(nf,nf),uc(nf/2+1,nf/2+1),uf(nf,nf)
CU    USES interp
      INTEGER i,j
      call interp(res,uc,nf)
      do 12 j=1,nf
        do 11 i=1,nf
          uf(i,j)=uf(i,j)+res(i,j)
11      continue
12    continue
      return
      END
