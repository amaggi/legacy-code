      SUBROUTINE psdes(lword,irword)
      INTEGER irword,lword,NITER
      PARAMETER (NITER=4)
      INTEGER i,ia,ib,iswap,itmph,itmpl,c1(4),c2(4)
      SAVE c1,c2
      DATA c1 /Z'BAA96887',Z'1E17D32C',Z'03BCDC3C',Z'0F33D1B2'/, c2 
     */Z'4B0F3B58',Z'E874F0C3',Z'6955C5A6', Z'55A7CA46'/
      do 11 i=1,NITER
        iswap=irword
        ia=ieor(irword,c1(i))
        itmpl=iand(ia,65535)
        itmph=iand(ishft(ia,-16),65535)
        ib=itmpl**2+not(itmph**2)
        ia=ior(ishft(ib,16),iand(ishft(ib,-16),65535))
        irword=ieor(lword,ieor(c2(i),ia)+itmpl*itmph)
        lword=iswap
11    continue
      return
      END
