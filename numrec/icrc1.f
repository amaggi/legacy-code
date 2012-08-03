      FUNCTION icrc1(crc,onech,ib1,ib2,ib3)
      INTEGER icrc1,ib1,ib2,ib3
      INTEGER i,ichr,ireg
      CHARACTER*1 onech,crc(4),creg(4)
      EQUIVALENCE (creg,ireg)
      ireg=0
      creg(ib1)=crc(ib1)
      creg(ib2)=char(ieor(ichar(crc(ib2)),ichar(onech)))
      do 11 i=1,8
        ichr=ichar(creg(ib2))
        ireg=ireg+ireg
        creg(ib3)=char(0)
        if(ichr.gt.127)ireg=ieor(ireg,4129)
11    continue
      icrc1=ireg
      return
      END
