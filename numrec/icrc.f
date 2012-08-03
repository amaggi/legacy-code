      FUNCTION icrc(crc,bufptr,len,jinit,jrev)
      INTEGER icrc,jinit,jrev,len
      CHARACTER*1 bufptr(*),crc(2)
CU    USES icrc1
      INTEGER ich,init,ireg,j,icrctb(0:255),it(0:15),icrc1,ib1,ib2,ib3
      CHARACTER*1 creg(4),rchr(0:255)
      SAVE icrctb,rchr,init,it,ib1,ib2,ib3
      EQUIVALENCE (creg,ireg)
      DATA it/0,8,4,12,2,10,6,14,1,9,5,13,3,11,7,15/, init /0/
      if (init.eq.0) then
        init=1
        ireg=256*(256*ichar('3')+ichar('2'))+ichar('1')
        do 11 j=1,4
          if (creg(j).eq.'1') ib1=j
          if (creg(j).eq.'2') ib2=j
          if (creg(j).eq.'3') ib3=j
11      continue
        do 12 j=0,255
          ireg=j*256
          icrctb(j)=icrc1(creg,char(0),ib1,ib2,ib3)
          ich=it(mod(j,16))*16+it(j/16)
          rchr(j)=char(ich)
12      continue
      endif
      if (jinit.ge.0) then
        crc(1)=char(jinit)
        crc(2)=char(jinit)
      else if (jrev.lt.0) then
        ich=ichar(crc(1))
        crc(1)=rchr(ichar(crc(2)))
        crc(2)=rchr(ich)
      endif
      do 13 j=1,len
        ich=ichar(bufptr(j))
        if(jrev.lt.0)ich=ichar(rchr(ich))
        ireg=icrctb(ieor(ich,ichar(crc(2))))
        crc(2)=char(ieor(ichar(creg(ib2)),ichar(crc(1))))
        crc(1)=creg(ib1)
13    continue
      if (jrev.ge.0) then
        creg(ib1)=crc(1)
        creg(ib2)=crc(2)
      else
        creg(ib2)=rchr(ichar(crc(1)))
        creg(ib1)=rchr(ichar(crc(2)))
        crc(1)=creg(ib1)
        crc(2)=creg(ib2)
      endif
      icrc=ireg
      return
      END
