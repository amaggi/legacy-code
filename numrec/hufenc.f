      SUBROUTINE hufenc(ich,code,lcode,nb)
      INTEGER ich,lcode,nb,MC,MQ
      PARAMETER (MC=512,MQ=2*MC-1)
      INTEGER k,l,n,nc,nch,nodemx,ntmp,ibset
      INTEGER icod(MQ),left(MQ),iright(MQ),ncod(MQ),nprob(MQ)
      LOGICAL btest
      CHARACTER*1 code(*)
      COMMON /hufcom/ icod,ncod,nprob,left,iright,nch,nodemx
      SAVE /hufcom/
      k=ich+1
      if(k.gt.nch.or.k.lt.1)pause 'ich out of range in hufenc.'
      do 11 n=ncod(k),1,-1
        nc=nb/8+1
        if (nc.gt.lcode) pause 'lcode too small in hufenc.'
        l=mod(nb,8)
        if (l.eq.0) code(nc)=char(0)
        if(btest(icod(k),n-1))then
          ntmp=ibset(ichar(code(nc)),l)
          code(nc)=char(ntmp)
        endif
        nb=nb+1
11    continue
      return
      END
