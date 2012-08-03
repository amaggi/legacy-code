      SUBROUTINE hufdec(ich,code,lcode,nb)
      INTEGER ich,lcode,nb,MC,MQ
      PARAMETER (MC=512,MQ=2*MC-1)
      INTEGER l,nc,nch,node,nodemx
      INTEGER icod(MQ),left(MQ),iright(MQ),ncod(MQ),nprob(MQ)
      LOGICAL btest
      CHARACTER*1 code(lcode)
      COMMON /hufcom/ icod,ncod,nprob,left,iright,nch,nodemx
      SAVE /hufcom/
      node=nodemx
1     continue
        nc=nb/8+1
        if (nc.gt.lcode)then
          ich=nch
          return
        endif
        l=mod(nb,8)
        nb=nb+1
        if(btest(ichar(code(nc)),l))then
          node=iright(node)
        else
          node=left(node)
        endif
        if(node.le.nch)then
          ich=node-1
          return
        endif
      goto 1
      END
