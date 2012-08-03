      SUBROUTINE hufmak(nfreq,nchin,ilong,nlong)
      INTEGER ilong,nchin,nlong,nfreq(nchin),MC,MQ
      PARAMETER (MC=512,MQ=2*MC-1)
CU    USES hufapp
      INTEGER ibit,j,k,n,nch,node,nodemx,nused,ibset,index(MQ),iup(MQ),
     *icod(MQ),left(MQ),iright(MQ),ncod(MQ),nprob(MQ)
      COMMON /hufcom/ icod,ncod,nprob,left,iright,nch,nodemx
      SAVE /hufcom/
      nch=nchin
      nused=0
      do 11 j=1,nch
        nprob(j)=nfreq(j)
        icod(j)=0
        ncod(j)=0
        if(nfreq(j).ne.0)then
          nused=nused+1
          index(nused)=j
        endif
11    continue
      do 12 j=nused,1,-1
        call hufapp(index,nprob,nused,j)
12    continue
      k=nch
1     if(nused.gt.1)then
      node=index(1)
        index(1)=index(nused)
        nused=nused-1
        call hufapp(index,nprob,nused,1)
        k=k+1
        nprob(k)=nprob(index(1))+nprob(node)
        left(k)=node
        iright(k)=index(1)
        iup(index(1)) = -k
        iup(node)=k
        index(1)=k
        call hufapp(index,nprob,nused,1)
      goto 1
      endif
      nodemx=k
      iup(nodemx)=0
      do 13 j=1,nch
        if(nprob(j).ne.0)then
          n=0
          ibit=0
          node=iup(j)
2         if(node.ne.0)then
            if(node.lt.0)then
              n=ibset(n,ibit)
              node = -node
            endif
            node=iup(node)
            ibit=ibit+1
          goto 2
          endif
          icod(j)=n
          ncod(j)=ibit
        endif
13    continue
      nlong=0
      do 14 j=1,nch
        if(ncod(j).gt.nlong)then
          nlong=ncod(j)
          ilong=j-1
        endif
14    continue
      return
      END
