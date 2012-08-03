      SUBROUTINE arcmak(nfreq,nchh,nradd)
      INTEGER nchh,nradd,nfreq(nchh),MC,NWK,MAXINT
      PARAMETER (MC=512,NWK=20,MAXINT=2147483647)
      INTEGER j,jdif,minint,nc,nch,nrad,ncum,ncumfq(MC+2),ilob(NWK),
     *iupb(NWK)
      COMMON /arccom/ ncumfq,iupb,ilob,nch,nrad,minint,jdif,nc,ncum
      SAVE /arccom/
      if(nchh.gt.MC)pause 'MC too small in arcmak'
      if(nradd.gt.256)pause 'nradd may not exceed 256 in arcmak'
      minint=MAXINT/nradd
      nch=nchh
      nrad=nradd
      ncumfq(1)=0
      do 11 j=2,nch+1
        ncumfq(j)=ncumfq(j-1)+max(nfreq(j-1),1)
11    continue
      ncumfq(nch+2)=ncumfq(nch+1)+1
      ncum=ncumfq(nch+2)
      return
      END
