      SUBROUTINE arcode(ich,code,lcode,lcd,isign)
      INTEGER ich,isign,lcd,lcode,MC,NWK
      CHARACTER*1 code(lcode)
      PARAMETER (MC=512,NWK=20)
CU    USES arcsum
      INTEGER ihi,j,ja,jdif,jh,jl,k,m,minint,nc,nch,nrad,ilob(NWK),
     *iupb(NWK),ncumfq(MC+2),ncum,JTRY
      COMMON /arccom/ ncumfq,iupb,ilob,nch,nrad,minint,jdif,nc,ncum
      SAVE /arccom/
      JTRY(j,k,m)=int((dble(k)*dble(j))/dble(m))
      if (isign.eq.0) then
        jdif=nrad-1
        do 11 j=NWK,1,-1
          iupb(j)=nrad-1
          ilob(j)=0
          nc=j
          if(jdif.gt.minint)return
          jdif=(jdif+1)*nrad-1
11      continue
        pause 'NWK too small in arcode'
      else
        if (isign.gt.0) then
          if(ich.gt.nch.or.ich.lt.0)pause 'bad ich in arcode'
        else
          ja=ichar(code(lcd))-ilob(nc)
          do 12 j=nc+1,NWK
            ja=ja*nrad+(ichar(code(j+lcd-nc))-ilob(j))
12        continue
          ich=0
          ihi=nch+1
1         if(ihi-ich.gt.1) then
            m=(ich+ihi)/2
            if (ja.ge.JTRY(jdif,ncumfq(m+1),ncum)) then
              ich=m
            else
              ihi=m
            endif
          goto 1
          endif
          if(ich.eq.nch)return
        endif
        jh=JTRY(jdif,ncumfq(ich+2),ncum)
        jl=JTRY(jdif,ncumfq(ich+1),ncum)
        jdif=jh-jl
        call arcsum(ilob,iupb,jh,NWK,nrad,nc)
        call arcsum(ilob,ilob,jl,NWK,nrad,nc)
        do 13 j=nc,NWK
          if(ich.ne.nch.and.iupb(j).ne.ilob(j))goto 2
          if(lcd.gt.lcode)pause 'lcode too small in arcode'
          if(isign.gt.0) code(lcd)=char(ilob(j))
          lcd=lcd+1
13      continue
        return
2       nc=j
        j=0
3       if (jdif.lt.minint) then
          j=j+1
          jdif=jdif*nrad
        goto 3
        endif
        if (nc-j.lt.1) pause 'NWK too small in arcode'
        if(j.ne.0)then
          do 14 k=nc,NWK
            iupb(k-j)=iupb(k)
            ilob(k-j)=ilob(k)
14        continue
        endif
        nc=nc-j
        do 15 k=NWK-j+1,NWK
          iupb(k)=0
          ilob(k)=0
15      continue
      endif
      return
      END
