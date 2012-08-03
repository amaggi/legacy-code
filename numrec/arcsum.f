      SUBROUTINE arcsum(iin,iout,ja,nwk,nrad,nc)
      INTEGER ja,nc,nrad,nwk,iin(*),iout(*)
      INTEGER j,jtmp,karry
      karry=0
      do 11 j=nwk,nc+1,-1
        jtmp=ja
        ja=ja/nrad
        iout(j)=iin(j)+(jtmp-ja*nrad)+karry
        if (iout(j).ge.nrad) then
          iout(j)=iout(j)-nrad
          karry=1
        else
          karry=0
        endif
11    continue
      iout(nc)=iin(nc)+ja+karry
      return
      END
