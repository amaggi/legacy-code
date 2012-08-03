      SUBROUTINE eulsum(sum,term,jterm,wksp)
      INTEGER jterm
      REAL sum,term,wksp(jterm)
      INTEGER j,nterm
      REAL dum,tmp
      SAVE nterm
      if(jterm.eq.1)then
        nterm=1
        wksp(1)=term
        sum=0.5*term
      else
        tmp=wksp(1)
        wksp(1)=term
        do 11 j=1,nterm-1
          dum=wksp(j+1)
          wksp(j+1)=0.5*(wksp(j)+tmp)
          tmp=dum
11      continue
        wksp(nterm+1)=0.5*(wksp(nterm)+tmp)
        if(abs(wksp(nterm+1)).le.abs(wksp(nterm)))then
          sum=sum+0.5*wksp(nterm+1)
          nterm=nterm+1
        else
          sum=sum+wksp(nterm+1)
        endif
      endif
      return
      END
