      FUNCTION selip(k,n,arr)
      INTEGER k,n,M
      REAL selip,arr(n),BIG
      PARAMETER (M=64,BIG=1.E30)
CU    USES shell
      INTEGER i,j,jl,jm,ju,kk,mm,nlo,nxtmm,isel(M+2)
      REAL ahi,alo,sum,sel(M+2)
      if(k.lt.1.or.k.gt.n.or.n.le.0) pause 'bad input to selip'
      kk=k
      ahi=BIG
      alo=-BIG
1     continue
        mm=0
        nlo=0
        sum=0.
        nxtmm=M+1
        do 11 i=1,n
          if(arr(i).ge.alo.and.arr(i).le.ahi)then
            mm=mm+1
            if(arr(i).eq.alo) nlo=nlo+1
            if(mm.le.M)then
              sel(mm)=arr(i)
            else if(mm.eq.nxtmm)then
              nxtmm=mm+mm/M
              sel(1+mod(i+mm+kk,M))=arr(i)
            endif
            sum=sum+arr(i)
          endif
11      continue
        if(kk.le.nlo)then
          selip=alo
          return
        else if(mm.le.M)then
          call shell(mm,sel)
          selip=sel(kk)
          return
        endif
        sel(M+1)=sum/mm
        call shell(M+1,sel)
        sel(M+2)=ahi
        do 12 j=1,M+2
          isel(j)=0
12      continue
        do 13 i=1,n
          if(arr(i).ge.alo.and.arr(i).le.ahi)then
            jl=0
            ju=M+2
2           if(ju-jl.gt.1)then
              jm=(ju+jl)/2
              if(arr(i).ge.sel(jm))then
                jl=jm
              else
                ju=jm
              endif
            goto 2
            endif
            isel(ju)=isel(ju)+1
          endif
13      continue
        j=1
3       if(kk.gt.isel(j))then
          alo=sel(j)
          kk=kk-isel(j)
          j=j+1
        goto 3
        endif
        ahi=sel(j)
      goto 1
      END
