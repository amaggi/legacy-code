      SUBROUTINE lfit(x,y,sig,ndat,a,ia,ma,covar,npc,chisq,funcs)
      INTEGER ma,ia(ma),npc,ndat,MMAX
      REAL chisq,a(ma),covar(npc,npc),sig(ndat),x(ndat),y(ndat)
      EXTERNAL funcs
      PARAMETER (MMAX=50)
CU    USES covsrt,gaussj
      INTEGER i,j,k,l,m,mfit
      REAL sig2i,sum,wt,ym,afunc(MMAX),beta(MMAX)
      mfit=0
      do 11 j=1,ma
        if(ia(j).ne.0) mfit=mfit+1
11    continue
      if(mfit.eq.0) pause 'lfit: no parameters to be fitted'
      do 13 j=1,mfit
        do 12 k=1,mfit
          covar(j,k)=0.
12      continue
        beta(j)=0.
13    continue
      do 17 i=1,ndat
        call funcs(x(i),afunc,ma)
        ym=y(i)
        if(mfit.lt.ma) then
          do 14 j=1,ma
            if(ia(j).eq.0) ym=ym-a(j)*afunc(j)
14        continue
        endif
        sig2i=1./sig(i)**2
        j=0
        do 16 l=1,ma
          if (ia(l).ne.0) then
            j=j+1
            wt=afunc(l)*sig2i
            k=0
            do 15 m=1,l
              if (ia(m).ne.0) then
                k=k+1
                covar(j,k)=covar(j,k)+wt*afunc(m)
              endif
15          continue
            beta(j)=beta(j)+ym*wt
          endif
16      continue
17    continue
      do 19 j=2,mfit
        do 18 k=1,j-1
          covar(k,j)=covar(j,k)
18      continue
19    continue
      call gaussj(covar,mfit,npc,beta,1,1)
      j=0
      do 21 l=1,ma
        if(ia(l).ne.0) then
          j=j+1
          a(l)=beta(j)
        endif
21    continue
      chisq=0.
      do 23 i=1,ndat
        call funcs(x(i),afunc,ma)
        sum=0.
        do 22 j=1,ma
          sum=sum+a(j)*afunc(j)
22      continue
        chisq=chisq+((y(i)-sum)/sig(i))**2
23    continue
      call covsrt(covar,npc,ma,ia,mfit)
      return
      END
