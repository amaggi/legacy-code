      SUBROUTINE cntab1(nn,ni,nj,chisq,df,prob,cramrv,ccc)
      INTEGER ni,nj,nn(ni,nj),MAXI,MAXJ
      REAL ccc,chisq,cramrv,df,prob,TINY
      PARAMETER (MAXI=100,MAXJ=100,TINY=1.e-30)
CU    USES gammq
      INTEGER i,j,nni,nnj
      REAL expctd,sum,sumi(MAXI),sumj(MAXJ),gammq
      sum=0
      nni=ni
      nnj=nj
      do 12 i=1,ni
        sumi(i)=0.
        do 11 j=1,nj
          sumi(i)=sumi(i)+nn(i,j)
          sum=sum+nn(i,j)
11      continue
        if(sumi(i).eq.0.)nni=nni-1
12    continue
      do 14 j=1,nj
        sumj(j)=0.
        do 13 i=1,ni
          sumj(j)=sumj(j)+nn(i,j)
13      continue
        if(sumj(j).eq.0.)nnj=nnj-1
14    continue
      df=nni*nnj-nni-nnj+1
      chisq=0.
      do 16 i=1,ni
        do 15 j=1,nj
          expctd=sumj(j)*sumi(i)/sum
          chisq=chisq+(nn(i,j)-expctd)**2/(expctd+TINY)
15      continue
16    continue
      prob=gammq(0.5*df,0.5*chisq)
      cramrv=sqrt(chisq/(sum*min(nni-1,nnj-1)))
      ccc=sqrt(chisq/(chisq+sum))
      return
      END
