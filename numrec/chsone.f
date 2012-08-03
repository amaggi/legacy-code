      SUBROUTINE chsone(bins,ebins,nbins,knstrn,df,chsq,prob)
      INTEGER knstrn,nbins
      REAL chsq,df,prob,bins(nbins),ebins(nbins)
CU    USES gammq
      INTEGER j
      REAL gammq
      df=nbins-knstrn
      chsq=0.
      do 11 j=1,nbins
        if(ebins(j).le.0.)pause 'bad expected number in chsone'
        chsq=chsq+(bins(j)-ebins(j))**2/ebins(j)
11    continue
      prob=gammq(0.5*df,0.5*chsq)
      return
      END
