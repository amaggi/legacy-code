      SUBROUTINE hufapp(index,nprob,m,l)
      INTEGER m,l,MC,MQ
      PARAMETER (MC=512,MQ=2*MC-1)
      INTEGER index(MQ),nprob(MQ)
      INTEGER i,j,k,n
      n=m
      i=l
      k=index(i)
2     if(i.le.n/2)then
        j=i+i
        if (j.lt.n.and.nprob(index(j)).gt.nprob(index(j+1))) j=j+1
        if (nprob(k).le.nprob(index(j))) goto 3
        index(i)=index(j)
        i=j
      goto 2
      endif
3     index(i)=k
      return
      END
