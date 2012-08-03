      SUBROUTINE pinvs(ie1,ie2,je1,jsf,jc1,k,c,nci,ncj,nck,s,nsi,nsj)
      INTEGER ie1,ie2,jc1,je1,jsf,k,nci,ncj,nck,nsi,nsj,NMAX
      REAL c(nci,ncj,nck),s(nsi,nsj)
      PARAMETER (NMAX=10)
      INTEGER i,icoff,id,ipiv,irow,j,jcoff,je2,jp,jpiv,js1,indxr(NMAX)
      REAL big,dum,piv,pivinv,pscl(NMAX)
      je2=je1+ie2-ie1
      js1=je2+1
      do 12 i=ie1,ie2
        big=0.
        do 11 j=je1,je2
          if(abs(s(i,j)).gt.big) big=abs(s(i,j))
11      continue
        if(big.eq.0.) pause 'singular matrix, row all 0 in pinvs'
        pscl(i)=1./big
        indxr(i)=0
12    continue
      do 18 id=ie1,ie2
        piv=0.
        do 14 i=ie1,ie2
          if(indxr(i).eq.0) then
            big=0.
            do 13 j=je1,je2
              if(abs(s(i,j)).gt.big) then
                jp=j
                big=abs(s(i,j))
              endif
13          continue
            if(big*pscl(i).gt.piv) then
              ipiv=i
              jpiv=jp
              piv=big*pscl(i)
            endif
          endif
14      continue
        if(s(ipiv,jpiv).eq.0.) pause 'singular matrix in pinvs'
        indxr(ipiv)=jpiv
        pivinv=1./s(ipiv,jpiv)
        do 15 j=je1,jsf
          s(ipiv,j)=s(ipiv,j)*pivinv
15      continue
        s(ipiv,jpiv)=1.
        do 17 i=ie1,ie2
          if(indxr(i).ne.jpiv) then
            if(s(i,jpiv).ne.0.) then
              dum=s(i,jpiv)
              do 16 j=je1,jsf
                s(i,j)=s(i,j)-dum*s(ipiv,j)
16            continue
              s(i,jpiv)=0.
            endif
          endif
17      continue
18    continue
      jcoff=jc1-js1
      icoff=ie1-je1
      do 21 i=ie1,ie2
        irow=indxr(i)+icoff
        do 19 j=js1,jsf
          c(irow,j+jcoff,k)=s(i,j)
19      continue
21    continue
      return
      END
