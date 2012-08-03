      SUBROUTINE eclass(nf,n,lista,listb,m)
      INTEGER m,n,lista(m),listb(m),nf(n)
      INTEGER j,k,l
      do 11 k=1,n
        nf(k)=k
11    continue
      do 12 l=1,m
        j=lista(l)
1       if(nf(j).ne.j)then
          j=nf(j)
        goto 1
        endif
        k=listb(l)
2       if(nf(k).ne.k)then
          k=nf(k)
        goto 2
        endif
        if(j.ne.k)nf(j)=k
12    continue
      do 13 j=1,n
3       if(nf(j).ne.nf(nf(j)))then
          nf(j)=nf(nf(j))
        goto 3
        endif
13    continue
      return
      END
