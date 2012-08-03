      SUBROUTINE sprsin(a,n,np,thresh,nmax,sa,ija)
      INTEGER n,nmax,np,ija(nmax)
      REAL thresh,a(np,np),sa(nmax)
      INTEGER i,j,k
      do 11 j=1,n
        sa(j)=a(j,j)
11    continue
      ija(1)=n+2
      k=n+1
      do 13 i=1,n
        do 12 j=1,n
          if(abs(a(i,j)).ge.thresh)then
            if(i.ne.j)then
              k=k+1
              if(k.gt.nmax)pause 'nmax too small in sprsin'
              sa(k)=a(i,j)
              ija(k)=j
            endif
          endif
12      continue
        ija(i+1)=k+1
13    continue
      return
      END
