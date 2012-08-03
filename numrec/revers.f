      SUBROUTINE revers(iorder,ncity,n)
      INTEGER ncity,iorder(ncity),n(6)
      INTEGER itmp,j,k,l,nn
      nn=(1+mod(n(2)-n(1)+ncity,ncity))/2
      do 11 j=1,nn
        k=1+mod((n(1)+j-2),ncity)
        l=1+mod((n(2)-j+ncity),ncity)
        itmp=iorder(k)
        iorder(k)=iorder(l)
        iorder(l)=itmp
11    continue
      return
      END
