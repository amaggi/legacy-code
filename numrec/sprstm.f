      SUBROUTINE sprstm(sa,ija,sb,ijb,thresh,nmax,sc,ijc)
      INTEGER nmax,ija(*),ijb(*),ijc(nmax)
      REAL thresh,sa(*),sb(*),sc(nmax)
      INTEGER i,ijma,ijmb,j,k,ma,mb,mbb
      REAL sum
      if (ija(1).ne.ijb(1)) pause 'sprstm sizes do not match'
      k=ija(1)
      ijc(1)=k
      do 14 i=1,ija(1)-2
        do 13 j=1,ijb(1)-2
          if(i.eq.j)then
            sum=sa(i)*sb(j)
          else
            sum=0.d0
          endif
          mb=ijb(j)
          do 11 ma=ija(i),ija(i+1)-1
            ijma=ija(ma)
            if(ijma.eq.j)then
              sum=sum+sa(ma)*sb(j)
            else
2             if(mb.lt.ijb(j+1))then
                ijmb=ijb(mb)
                if(ijmb.eq.i)then
                  sum=sum+sa(i)*sb(mb)
                  mb=mb+1
                  goto 2
                else if(ijmb.lt.ijma)then
                  mb=mb+1
                  goto 2
                else if(ijmb.eq.ijma)then
                  sum=sum+sa(ma)*sb(mb)
                  mb=mb+1
                  goto 2
                endif
              endif
            endif
11        continue
          do 12 mbb=mb,ijb(j+1)-1
            if(ijb(mbb).eq.i)then
              sum=sum+sa(i)*sb(mbb)
            endif
12        continue
          if(i.eq.j)then
            sc(i)=sum
          else if(abs(sum).gt.thresh)then
            if(k.gt.nmax)pause 'sprstm: nmax to small'
            sc(k)=sum
            ijc(k)=j
            k=k+1
          endif
13      continue
        ijc(i+1)=k
14    continue
      return
      END
