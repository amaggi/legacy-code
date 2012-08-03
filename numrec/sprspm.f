      SUBROUTINE sprspm(sa,ija,sb,ijb,sc,ijc)
      INTEGER ija(*),ijb(*),ijc(*)
      REAL sa(*),sb(*),sc(*)
      INTEGER i,ijma,ijmb,j,m,ma,mb,mbb,mn
      REAL sum
      if (ija(1).ne.ijb(1).or.ija(1).ne.ijc(1))pause
     *'sprspm sizes do not match'
      do 13 i=1,ijc(1)-2
        j=i
        m=i
        mn=ijc(i)
        sum=sa(i)*sb(i)
1       continue
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
          sc(m)=sum
          sum=0.e0
          if(mn.ge.ijc(i+1))goto 3
          m=mn
          mn=mn+1
          j=ijc(m)
        goto 1
3       continue
13    continue
      return
      END
