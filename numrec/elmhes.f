      SUBROUTINE elmhes(a,n,np)
      INTEGER n,np
      REAL a(np,np)
      INTEGER i,j,m
      REAL x,y
      do 17 m=2,n-1
        x=0.
        i=m
        do 11 j=m,n
          if(abs(a(j,m-1)).gt.abs(x))then
            x=a(j,m-1)
            i=j
          endif
11      continue
        if(i.ne.m)then
          do 12 j=m-1,n
            y=a(i,j)
            a(i,j)=a(m,j)
            a(m,j)=y
12        continue
          do 13 j=1,n
            y=a(j,i)
            a(j,i)=a(j,m)
            a(j,m)=y
13        continue
        endif
        if(x.ne.0.)then
          do 16 i=m+1,n
            y=a(i,m-1)
            if(y.ne.0.)then
              y=y/x
              a(i,m-1)=y
              do 14 j=m,n
                a(i,j)=a(i,j)-y*a(m,j)
14            continue
              do 15 j=1,n
                a(j,m)=a(j,m)+y*a(j,i)
15            continue
            endif
16        continue
        endif
17    continue
      return
      END
