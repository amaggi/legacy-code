      SUBROUTINE balanc(a,n,np)
      INTEGER n,np
      REAL a(np,np),RADIX,SQRDX
      PARAMETER (RADIX=2.,SQRDX=RADIX**2)
      INTEGER i,j,last
      REAL c,f,g,r,s
1     continue
        last=1
        do 14 i=1,n
          c=0.
          r=0.
          do 11 j=1,n
            if(j.ne.i)then
              c=c+abs(a(j,i))
              r=r+abs(a(i,j))
            endif
11        continue
          if(c.ne.0..and.r.ne.0.)then
            g=r/RADIX
            f=1.
            s=c+r
2           if(c.lt.g)then
              f=f*RADIX
              c=c*SQRDX
            goto 2
            endif
            g=r*RADIX
3           if(c.gt.g)then
              f=f/RADIX
              c=c/SQRDX
            goto 3
            endif
            if((c+r)/f.lt.0.95*s)then
              last=0
              g=1./f
              do 12 j=1,n
                a(i,j)=a(i,j)*g
12            continue
              do 13 j=1,n
                a(j,i)=a(j,i)*f
13            continue
            endif
          endif
14      continue
      if(last.eq.0)goto 1
      return
      END
