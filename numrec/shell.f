      SUBROUTINE shell(n,a)
      INTEGER n
      REAL a(n)
      INTEGER i,j,inc
      REAL v
      inc=1
1     inc=3*inc+1
      if(inc.le.n)goto 1
2     continue
        inc=inc/3
        do 11 i=inc+1,n
          v=a(i)
          j=i
3         if(a(j-inc).gt.v)then
            a(j)=a(j-inc)
            j=j-inc
            if(j.le.inc)goto 4
          goto 3
          endif
4         a(j)=v
11      continue
      if(inc.gt.1)goto 2
      return
      END
