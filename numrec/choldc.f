      SUBROUTINE choldc(a,n,np,p)
      INTEGER n,np
      REAL a(np,np),p(n)
      INTEGER i,j,k
      REAL sum
      do 13 i=1,n
        do 12 j=i,n
          sum=a(i,j)
          do 11 k=i-1,1,-1
            sum=sum-a(i,k)*a(j,k)
11        continue
          if(i.eq.j)then
            if(sum.le.0.)pause 'choldc failed'
            p(i)=sqrt(sum)
          else
            a(j,i)=sum/p(i)
          endif
12      continue
13    continue
      return
      END
