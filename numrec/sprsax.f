      SUBROUTINE sprsax(sa,ija,x,b,n)
      INTEGER n,ija(*)
      REAL b(n),sa(*),x(n)
      INTEGER i,k
      if (ija(1).ne.n+2) pause 'mismatched vector and matrix in sprsax'
      do 12 i=1,n
        b(i)=sa(i)*x(i)
        do 11 k=ija(i),ija(i+1)-1
          b(i)=b(i)+sa(k)*x(ija(k))
11      continue
12    continue
      return
      END
