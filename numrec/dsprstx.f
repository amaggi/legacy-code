      SUBROUTINE dsprstx(sa,ija,x,b,n)
      INTEGER n,ija(*)
      DOUBLE PRECISION b(n),sa(*),x(n)
      INTEGER i,j,k
      if (ija(1).ne.n+2) pause 'mismatched vector and matrix in sprstx'
      do 11 i=1,n
        b(i)=sa(i)*x(i)
11    continue
      do 13 i=1,n
        do 12 k=ija(i),ija(i+1)-1
          j=ija(k)
          b(j)=b(j)+sa(k)*x(i)
12      continue
13    continue
      return
      END
