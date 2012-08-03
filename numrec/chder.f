      SUBROUTINE chder(a,b,c,cder,n)
      INTEGER n
      REAL a,b,c(n),cder(n)
      INTEGER j
      REAL con
      cder(n)=0.
      cder(n-1)=2*(n-1)*c(n)
      if(n.ge.3)then
        do 11 j=n-2,1,-1
          cder(j)=cder(j+2)+2*j*c(j+1)
11      continue
      endif
      con=2./(b-a)
      do 12 j=1,n
        cder(j)=cder(j)*con
12    continue
      return
      END
