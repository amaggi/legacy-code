      SUBROUTINE orthog(n,anu,alpha,beta,a,b)
      INTEGER n,NMAX
      REAL a(n),alpha(2*n-1),anu(2*n),b(n),beta(2*n-1)
      PARAMETER (NMAX=64)
      INTEGER k,l
      REAL sig(2*NMAX+1,2*NMAX+1)
      do 11 l=3,2*n
        sig(1,l)=0.
11    continue
      do 12 l=2,2*n+1
        sig(2,l)=anu(l-1)
12    continue
      a(1)=alpha(1)+anu(2)/anu(1)
      b(1)=0.
      do 14 k=3,n+1
        do 13 l=k,2*n-k+3
          sig(k,l)=sig(k-1,l+1)+(alpha(l-1)-a(k-2))*sig(k-1,l)-b(k-2)*
     *sig(k-2,l)+beta(l-1)*sig(k-1,l-1)
13      continue
        a(k-1)=alpha(k-1)+sig(k,k+1)/sig(k,k)-sig(k-1,k)/sig(k-1,k-1)
        b(k-1)=sig(k,k)/sig(k-1,k-1)
14    continue
      return
      END
