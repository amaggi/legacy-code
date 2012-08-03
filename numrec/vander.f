      SUBROUTINE vander(x,w,q,n)
      INTEGER n,NMAX
      DOUBLE PRECISION q(n),w(n),x(n)
      PARAMETER (NMAX=100)
      INTEGER i,j,k,k1
      DOUBLE PRECISION b,s,t,xx,c(NMAX)
      if(n.eq.1)then
        w(1)=q(1)
      else
        do 11 i=1,n
          c(i)=0.d0
11      continue
        c(n)=-x(1)
        do 13 i=2,n
          xx=-x(i)
          do 12 j=n+1-i,n-1
            c(j)=c(j)+xx*c(j+1)
12        continue
          c(n)=c(n)+xx
13      continue
        do 15 i=1,n
          xx=x(i)
          t=1.d0
          b=1.d0
          s=q(n)
          k=n
          do 14 j=2,n
            k1=k-1
            b=c(k)+xx*b
            s=s+q(k1)*b
            t=xx*t+b
            k=k1
14        continue
          w(i)=s/t
15      continue
      endif
      return
      END
