      SUBROUTINE poldiv(u,n,v,nv,q,r)
      INTEGER n,nv
      REAL q(n),r(n),u(n),v(nv)
      INTEGER j,k
      do 11 j=1,n
        r(j)=u(j)
        q(j)=0.
11    continue
      do 13 k=n-nv,0,-1
        q(k+1)=r(nv+k)/v(nv)
        do 12 j=nv+k-1,k+1,-1
          r(j)=r(j)-q(k+1)*v(j-k)
12      continue
13    continue
      do 14 j=nv,n
        r(j)=0.
14    continue
      return
      END
