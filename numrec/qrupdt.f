      SUBROUTINE qrupdt(r,qt,n,np,u,v)
      INTEGER n,np
      REAL r(np,np),qt(np,np),u(np),v(np)
CU    USES rotate
      INTEGER i,j,k
      do 11 k=n,1,-1
        if(u(k).ne.0.)goto 1
11    continue
      k=1
1     do 12 i=k-1,1,-1
        call rotate(r,qt,n,np,i,u(i),-u(i+1))
        if(u(i).eq.0.)then
          u(i)=abs(u(i+1))
        else if(abs(u(i)).gt.abs(u(i+1)))then
          u(i)=abs(u(i))*sqrt(1.+(u(i+1)/u(i))**2)
        else
          u(i)=abs(u(i+1))*sqrt(1.+(u(i)/u(i+1))**2)
        endif
12    continue
      do 13 j=1,n
        r(1,j)=r(1,j)+u(1)*v(j)
13    continue
      do 14 i=1,k-1
        call rotate(r,qt,n,np,i,r(i,i),-r(i+1,i))
14    continue
      return
      END
