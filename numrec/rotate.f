      SUBROUTINE rotate(r,qt,n,np,i,a,b)
      INTEGER n,np,i
      REAL a,b,r(np,np),qt(np,np)
      INTEGER j
      REAL c,fact,s,w,y
      if(a.eq.0.)then
        c=0.
        s=sign(1.,b)
      else if(abs(a).gt.abs(b))then
        fact=b/a
        c=sign(1./sqrt(1.+fact**2),a)
        s=fact*c
      else
        fact=a/b
        s=sign(1./sqrt(1.+fact**2),b)
        c=fact*s
      endif
      do 11 j=i,n
        y=r(i,j)
        w=r(i+1,j)
        r(i,j)=c*y-s*w
        r(i+1,j)=s*y+c*w
11    continue
      do 12 j=1,n
        y=qt(i,j)
        w=qt(i+1,j)
        qt(i,j)=c*y-s*w
        qt(i+1,j)=s*y+c*w
12    continue
      return
      END
