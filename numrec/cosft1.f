      SUBROUTINE cosft1(y,n)
      INTEGER n
      REAL y(n+1)
CU    USES realft
      INTEGER j
      REAL sum,y1,y2
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      theta=3.141592653589793d0/n
      wr=1.0d0
      wi=0.0d0
      wpr=-2.0d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      sum=0.5*(y(1)-y(n+1))
      y(1)=0.5*(y(1)+y(n+1))
      do 11 j=1,n/2-1
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
        y1=0.5*(y(j+1)+y(n-j+1))
        y2=(y(j+1)-y(n-j+1))
        y(j+1)=y1-wi*y2
        y(n-j+1)=y1+wi*y2
        sum=sum+wr*y2
11    continue
      call realft(y,n,+1)
      y(n+1)=y(2)
      y(2)=sum
      do 12 j=4,n,2
        sum=sum+y(j)
        y(j)=sum
12    continue
      return
      END
