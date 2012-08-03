      SUBROUTINE cosft2(y,n,isign)
      INTEGER isign,n
      REAL y(n)
CU    USES realft
      INTEGER i
      REAL sum,sum1,y1,y2,ytemp
      DOUBLE PRECISION theta,wi,wi1,wpi,wpr,wr,wr1,wtemp,PI
      PARAMETER (PI=3.141592653589793d0)
      theta=0.5d0*PI/n
      wr=1.0d0
      wi=0.0d0
      wr1=cos(theta)
      wi1=sin(theta)
      wpr=-2.0d0*wi1**2
      wpi=sin(2.d0*theta)
      if(isign.eq.1)then
        do 11 i=1,n/2
          y1=0.5*(y(i)+y(n-i+1))
          y2=wi1*(y(i)-y(n-i+1))
          y(i)=y1+y2
          y(n-i+1)=y1-y2
          wtemp=wr1
          wr1=wr1*wpr-wi1*wpi+wr1
          wi1=wi1*wpr+wtemp*wpi+wi1
11      continue
        call realft(y,n,1)
        do 12 i=3,n,2
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
          y1=y(i)*wr-y(i+1)*wi
          y2=y(i+1)*wr+y(i)*wi
          y(i)=y1
          y(i+1)=y2
12      continue
        sum=0.5*y(2)
        do 13 i=n,2,-2
          sum1=sum
          sum=sum+y(i)
          y(i)=sum1
13      continue
      else if(isign.eq.-1)then
        ytemp=y(n)
        do 14 i=n,4,-2
          y(i)=y(i-2)-y(i)
14      continue
        y(2)=2.0*ytemp
        do 15 i=3,n,2
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
          y1=y(i)*wr+y(i+1)*wi
          y2=y(i+1)*wr-y(i)*wi
          y(i)=y1
          y(i+1)=y2
15      continue
        call realft(y,n,-1)
        do 16 i=1,n/2
          y1=y(i)+y(n-i+1)
          y2=(0.5/wi1)*(y(i)-y(n-i+1))
          y(i)=0.5*(y1+y2)
          y(n-i+1)=0.5*(y1-y2)
          wtemp=wr1
          wr1=wr1*wpr-wi1*wpi+wr1
          wi1=wi1*wpr+wtemp*wpi+wi1
16      continue
      endif
      return
      END
