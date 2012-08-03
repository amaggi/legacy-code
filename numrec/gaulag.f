      SUBROUTINE gaulag(x,w,n,alf)
      INTEGER n,MAXIT
      REAL alf,w(n),x(n)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=3.D-14,MAXIT=10)
CU    USES gammln
      INTEGER i,its,j
      REAL ai,gammln
      DOUBLE PRECISION p1,p2,p3,pp,z,z1
      do 13 i=1,n
        if(i.eq.1)then
          z=(1.+alf)*(3.+.92*alf)/(1.+2.4*n+1.8*alf)
        else if(i.eq.2)then
          z=z+(15.+6.25*alf)/(1.+.9*alf+2.5*n)
        else
          ai=i-2
          z=z+((1.+2.55*ai)/(1.9*ai)+1.26*ai*alf/(1.+3.5*ai))*
     *(z-x(i-2))/(1.+.3*alf)
        endif
        do 12 its=1,MAXIT
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j
11        continue
          pp=(n*p1-(n+alf)*p2)/z
          z1=z
          z=z1-p1/pp
          if(abs(z-z1).le.EPS)goto 1
12      continue
        pause 'too many iterations in gaulag'
1       x(i)=z
        w(i)=-exp(gammln(alf+n)-gammln(float(n)))/(pp*n*p2)
13    continue
      return
      END
