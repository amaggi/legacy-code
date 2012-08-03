      SUBROUTINE toeplz(r,x,y,n)
      INTEGER n,NMAX
      REAL r(2*n-1),x(n),y(n)
      PARAMETER (NMAX=100)
      INTEGER j,k,m,m1,m2
      REAL pp,pt1,pt2,qq,qt1,qt2,sd,sgd,sgn,shn,sxn,g(NMAX),h(NMAX)
      if(r(n).eq.0.) goto 99
      x(1)=y(1)/r(n)
      if(n.eq.1)return
      g(1)=r(n-1)/r(n)
      h(1)=r(n+1)/r(n)
      do 15 m=1,n
        m1=m+1
        sxn=-y(m1)
        sd=-r(n)
        do 11 j=1,m
          sxn=sxn+r(n+m1-j)*x(j)
          sd=sd+r(n+m1-j)*g(m-j+1)
11      continue
        if(sd.eq.0.)goto 99
        x(m1)=sxn/sd
        do 12 j=1,m
          x(j)=x(j)-x(m1)*g(m-j+1)
12      continue
        if(m1.eq.n)return
        sgn=-r(n-m1)
        shn=-r(n+m1)
        sgd=-r(n)
        do 13 j=1,m
          sgn=sgn+r(n+j-m1)*g(j)
          shn=shn+r(n+m1-j)*h(j)
          sgd=sgd+r(n+j-m1)*h(m-j+1)
13      continue
        if(sd.eq.0..or.sgd.eq.0.)goto 99
        g(m1)=sgn/sgd
        h(m1)=shn/sd
        k=m
        m2=(m+1)/2
        pp=g(m1)
        qq=h(m1)
        do 14 j=1,m2
          pt1=g(j)
          pt2=g(k)
          qt1=h(j)
          qt2=h(k)
          g(j)=pt1-pp*qt2
          g(k)=pt2-pp*qt1
          h(j)=qt1-qq*pt2
          h(k)=qt2-qq*pt1
          k=k-1
14      continue
15    continue
      pause 'never get here in toeplz'
99    pause 'singular principal minor in toeplz'
      END
