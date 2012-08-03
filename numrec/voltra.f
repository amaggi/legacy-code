      SUBROUTINE voltra(n,m,t0,h,t,f,g,ak)
      INTEGER m,n,MMAX
      REAL h,t0,f(m,n),t(n),g,ak
      EXTERNAL ak,g
      PARAMETER (MMAX=5)
CU    USES ak,g,lubksb,ludcmp
      INTEGER i,j,k,l,indx(MMAX)
      REAL d,sum,a(MMAX,MMAX),b(MMAX)
      t(1)=t0
      do 11 k=1,m
        f(k,1)=g(k,t(1))
11    continue
      do 16 i=2,n
        t(i)=t(i-1)+h
        do 14 k=1,m
          sum=g(k,t(i))
          do 13 l=1,m
            sum=sum+0.5*h*ak(k,l,t(i),t(1))*f(l,1)
            do 12 j=2,i-1
              sum=sum+h*ak(k,l,t(i),t(j))*f(l,j)
12          continue
            if(k.eq.l)then
              a(k,l)=1.
            else
              a(k,l)=0.
            endif
            a(k,l)=a(k,l)-0.5*h*ak(k,l,t(i),t(i))
13        continue
          b(k)=sum
14      continue
        call ludcmp(a,m,MMAX,indx,d)
        call lubksb(a,m,MMAX,indx,b)
        do 15 k=1,m
          f(k,i)=b(k)
15      continue
16    continue
      return
      END
