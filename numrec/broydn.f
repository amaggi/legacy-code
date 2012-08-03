      SUBROUTINE broydn(x,n,check)
      INTEGER n,nn,NP,MAXITS
      REAL x(n),fvec,EPS,TOLF,TOLMIN,TOLX,STPMX
      LOGICAL check
      PARAMETER (NP=40,MAXITS=200,EPS=1.e-7,TOLF=1.e-4,TOLMIN=1.e-6,
     *TOLX=EPS,STPMX=100.)
      COMMON /newtv/ fvec(NP),nn
CU    USES fdjac,fmin,lnsrch,qrdcmp,qrupdt,rsolv
      INTEGER i,its,j,k
      REAL den,f,fold,stpmax,sum,temp,test,c(NP),d(NP),fvcold(NP),g(NP),
     *p(NP),qt(NP,NP),r(NP,NP),s(NP),t(NP),w(NP),xold(NP),fmin
      LOGICAL restrt,sing,skip
      EXTERNAL fmin
      nn=n
      f=fmin(x)
      test=0.
      do 11 i=1,n
        if(abs(fvec(i)).gt.test)test=abs(fvec(i))
11    continue
      if(test.lt..01*TOLF)return
      sum=0.
      do 12 i=1,n
        sum=sum+x(i)**2
12    continue
      stpmax=STPMX*max(sqrt(sum),float(n))
      restrt=.true.
      do 44 its=1,MAXITS
        if(restrt)then
          call fdjac(n,x,fvec,NP,r)
          call qrdcmp(r,n,NP,c,d,sing)
          if(sing) pause 'singular Jacobian in broydn'
          do 14 i=1,n
            do 13 j=1,n
              qt(i,j)=0.
13          continue
            qt(i,i)=1.
14        continue
          do 18 k=1,n-1
            if(c(k).ne.0.)then
              do 17 j=1,n
                sum=0.
                do 15 i=k,n
                  sum=sum+r(i,k)*qt(i,j)
15              continue
                sum=sum/c(k)
                do 16 i=k,n
                  qt(i,j)=qt(i,j)-sum*r(i,k)
16              continue
17            continue
            endif
18        continue
          do 21 i=1,n
            r(i,i)=d(i)
            do 19 j=1,i-1
              r(i,j)=0.
19          continue
21        continue
        else
          do 22 i=1,n
            s(i)=x(i)-xold(i)
22        continue
          do 24 i=1,n
            sum=0.
            do 23 j=i,n
              sum=sum+r(i,j)*s(j)
23          continue
            t(i)=sum
24        continue
          skip=.true.
          do 26 i=1,n
            sum=0.
            do 25 j=1,n
              sum=sum+qt(j,i)*t(j)
25          continue
            w(i)=fvec(i)-fvcold(i)-sum
            if(abs(w(i)).ge.EPS*(abs(fvec(i))+abs(fvcold(i))))then
              skip=.false.
            else
              w(i)=0.
            endif
26        continue
          if(.not.skip)then
            do 28 i=1,n
              sum=0.
              do 27 j=1,n
                sum=sum+qt(i,j)*w(j)
27            continue
              t(i)=sum
28          continue
            den=0.
            do 29 i=1,n
              den=den+s(i)**2
29          continue
            do 31 i=1,n
              s(i)=s(i)/den
31          continue
            call qrupdt(r,qt,n,NP,t,s)
            do 32 i=1,n
              if(r(i,i).eq.0.) pause 'r singular in broydn'
              d(i)=r(i,i)
32          continue
          endif
        endif
        do 34 i=1,n
          sum=0.
          do 33 j=1,n
            sum=sum+qt(i,j)*fvec(j)
33        continue
          g(i)=sum
34      continue
        do 36 i=n,1,-1
          sum=0.
          do 35 j=1,i
            sum=sum+r(j,i)*g(j)
35        continue
          g(i)=sum
36      continue
        do 37 i=1,n
          xold(i)=x(i)
          fvcold(i)=fvec(i)
37      continue
        fold=f
        do 39 i=1,n
          sum=0.
          do 38 j=1,n
            sum=sum+qt(i,j)*fvec(j)
38        continue
          p(i)=-sum
39      continue
        call rsolv(r,n,NP,d,p)
        call lnsrch(n,xold,fold,g,p,x,f,stpmax,check,fmin)
        test=0.
        do 41 i=1,n
          if(abs(fvec(i)).gt.test)test=abs(fvec(i))
41      continue
        if(test.lt.TOLF)then
          check=.false.
          return
        endif
        if(check)then
          if(restrt)then
            return
          else
            test=0.
            den=max(f,.5*n)
            do 42 i=1,n
              temp=abs(g(i))*max(abs(x(i)),1.)/den
              if(temp.gt.test)test=temp
42          continue
            if(test.lt.TOLMIN)then
              return
            else
              restrt=.true.
            endif
          endif
        else
          restrt=.false.
          test=0.
          do 43 i=1,n
            temp=(abs(x(i)-xold(i)))/max(abs(x(i)),1.)
            if(temp.gt.test)test=temp
43        continue
          if(test.lt.TOLX)return
        endif
44    continue
      pause 'MAXITS exceeded in broydn'
      END
