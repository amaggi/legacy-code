      SUBROUTINE dfpmin(p,n,gtol,iter,fret,func,dfunc)
      INTEGER iter,n,NMAX,ITMAX
      REAL fret,gtol,p(n),func,EPS,STPMX,TOLX
      PARAMETER (NMAX=50,ITMAX=200,STPMX=100.,EPS=3.e-8,TOLX=4.*EPS)
      EXTERNAL dfunc,func
CU    USES dfunc,func,lnsrch
      INTEGER i,its,j
      LOGICAL check
      REAL den,fac,fad,fae,fp,stpmax,sum,sumdg,sumxi,temp,test,dg(NMAX),
     *g(NMAX),hdg(NMAX),hessin(NMAX,NMAX),pnew(NMAX),xi(NMAX)
      fp=func(p)
      call dfunc(p,g)
      sum=0.
      do 12 i=1,n
        do 11 j=1,n
          hessin(i,j)=0.
11      continue
        hessin(i,i)=1.
        xi(i)=-g(i)
        sum=sum+p(i)**2
12    continue
      stpmax=STPMX*max(sqrt(sum),float(n))
      do 27 its=1,ITMAX
        iter=its
        call lnsrch(n,p,fp,g,xi,pnew,fret,stpmax,check,func)
        fp=fret
        do 13 i=1,n
          xi(i)=pnew(i)-p(i)
          p(i)=pnew(i)
13      continue
        test=0.
        do 14 i=1,n
          temp=abs(xi(i))/max(abs(p(i)),1.)
          if(temp.gt.test)test=temp
14      continue
        if(test.lt.TOLX)return
        do 15 i=1,n
          dg(i)=g(i)
15      continue
        call dfunc(p,g)
        test=0.
        den=max(fret,1.)
        do 16 i=1,n
          temp=abs(g(i))*max(abs(p(i)),1.)/den
          if(temp.gt.test)test=temp
16      continue
        if(test.lt.gtol)return
        do 17 i=1,n
          dg(i)=g(i)-dg(i)
17      continue
        do 19 i=1,n
          hdg(i)=0.
          do 18 j=1,n
            hdg(i)=hdg(i)+hessin(i,j)*dg(j)
18        continue
19      continue
        fac=0.
        fae=0.
        sumdg=0.
        sumxi=0.
        do 21 i=1,n
          fac=fac+dg(i)*xi(i)
          fae=fae+dg(i)*hdg(i)
          sumdg=sumdg+dg(i)**2
          sumxi=sumxi+xi(i)**2
21      continue
        if(fac**2.gt.EPS*sumdg*sumxi)then
          fac=1./fac
          fad=1./fae
          do 22 i=1,n
            dg(i)=fac*xi(i)-fad*hdg(i)
22        continue
          do 24 i=1,n
            do 23 j=1,n
              hessin(i,j)=hessin(i,j)+fac*xi(i)*xi(j)-fad*hdg(i)*hdg(j)+
     *fae*dg(i)*dg(j)
23          continue
24        continue
        endif
        do 26 i=1,n
          xi(i)=0.
          do 25 j=1,n
            xi(i)=xi(i)-hessin(i,j)*g(j)
25        continue
26      continue
27    continue
      pause 'too many iterations in dfpmin'
      return
      END
