      SUBROUTINE frenel(x,s,c)
      INTEGER MAXIT
      REAL c,s,x,EPS,FPMIN,PI,PIBY2,XMIN
      PARAMETER (EPS=6.e-8,MAXIT=100,FPMIN=1.e-30,XMIN=1.5,PI=3.1415927,
     *PIBY2=1.5707963)
      INTEGER k,n
      REAL a,absc,ax,fact,pix2,sign,sum,sumc,sums,term,test
      COMPLEX b,cc,d,h,del,cs
      LOGICAL odd
      absc(h)=abs(real(h))+abs(aimag(h))
      ax=abs(x)
      if(ax.lt.sqrt(FPMIN))then
        s=0.
        c=ax
      else if(ax.le.XMIN)then
        sum=0.
        sums=0.
        sumc=ax
        sign=1.
        fact=PIBY2*ax*ax
        odd=.true.
        term=ax
        n=3
        do 11 k=1,MAXIT
          term=term*fact/k
          sum=sum+sign*term/n
          test=abs(sum)*EPS
          if(odd)then
            sign=-sign
            sums=sum
            sum=sumc
          else
            sumc=sum
            sum=sums
          endif
          if(term.lt.test)goto 1
          odd=.not.odd
          n=n+2
11      continue
        pause 'series failed in frenel'
1       s=sums
        c=sumc
      else
        pix2=PI*ax*ax
        b=cmplx(1.,-pix2)
        cc=1./FPMIN
        d=1./b
        h=d
        n=-1
        do 12 k=2,MAXIT
          n=n+2
          a=-n*(n+1)
          b=b+4.
          d=1./(a*d+b)
          cc=b+a/cc
          del=cc*d
          h=h*del
          if(absc(del-1.).lt.EPS)goto 2
12      continue
        pause 'cf failed in frenel'
2       h=h*cmplx(ax,-ax)
        cs=cmplx(.5,.5)*(1.-cmplx(cos(.5*pix2),sin(.5*pix2))*h)
        c=real(cs)
        s=aimag(cs)
      endif
      if(x.lt.0.)then
        c=-c
        s=-s
      endif
      return
      END
