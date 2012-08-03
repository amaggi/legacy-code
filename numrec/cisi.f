      SUBROUTINE cisi(x,ci,si)
      INTEGER MAXIT
      REAL ci,si,x,EPS,EULER,PIBY2,FPMIN,TMIN
      PARAMETER (EPS=6.e-8,EULER=.57721566,MAXIT=100,PIBY2=1.5707963,
     *FPMIN=1.e-30,TMIN=2.)
      INTEGER i,k
      REAL a,err,fact,sign,sum,sumc,sums,t,term,absc
      COMPLEX h,b,c,d,del
      LOGICAL odd
      absc(h)=abs(real(h))+abs(aimag(h))
      t=abs(x)
      if(t.eq.0.)then
        si=0.
        ci=-1./FPMIN
        return
      endif
      if(t.gt.TMIN)then
        b=cmplx(1.,t)
        c=1./FPMIN
        d=1./b
        h=d
        do 11 i=2,MAXIT
          a=-(i-1)**2
          b=b+2.
          d=1./(a*d+b)
          c=b+a/c
          del=c*d
          h=h*del
          if(absc(del-1.).lt.EPS)goto 1
11      continue
        pause 'cf failed in cisi'
1       continue
        h=cmplx(cos(t),-sin(t))*h
        ci=-real(h)
        si=PIBY2+aimag(h)
      else
        if(t.lt.sqrt(FPMIN))then
          sumc=0.
          sums=t
        else
          sum=0.
          sums=0.
          sumc=0.
          sign=1.
          fact=1.
          odd=.true.
          do 12 k=1,MAXIT
            fact=fact*t/k
            term=fact/k
            sum=sum+sign*term
            err=term/abs(sum)
            if(odd)then
              sign=-sign
              sums=sum
              sum=sumc
            else
              sumc=sum
              sum=sums
            endif
            if(err.lt.EPS)goto 2
            odd=.not.odd
12        continue
          pause 'maxits exceeded in cisi'
        endif
2       si=sums
        ci=sumc+log(t)+EULER
      endif
      if(x.lt.0.)si=-si
      return
      END
