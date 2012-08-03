      SUBROUTINE spread(y,yy,n,x,m)
      INTEGER m,n
      REAL x,y,yy(n)
      INTEGER ihi,ilo,ix,j,nden,nfac(10)
      REAL fac
      SAVE nfac
      DATA nfac /1,1,2,6,24,120,720,5040,40320,362880/
      if(m.gt.10) pause 'factorial table too small in spread'
      ix=x
      if(x.eq.float(ix))then
        yy(ix)=yy(ix)+y
      else
        ilo=min(max(int(x-0.5*m+1.0),1),n-m+1)
        ihi=ilo+m-1
        nden=nfac(m)
        fac=x-ilo
        do 11 j=ilo+1,ihi
          fac=fac*(x-j)
11      continue
        yy(ihi)=yy(ihi)+y*fac/(nden*(x-ihi))
        do 12 j=ihi-1,ilo,-1
          nden=(nden/(j+1-ilo))*(j-ihi)
          yy(j)=yy(j)+y*fac/(nden*(x-j))
12      continue
      endif
      return
      END
