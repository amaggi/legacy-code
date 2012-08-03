      FUNCTION chixy(bang)
      REAL chixy,bang,BIG
      INTEGER NMAX
      PARAMETER (NMAX=1000,BIG=1.E30)
      INTEGER nn,j
      REAL xx(NMAX),yy(NMAX),sx(NMAX),sy(NMAX),ww(NMAX),aa,offs,avex,
     *avey,sumw,b
      COMMON /fitxyc/ xx,yy,sx,sy,ww,aa,offs,nn
      b=tan(bang)
      avex=0.
      avey=0.
      sumw=0.
      do 11 j=1,nn
        ww(j)=(b*sx(j))**2+sy(j)**2
        if(ww(j).eq.0.) then
          ww(j)=BIG
        else
          ww(j)=1./ww(j)
        endif
        sumw=sumw+ww(j)
        avex=avex+ww(j)*xx(j)
        avey=avey+ww(j)*yy(j)
11    continue
      avex=avex/sumw
      avey=avey/sumw
      aa=avey-b*avex
      chixy=-offs
      do 12 j=1,nn
        chixy=chixy+ww(j)*(yy(j)-aa-b*xx(j))**2
12    continue
      return
      END
