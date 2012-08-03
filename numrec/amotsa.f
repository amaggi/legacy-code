      FUNCTION amotsa(p,y,psum,mp,np,ndim,pb,yb,funk,ihi,yhi,fac)
      INTEGER ihi,mp,ndim,np,NMAX
      REAL amotsa,fac,yb,yhi,p(mp,np),pb(np),psum(np),y(mp),funk
      PARAMETER (NMAX=200)
      EXTERNAL funk
CU    USES funk,ran1
      INTEGER idum,j
      REAL fac1,fac2,tt,yflu,ytry,ptry(NMAX),ran1
      COMMON /ambsa/ tt,idum
      fac1=(1.-fac)/ndim
      fac2=fac1-fac
      do 11 j=1,ndim
        ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
11    continue
      ytry=funk(ptry)
      if (ytry.le.yb) then
        do 12 j=1,ndim
          pb(j)=ptry(j)
12      continue
        yb=ytry
      endif
      yflu=ytry-tt*log(ran1(idum))
      if (yflu.lt.yhi) then
        y(ihi)=ytry
        yhi=yflu
        do 13 j=1,ndim
          psum(j)=psum(j)-p(ihi,j)+ptry(j)
          p(ihi,j)=ptry(j)
13      continue
      endif
      amotsa=yflu
      return
      END
