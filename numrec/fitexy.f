      SUBROUTINE fitexy(x,y,ndat,sigx,sigy,a,b,siga,sigb,chi2,q)
      INTEGER ndat,NMAX
      REAL x(ndat),y(ndat),sigx(ndat),sigy(ndat),a,b,siga,sigb,chi2,q,
     *POTN,PI,BIG,ACC
      PARAMETER (NMAX=1000,POTN=1.571000,BIG=1.e30,PI=3.14159265,
     *ACC=1.e-3)
CU    USES avevar,brent,chixy,fit,gammq,mnbrak,zbrent
      INTEGER j,nn
      REAL xx(NMAX),yy(NMAX),sx(NMAX),sy(NMAX),ww(NMAX),swap,amx,amn,
     *varx,vary,aa,offs,ang(6),ch(6),scale,bmn,bmx,d1,d2,r2,dum1,dum2,
     *dum3,dum4,dum5,brent,chixy,gammq,zbrent
      COMMON /fitxyc/ xx,yy,sx,sy,ww,aa,offs,nn
      EXTERNAL chixy
      if (ndat.gt.NMAX) pause 'NMAX too small in fitexy'
      call avevar(x,ndat,dum1,varx)
      call avevar(y,ndat,dum1,vary)
      scale=sqrt(varx/vary)
      nn=ndat
      do 11 j=1,ndat
        xx(j)=x(j)
        yy(j)=y(j)*scale
        sx(j)=sigx(j)
        sy(j)=sigy(j)*scale
        ww(j)=sqrt(sx(j)**2+sy(j)**2)
11    continue
      call fit(xx,yy,nn,ww,1,dum1,b,dum2,dum3,dum4,dum5)
      offs=0.
      ang(1)=0.
      ang(2)=atan(b)
      ang(4)=0.
      ang(5)=ang(2)
      ang(6)=POTN
      do 12 j=4,6
        ch(j)=chixy(ang(j))
12    continue
      call mnbrak(ang(1),ang(2),ang(3),ch(1),ch(2),ch(3),chixy)
      chi2=brent(ang(1),ang(2),ang(3),chixy,ACC,b)
      chi2=chixy(b)
      a=aa
      q=gammq(0.5*(nn-2),0.5*chi2)
      r2=0.
      do 13 j=1,nn
        r2=r2+ww(j)
13    continue
      r2=1./r2
      bmx=BIG
      bmn=BIG
      offs=chi2+1.
      do 14 j=1,6
        if (ch(j).gt.offs) then
          d1=mod(abs(ang(j)-b),PI)
          d2=PI-d1
          if(ang(j).lt.b)then
            swap=d1
            d1=d2
            d2=swap
          endif
          if (d1.lt.bmx) bmx=d1
          if (d2.lt.bmn) bmn=d2
        endif
14    continue
      if (bmx.lt. BIG) then
        bmx=zbrent(chixy,b,b+bmx,ACC)-b
        amx=aa-a
        bmn=zbrent(chixy,b,b-bmn,ACC)-b
        amn=aa-a
        sigb=sqrt(0.5*(bmx**2+bmn**2))/(scale*cos(b)**2)
        siga=sqrt(0.5*(amx**2+amn**2)+r2)/scale
      else
        sigb=BIG
        siga=BIG
      endif
      a=a/scale
      b=tan(b)/scale
      return
      END
