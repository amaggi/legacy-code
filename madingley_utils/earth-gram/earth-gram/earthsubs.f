c
c $Id: earthsubs.f,v 1.1.1.1 2002/07/12 11:15:19 maggi Exp $
c $Log: earthsubs.f,v $
c Revision 1.1.1.1  2002/07/12 11:15:19  maggi
c
c
c Revision 1.1  2002/05/23 10:28:29  maggi
c Initial revision
c
c
      subroutine detk(cc,w,kount,de,jcom)
c  routine computes rayleigh minor vector x or y and
c  propagates it from the bottom upwards. keeps track
c  of the determinant (x(5)) and increments the counter
c  when a branch is crossed
      implicit real*8(a-h,o-z)
      include '../include/sizes.inc'
      common/y/y1,y2,y3,y4,y5
      common/m/d(lyrs),ro(lyrs),vp2(lyrs),vs2(lyrs),fu(lyrs),n,noc,ist
      dimension y(5),x(5)
      equivalence (y,y1)
      data tpi/6.2831853071796d0/
 
      if(jcom.eq.1) goto 1
      call lovknt(cc,w,kount,de)
      return
    1 kount=0
      p=1.d0/cc
      psq=p*p
c  compute the start depth
      call strtdp(psq,w)
c  compute the starting solution
      r2=2.d0*fu(ist)*p
      y3=dsqrt(psq-1.d0/vp2(ist))
      y4=-dsqrt(psq-1.d0/vs2(ist))
      y1=-(y3*y4+psq)/ro(ist)
      y2=r2*y1+p
      y5=ro(ist)-r2*(p+y2)
      if(y5.gt.0.d0) kount=1
      i=ist
c**** propagate up layers ****
  100 i=i-1
      r1=1.d0/ro(i)
      r2=2.d0*fu(i)*p
      b1=r2*y1-y2
      g3=(y5+r2*(y2-b1))*r1
      g1=b1+p*g3
      g2=ro(i)*y1-p*(g1+b1)
      ha=psq-1.d0/vp2(i)
      hb=psq-1.d0/vs2(i)
      sum=0.d0
      dels=0.d0
      delm=d(i)
      if(hb.lt.0.d0) delm=tpi*.1d0/(w*dsqrt(-hb))
      if(ha.lt.0.d0) delm=tpi*.1d0/(w*(dsqrt(-hb)+dsqrt(-ha)))
   49 del=delm
      isplit=0
   50 if(sum+del.gt.d(i)) del=d(i)-sum
      if(dels.eq.del) goto 5
      dels=del
      wd=w*del
      call arg(wd,ha,ca,sa,faca)
      call arg(wd,hb,cb,sb,facb)
      hbs=hb*sb
      has=ha*sa
    5 g1=g1*dexp(-faca-facb)
      e1=cb*g2-hbs*y3
      e2=-sb*g2+cb*y3
      e3=cb*y4+hbs*g3
      e4=sb*y4+cb*g3
      x(3)=ca*e2-has*e4
      x(4)=sa*e1+ca*e3
      g2=ca*e1+has*e3
      g3=ca*e4-sa*e2
      b1=g1-p*g3
      x(1)=(g2+p*(b1+g1))*r1
      x(2)=r2*x(1)-b1
      x(5)=ro(i)*g3-r2*(x(2)-b1)
      if(x(5)*y(5).le.0.d0) goto 11
      x5p=y5-wd*(ro(i)*(y4-y3)-4.d0*fu(i)*psq*y4*(1.d0-vs2(i)/vp2(i)))
      if(x(5)*x5p.lt.0.d0)go to 12
      go to 15
   11 t1=y(3)-y(4)
      t2=x(3)-x(4)
      if(t1*t2.gt.0.d0) goto 10
   12 del=0.5d0*del
      isplit=isplit+1
      b1=r2*y1-y2
      g3=(y5+r2*(y2-b1))*r1
      g1=b1+p*g3
      g2=ro(i)*y1-p*(g1+b1)
      if(isplit.lt.100)go to 50
      write(*,998) i
  998 format(' count problem in layer ',i5)
      stop
   10 tes=t1*(y(5)-x(5))
      if(tes.lt.0.d0) kount=kount+1
      if(tes.gt.0.d0) kount=kount-1
   15 do 20 j=1,5
   20 y(j)=x(j)
      sum=sum+del
      if(sum.lt.d(i)) goto 49
      if(i.gt.noc) goto 100
      if(noc.eq.1) goto 25
c     propagate through ocean layer
      y1=y3
      y2=y5
      r1=1.d0/ro(1)
      ha=psq-1.d0/vp2(1)
      sum=0.d0
      dels=0.d0
      del=d(1)
      if(ha.lt.0.d0) del=tpi*0.1d0/(w*dsqrt(-ha))
   30 if(sum+del.gt.d(1))del=d(1)-sum
      if(dels.eq.del) goto 35
      dels=del
      wd=w*del
      call arg(wd,ha,ca,sa,faca)
      rsa=ro(1)*sa
      hsa=ha*sa*r1
   35 x1=ca*y1-hsa*y2
      x2=ca*y2-rsa*y1
      if(x2*y2.gt.0.d0) goto 45
      if(x1*y1.gt.0.d0) goto 40
      del=0.5d0*del
      goto 30
   40 tes=x1*(y2-x2)
      if(tes.lt.0.d0) kount=kount+1
      if(tes.gt.0.d0) kount=kount-1
   45 y1=x1
      y2=x2
      sum=sum+del
      if(sum.lt.d(1)) goto 30
      y5=y2
   25 de=y5/dsqrt(y1*y1+y2*y2)
      return
      end

      subroutine arg(wd,h,c,s,fac)
      implicit real*8(a-h,o-z)
      hh=dsqrt(dabs(h))
      th=wd*hh
      if(th.lt.1.5d-14) goto 10
      if(h.gt.0.d0) goto 5
      c=dcos(th)
      s=-dsin(th)/hh
      fac=0.d0
      return
    5 d=dexp(-2.d0*th)
      c=0.5d0*(1.d0+d)
      s=-0.5d0*(1.d0-d)/hh
      fac=th
      return
   10 c=1.d0
      s=-wd
      fac=0.d0
      return
      end

      subroutine lovknt(cc,w,kount,de)
c  compute the love stress-displacement vector at the bottom
c  and propagate up. determines mode count.
      implicit real*8(a-h,o-z)
      include '../include/sizes.inc'
      include '../include/units.inc'
      common/m/d(lyrs),ro(lyrs),vp2(lyrs),vs2(lyrs),fu(lyrs),n,noc,ist
      data pi/3.14159265358979d0/
      kount=0
      p=1.d0/cc
      psq=p*p
      call strtdp(psq,w)
      y1=1.d0
      y2=-fu(ist)*dsqrt(psq-1.d0/vs2(ist))
      i=ist
  100 i=i-1
      hb=psq-1.d0/vs2(i)
      if(hb.ge.0.d0) goto 25
      hh=dsqrt(-hb)
      wh=w*hh
      a1=fu(i)*hb*y1
      c1=0.5d0*pi/wh
      if(a1.ne.0.d0) c1=datan(y2*hh/a1)/wh
      pih=pi/wh
      k=-1
   10 k=k+1
      xtry=c1+k*pih
      if(xtry.gt.d(i)) goto 12
      if(xtry.le.0.d0) goto 10
      kount=kount+1
      goto 10
   12 th=wh*d(i)
      cb=dcos(th)
      sb=-dsin(th)/hh
      x1=cb*y1+sb*y2/fu(i)
      x2=hb*sb*fu(i)*y1 +cb*y2
      goto 15
   25 wd=w*d(i)
      call arg(wd,hb,cb,sb,fac)
      x1=cb*y1+sb*y2/fu(i)
      x2=hb*sb*fu(i)*y1+cb*y2
      if(x2*y2.lt.0.d0) kount=kount-1
   15 y1=x1
      y2=x2
      if(i.gt.noc) goto 100
      de=y2/(dabs(y2)+dabs(y1))
      return
      end

      subroutine cex(om,nb,jcom,nev)
c  makes sure branch nb is isolated given ctry+-ceps as trial c
      implicit real*8(a-h,o-z)
      include '../include/sizes.inc'
      common/bran/ce(2),ke(2),de(2),ctry,ceps,um,cm,cmx,cmn
      data tol/2.d-9/
c*** bracket this mode ***
      nev=0
      cepm=20.d0*ceps
      nup=nb+1
      cx=ctry
      ce(1)=0.5d0*cmn
      ce(2)=2.d0*cmx
      ke(1)=10000
      ke(2)=-10000
      i=0
    5 if(ke(2)-ke(1).eq.1) goto 50
        i=i+1
        cx=dmin1(dmax1(cx,cmn),cmx)
        call detk(cx,om,kx,dx,jcom)
        if(kx.ge.nup) goto 30
        if(cx.ge.cmx) return
        if(kx.ne.nb.or.ke(1).ne.nb) goto 10
        delc=-2.d0*dx*(cx-ce(1))/(dx-de(1))
        ceps=dmax1(dmin1(cepm,delc),ceps)
   10   ce(1)=cx
        ke(1)=kx
        de(1)=dx
        cx=cx+ceps
        if(cx.ge.ce(2)) cx=0.5d0*(ce(1)+ce(2))
        goto 5
   30   if(cx.le.cmn) return
        if(kx.ne.nup.or.ke(2).ne.nup) goto 35
        delc=2.d0*dx*(cx-ce(2))/(dx-de(2))
        ceps=dmax1(dmin1(cepm,delc),ceps)
   35   ce(2)=cx
        ke(2)=kx
        de(2)=dx
        cx=cx-ceps
        if(cx.le.ce(1)) cx=0.5d0*(ce(1)+ce(2))
        goto 5
   50   continue
      nev=1
      return
      end

      subroutine flat(jcom,iefl)
c  subroutine applies earth-flattening corrections by
c  scaling the model.for reference see biswas 1972.
      implicit real*8(a-h,o-z)
      include '../include/sizes.inc'
      common/bits/u,nsrce,idep(lsd),nord,tpi,sdep(lsd),ig,idisc,irdep
      common/m/h(lyrs),ro(lyrs),vp(lyrs),vs(lyrs),z(lyrs),n,noc,ist
      common/mq/vps(lyrs),vss(lyrs)
      common/files/iofl
      data a/6371.d0/
      if(iefl.eq.0)return
      pwr=2.275d0
      if(iabs(jcom).gt.1)pwr=5.d0
      atp=a**pwr
c  flatten the  model
      nm=n-1
      hs=0.d0
      do 3 i=1,n
      ht=hs
      hs=hs+h(i)
    3 h(i)=a-ht
      do 13 i=1,nm
      ii=i+1
      fltd=dlog(h(i)/h(ii))
      dif=(1.d0/h(ii)-1.d0/h(i))*a/fltd
      difr=h(i)**pwr-h(ii)**pwr
      ro(i)=ro(i)*difr/(fltd*(a**pwr)*pwr)
      vps(i)=vps(i)*dif
   13 vss(i)=vss(i)*dif
c half space scaling
      fact=a/h(n)
      facti=(h(n)**pwr)/atp
      vps(n)=vps(n)*fact
      vss(n)=vss(n)*fact
      ro(n)=ro(n)*facti
      z0=0.d0
      do 14 i=2,n
      z1=a*dlog(a/h(i))
      h(i-1)=z1-z0
   14 z0=z1
      h(n)=0.d0
      return
      end

      subroutine strtdp(psq,w)
c  routine computes the start depth where the solution
c  has decayed by 10e-20 from the turning point depth
      implicit real*8(a-h,o-z)
      include '../include/sizes.inc'
      common/m/d(lyrs),ro(lyrs),vp2(lyrs),vs2(lyrs),fu(lyrs),n,noc,ist
      nm1=n-1
      do 5 k=noc,n
      j=n+noc-k
    5 if(psq-1.d0/vs2(j).le.0.d0) goto 10
   10 if(j.lt.nm1) goto 15
      ist=n
      return
   15 test=0.d0
      j=j+1
      knt=j-1
      do 20 k=j,nm1
      knt=knt+1
      hb=dsqrt(psq-1.d0/vs2(k))
      test=test+d(k)*hb*w
   20 if(test.gt.40.d0) goto 25
   25 ist=knt+1
      return
      end

      subroutine shell(a,ncol)
      implicit real*8(a-h,o-z)
      include '../include/sizes.inc'
c  shell sorts a(nc) into increasing order
      dimension a(ncol)
      ig=ncol
    5 if(ig.le.1) return
      ig=ig/2
      im=ncol-ig
   10 iex=0
      do 20 i=1,im
      ipg=i+ig
      if(a(i).le.a(ipg)) go to 20
      te=a(i)
      a(i)=a(ipg)
      a(ipg)=te
      iex=iex+1
   20 continue
      if(iex.gt.0) go to 10
      go to 5
      end

      subroutine fixray(ysav,w,p)
      implicit real*8(a-h,o-z)
      include '../include/sizes.inc'
      common/x/x(4,lyrs),scale(lyrs),der(3,lyrs),dummy(lyrs*5)
      common/m/d(lyrs),ro(lyrs),vp2(lyrs),vs2(lyrs),fu(lyrs),n,noc,ist
      common/fix/b(2,4)
      dimension xx(4)
      psq=p*p
      lvz=ist+1
      do 30 i=noc,ist
      lvz=lvz-1
      hb=psq-1.d0/vs2(lvz)
      if(hb.lt.0.d0) goto 35
   30 continue
   35 write(*,999) lvz
  999 format(' **fixray entered-redo propagation to layer ',i4,'**')
      lvz=max0(noc,lvz-1)
      do 40 i=1,4
   40 xx(i)=x(i,lvz)
      do 1 i=1,2
      do 1 j=1,4
    1 b(i,j)=0.d0
      b(1,1)=1.d0
      b(2,2)=1.d0
      b(1,3)=ysav
      cosav=0.d0
      i=noc
    2 ii=i+1
      r2=2.d0*fu(i)*p
      wd=w*d(i)
      ha=psq-1.d0/vp2(i)
      hb=psq-1.d0/vs2(i)
      call arg(wd,ha,ca,sa,faca)
      call arg(wd,hb,cb,sb,facb)
      dfac=dexp(faca)
      sa=sa*dfac
      ca=ca*dfac
      dfac=dexp(facb)
      sb=sb*dfac
      cb=cb*dfac
      do 3 k=1,2
      e2=r2*b(k,2)-b(k,3)
      e3=ro(i)*b(k,2)-p*e2
      e4=r2*b(k,1)-b(k,4)
      e1=ro(i)*b(k,1)-p*e4
      e6=ca*e2-sa*e1
      e8=cb*e4-sb*e3
      b(k,1)=(ca*e1-ha*sa*e2+p*e8)/ro(i)
      b(k,2)=(cb*e3-hb*sb*e4+p*e6)/ro(i)
      b(k,3)=r2*b(k,2)-e6
    3 b(k,4)=r2*b(k,1)-e8
      call solve(ii,coef)
      dco=dabs((coef-cosav)/coef)
      cosav=coef
      i=ii
      if(i.eq.lvz)go to 99
      if(dco.gt.1.d-9) goto 2
   99 x(1,noc)=1.d0
      x(2,noc)=coef
      x(3,noc)=ysav
      x(4,noc)=0.d0
      i=noc
    4 ii=i+1
      r2=2.d0*fu(i)*p
      wd=w*d(i)
      ha=psq-1.d0/vp2(i)
      hb=psq-1.d0/vs2(i)
      call arg(wd,ha,ca,sa,faca)
      call arg(wd,hb,cb,sb,facb)
      dfac=dexp(faca)
      sa=sa*dfac
      ca=ca*dfac
      dfac=dexp(facb)
      sb=sb*dfac
      cb=cb*dfac
      e2=r2*x(2,i)-x(3,i)
      e3=ro(i)*x(2,i)-p*e2
      e4=r2*x(1,i)-x(4,i)
      e1=ro(i)*x(1,i)-p*e4
      e6=ca*e2-sa*e1
      e8=cb*e4-sb*e3
      x(1,ii)=(ca*e1-ha*sa*e2+p*e8)/ro(i)
      x(2,ii)=(cb*e3-hb*sb*e4+p*e6)/ro(i)
      x(3,ii)=r2*x(2,ii)-e6
      x(4,ii)=r2*x(1,ii)-e8
      i=ii
      if(i.lt.lvz) goto 4
      xtxx=0.d0
      xtx=0.d0
      do 5 j=1,4
      xtxx=xtxx+x(j,lvz)*xx(j)
    5 xtx=xtx+xx(j)*xx(j)
      scl=xtxx/xtx
      lvz1=lvz+1
      do 8 i=lvz1,ist
      do 8 j=1,4
    8 x(j,i)=scl*x(j,i)
      return
      end
	  
      subroutine solve(ii,coef)
      implicit real*8(a-h,o-z)
      include '../include/sizes.inc'
      common/x/x(4,lyrs),scale(lyrs),der(3,lyrs),dummy(lyrs*5)
      common/fix/b(2,4)
      xx=0.d0
      b22=0.d0
      do 5 k=1,4
      xx=xx+x(k,ii)*x(k,ii)
    5 b22=b22+b(2,k)*b(2,k)
      f=dsqrt(b22/xx)
      xb1=0.d0
      xb2=0.d0
      b12=0.d0
      do 10 k=1,4
      xb1=xb1+x(k,ii)*b(1,k)*f
      xb2=xb2+x(k,ii)*b(2,k)*f
   10 b12=b12+b(1,k)*b(2,k)
      det=b22*b22-xb2*xb2
      coef=(xb2*xb1-b22*b12)/det
      return
      end

      subroutine fixlov(w,p)
      implicit real*8(a-h,o-z)
      include '../include/sizes.inc'
      common/x/x(4,lyrs),scale(lyrs),der(3,lyrs),dummy(lyrs*5)
      common/m/d(lyrs),ro(lyrs),vp2(lyrs),vs2(lyrs),fu(lyrs),n,noc,ist
      psq=p*p
      lvz=ist+1
      do 30 i=noc,ist
      lvz=lvz-1
      hb=psq-1.d0/vs2(lvz)
      if(hb.lt.0.d0) goto 35
   30 continue
   35 write(*,999) lvz
  999 format(' **fixlov entered-redo propagation to layer ',i4,'**')
      lvz=max0(noc,lvz-1)
      x1s=x(1,lvz)
      x2s=x(2,lvz)
      x(1,noc)=1.d0
      x(2,noc)=0.d0
      i=noc
    2 ii=i+1
      hb=psq-1.d0/vs2(i)
      wd=w*d(i)
      call arg(wd,hb,cb,sb,facb)
      dfac=dexp(facb)
      sb=sb*dfac
      cb=cb*dfac
      x(1,ii)=cb*x(1,i)-sb*x(2,i)/fu(i)
      x(2,ii)=cb*x(2,i)-hb*sb*fu(i)*x(1,i)
      i=ii
      if(i.lt.lvz) goto 2
      scl=(x(1,lvz)*x1s+x(2,lvz)*x2s)/(x1s*x1s+x2s*x2s)
      lvz1=lvz+1
      do 4 i=lvz1,ist
      do 4 j=1,2
    4 x(j,i)=scl*x(j,i)
      return
      end
