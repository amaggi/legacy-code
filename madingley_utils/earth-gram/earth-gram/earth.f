c
c $Id: earth.f,v 1.1.1.1 2002/07/12 11:15:19 maggi Exp $
c $Log: earth.f,v $
c Revision 1.1.1.1  2002/07/12 11:15:19  maggi
c
c
c Revision 1.1  2002/05/23 10:28:29  maggi
c Initial revision
c
c
c*** locked mode program ***
c input file consists of :
c    n0,iefl,tref   # layers in model(including half space) and earth
c                  flattening control variable ( > 0 applies correction),
c                  reference period for material dispersion correction
c                  (0 for none)
c    d,vp,vs,rho,qbeta,qalpha   model,d is layer thickness. model can
c                  include a one layer ocean (signalled by setting vs=0
c                  in the top layer). half space can have any thickness
c                  assocoated with it ( 0 is ok)
c    jcom    =0 quits ;=1 rayleigh waves ; <> 1 love waves
c    output file name
c    c1,c2,nbran1,nbran2  min and max phase velocities and min and max
c                  branch numbers. note that c1=c2=0 causes program to
c                  choose phase velocity range for itself
c    nsrce,npts,dt  # source depths; # points in seismogram and its sample
c                  interval (in secs)
c    sdep          source depths
c    rdep          receiver depth
c
      program earthstuf
      implicit real*8(a-h,o-z)
      include '../include/sizes.inc'
      include '../include/units.inc'
      common/m/d(lyrs),ro(lyrs),vp(lyrs),vs(lyrs),fu(lyrs),n,noc,ist
      common/m0/d0(lyrs),ro0(lyrs),vp0(lyrs),vs0(lyrs),
     &           qb0(lyrs),qa0(lyrs),n0
      common/mq/vps(lyrs),vss(lyrs)
      common/q/qb(lyrs),qa(lyrs)
      common/bran/ce(2),ke(2),de(2),ctry,ceps,um,cm,cmx,cmn
      common/bits/u,nsrce,idep(lsd),nord,tpi,sdep(lsd),ig,idisc,irdep
      character*(256) infil
      istop=-1
      tpi=6.2831853071796d0
      do 5 i=1,256
    5 infil(i:i)=' '
      write(*,'(a$)') ' input file : '
      read(*,'(a256)') infil
      open(iinf1,file=infil,status='old')
	  
      write(*,'(a$)') ' output dispersion file : '
      read(*,'(a256)') infil
      open(iouf2,file=infil,status='unknown')
	  
      read(iinf1,*,end=777) n0,iefl,tref
      omref=0.d0
      if(tref.ne.0.d0) omref=tpi/tref
      do 10 i=1,n0
      read(iinf1,*) d0(i),vp0(i),vs0(i),ro0(i),qb0(i),qa0(i)
      if(qb0(i).ne.0.d0) qb0(i)=1.d0/qb0(i)
      if(qa0(i).ne.0.d0) qa0(i)=1.d0/qa0(i)
   10 continue
   45 read(iinf1,*)jcom
      if(jcom.eq.0)go to 777
      if(jcom.ne.1) jcom=2
      do 15 i=1,256
   15 infil(i:i)=' '
      read(iinf1,'(a256)') infil
      open(iouf1,file=infil,status='unknown',form='unformatted')
      read(iinf1,*) c1,c2,nbran1,nbran2
      n=n0
      do 20 i=1,n
      d(i)=d0(i)
      ro(i)=ro0(i)
      vps(i)=vp0(i)
      vss(i)=vs0(i)
      qa(i)=qa0(i)
   20 qb(i)=qb0(i)
      read(iinf1,*) nsrce,npts,dt
      npts=2*nfac(npts/2)
      read(iinf1,*) (sdep(i),i=1,nsrce)
      call shell(sdep,nsrce)
      read(iinf1,*) rdep
      call rsplit(rdep)
      call split(rdep)
	  write(iouf2,'(/,a,i3)') ' receiver index=',irdep
	  write(iouf2,'(/,a)') ' source indices=' 
	  do k=1,nsrce
      write(iouf2,'(2x,i3,1x,f10.4)') idep(k),sdep(k)
	  enddo

	  write(iouf2,'(/,a,1x,i3)') ' model ', n
	  do i=1,n
	  write(iouf2,'(4(1x,f7.3))') d(i),ro(i),vps(i),vss(i)
	  enddo
      call flat(jcom,iefl)
	  write(iouf2,'(/,a)') ' flattened model '
	  do i=1,n
	  write(iouf2,'(4(1x,f7.3))') d(i),ro(i),vps(i),vss(i)
	  enddo
      noc=1
      if(vss(1).le.0.d0) noc=2
      cmin=0.6d0*vss(noc)
      cmin=dmin1(cmin,vps(1))
      cmn=dmax1(c1,cmin)
      do 25 i=1,n
      vss(i)=vss(i)*vss(i)
   25 vps(i)=vps(i)*vps(i)
      trec=npts*dt
      dom=tpi/trec
      nh=npts/2
      omax=(nh+1)*dom
      om=omax-dom
      call qcor(om,omref)
      cmax=dsqrt(vs(n))-1.d-8
      ctst=c2
      if(ctst.le.0.d0) ctst=cmax
      cmx=dmin1(ctst,cmax)
      cmaxi=cmx
      call detk(cmaxi,om,kei,dei,jcom)
      if(nbran2.lt.0.or.nbran2.ge.kei) nbran2=kei-1
      nbran=max0((nbran2-nbran1+1),1)
      write(iouf2,910) nbran
  910 format(' number of modes to include : ',i5)
      write(iouf1)nsrce,npts,dt,jcom,nbran
      do 123 i=1,nsrce
      write(iouf1)sdep(i)
  123 continue
      nb=nbran1
   30 ctry=0.5d0*(cmn+cmx)
      ceps=0.5d0*(cmx-ctry)
      cm=0.d0
      write(iouf1) nb
      write(*,*) nb
      do 35 i=1,nh
      om=omax-i*dom
      call qcor(om,omref)
      cmax=dsqrt(vs(n))-1.d-8
      ctst=c2
      if(ctst.le.0.d0) ctst=cmax
      cmx=dmin1(ctst,cmax)
      call cex(om,nb,jcom,nev)
      if(nev.eq.0)go to 40
   35 call intrp(om,dom,jcom)
   40 write(iouf2,905) nb
      write(*,905) nb
  905 format(' ******** mode number : ',i5,' done. ********')
      nb=nb+1
      write(iouf1)istop
      write(*,*)istop
      if(nb.le.nbran2) goto 30
      close(iouf1)
      go to 45
  777 close(iinf1)
      close(iouf2)
      stop
      end

      subroutine split(rdep)
c  subroutine splits a layer at the source depth.the source
c  must be above the half space and not at the surface.
      implicit real*8(a-h,o-z)
      include '../include/sizes.inc'
      common/m/d(lyrs),ro(lyrs),vp(lyrs),vs(lyrs),fu(lyrs),n,noc,ist
      common/q/qb(lyrs),qa(lyrs)
      common/mq/vps(lyrs),vss(lyrs)
      common/bits/u,nsrce,idep(lsd),nord,tpi,sdep(lsd),ig,idisc,irdep
      do 10 ndp=1,nsrce
      thk=0.d0
      iss=0
      if(sdep(ndp).le.0.d0) go to 10
      do 1 is=1,n
      iss=iss+1
      thk=thk+d(is)
    1 if(sdep(ndp).le.thk)go to 2
    2 if(sdep(ndp).eq.thk)go to 10
c  if the source is already at an interface don't do anything
c  split the layer containing the source
      splt=thk-sdep(ndp)
c  shift the model down
      is1=iss+1
      nsf=n-iss
      do 5 l=1,nsf
      k=n+1-l
      j=k+1
      d(j)=d(k)
      ro(j)=ro(k)
      vps(j)=vps(k)
      vss(j)=vss(k)
      qa(j)=qa(k)
    5 qb(j)=qb(k)
c  split the source layer
      d(is1)=splt
      d(iss)=d(iss)-splt
      ro(is1)=ro(iss)
      vps(is1)=vps(iss)
      vss(is1)=vss(iss)
      qa(is1)=qa(iss)
      qb(is1)=qb(iss)
      n=n+1
      if(sdep(ndp).lt.rdep) irdep=irdep+1
   10 idep(ndp)=iss+1
      return
      end

      subroutine rsplit(rdep)
c  subroutine splits a layer at the receiver depth.
      implicit real*8(a-h,o-z)
      include '../include/sizes.inc'
      common/m/d(lyrs),ro(lyrs),vp(lyrs),vs(lyrs),fu(lyrs),n,noc,ist
      common/q/qb(lyrs),qa(lyrs)
      common/mq/vps(lyrs),vss(lyrs)
      common/bits/u,nsrce,idep(lsd),nord,tpi,sdep(lsd),ig,idisc,irdep

      if(rdep.le.0.d0) then
	    irdep=1
		return
	  endif

      thk=0.d0
      iss=0	  
      do 1 is=1,n
      iss=iss+1
      thk=thk+d(is)
    1 if(rdep.le.thk)go to 2
    2 if(rdep.eq.thk) then
c     if the receiver is already at an interface don't do anything
	    irdep=iss+1
		return
	  endif
	  irdep=iss+1
c  split the layer containing the source
      splt=thk-rdep
c  shift the model down
      is1=iss+1
      nsf=n-iss
      do 5 l=1,nsf
      k=n+1-l
      j=k+1
      d(j)=d(k)
      ro(j)=ro(k)
      vps(j)=vps(k)
      vss(j)=vss(k)
      qa(j)=qa(k)
    5 qb(j)=qb(k)
c  split the receiver layer
      d(is1)=splt
      d(iss)=d(iss)-splt
      ro(is1)=ro(iss)
      vps(is1)=vps(iss)
      vss(is1)=vss(iss)
      qa(is1)=qa(iss)
      qb(is1)=qb(iss)
      n=n+1
      return
      end
	  
      subroutine deriv(cc,w,ls)
c  deriv analytically calculates the rayleigh layer integrals required for the
c  group velocity and for the phase velocity derivatives.it is assumed
c  that uz,ur,tz,tr have been calculated by detray and stored in array x.
c  log derivatives are stored in array der as dc/drho,dc/dalf,dc/dbet,
c  spectra are computed at a distance of 1000km for source size 10.**27
c  cf mendiguren j.g.r. 1977 for excitation functions
      implicit real*8(a-h,o-z)
      include '../include/sizes.inc'
      include '../include/units.inc'
      common/x/x(4,lyrs),scale(lyrs),der(3,lyrs),dummy(lyrs*5)
      common/q/qb(lyrs),qa(lyrs)
      common/bits/u,nsrce,idep(lsd),nord,tpi,sdep(lsd),ig,idisc,irdep
      common/m/d(lyrs),ro(lyrs),vp2(lyrs),vs2(lyrs),fu(lyrs),n,noc,ist

      p=1.d0/cc
      psq=p*p
      ha=dsqrt(psq-1.d0/vp2(ls))
      hb=dsqrt(psq-1.d0/vs2(ls))
      c1=(hb*x(1,ls)+p*x(2,ls))/(psq-ha*hb)
      c2=(ha*x(2,ls)+p*x(1,ls))/(psq-ha*hb)
      c3=ro(ls)*c1*c2*(hb/vp2(ls)+ha/vs2(ls))/(p*(ha+hb))
      t1=ro(ls)*(ha*c1*c1+hb*c2*c2-2.d0*p*c1*c2)
      t2=0.5d0*ro(ls)*c1*c1/(vp2(ls)*ha)
      t3=0.5d0*ro(ls)*c2*c2/(vs2(ls)*hb)+t2
c  si1,si2,si3 are the energy integrals
      si1=t1+t3
      si2=4.d0*vs2(ls)*psq*(t1+c3)+t3
      si3=0.5d0*(si2+si1)
      der(1,ls)=0.5d0*(si2-si1)
      der(2,ls)=t2
      der(3,ls)=si2-t2
      i=ls
  100 i1=i
      i=i-1
      r2=2.d0*fu(i)*p
      e4=r2*x(1,i1)-x(4,i1)
      e1=ro(i)*x(1,i1)-p*e4
      e2=r2*x(2,i1)-x(3,i1)
      e3=ro(i)*x(2,i1)-p*e2
      f4=r2*x(1,i)-x(4,i)
      f1=ro(i)*x(1,i)-p*f4
      f2=r2*x(2,i)-x(3,i)
      f3=ro(i)*x(2,i)-p*f2
      dh=0.5d0*d(i)*w
      ha=psq-1.d0/vp2(i)
      hb=psq-1.d0/vs2(i)
      c1=dh*(e1*e1-ha*e2*e2)
      c2=0.5d0*(e1*e2-f1*f2)
      c3=dh*(e3*e3-hb*e4*e4)
      c4=0.5d0*(e3*e4-f3*f4)
      c5=p*(e2*e4-f2*f4)
      c6=cc*(e1*e3-f1*f3)-c5
      t1=2.d0*(c5+c2+c4)/ro(i)
      t2=(c2-c1)/(ro(i)*vp2(i)*ha)
      t3=(c4-c3)/(fu(i)*hb)+t2
      sj1=t1+t3
      sj2=4.d0*vs2(i)*psq*(t1+c6/ro(i))+t3
      si1=si1+sj1
      si2=si2+sj2
      si3=si3+0.5d0*(sj2+sj1)-(c1+c3)/ro(i)
      der(1,i)=0.5d0*(sj2-sj1)
      der(2,i)=t2
      der(3,i)=sj2-t2
      if(i.gt.noc) goto 100
      if(noc.eq.1) goto 30
c  energy integrals and partials in ocean layer
      ha=(psq-1.d0/vp2(1))/ro(1)
      c1=0.5d0*d(1)*w*(x(1,noc)*x(1,noc)-ha*x(3,noc)*x(3,noc)/ro(1))
      c2=-0.5d0*x(1,noc)*x(3,noc)/ro(1)
      t2=(c2-c1)/(ha*vp2(1))
      t1=2.d0*ro(1)*c2+t2
      si1=si1+t1
      si2=si2+t2
      si3=si3+0.5d0*(t1+t2)-ro(1)*c1
      der(1,1)=0.5d0*(t2-t1)
      der(2,1)=t2
      der(3,1)=0.d0
   30 u=cc*si3/si1
c  flan should be almost zero if the eigenvector is accurate
      flan=si1/si2-1.d0
      if(dabs(flan).ge.1.e-8) goto 999
      fnorm=1.d0/si3
      do 40 i=1,ls
      do 40 j=1,3
   40 der(j,i)=der(j,i)*fnorm
c    ignore attenuation in the ocean
      q=0.d0
      do 50 i=noc,ls
   50 q=q+der(2,i)*qa(i)+der(3,i)*qb(i)
      gam=0.5d0*w*q/cc
 	  gb2=w*(-x(1,irdep)/cc + x(4,irdep)/fu(irdep))
	  gb1=w*( (1.d0-2.d0*vs2(irdep)/vp2(irdep))*x(2,irdep)/cc 
     1          + x(3,irdep)/(ro(irdep)*vp2(irdep)) )
      write(iouf1) w,cc,gam,x(1,irdep),x(2,irdep),gb1,gb2
c  write out the excitations.py1,py2,py3 correspond to a,b,c in mendiguren
      fact=p*dsqrt(p*w)/(si3*15.853309d-6)
      do 60 i=1,nsrce
      id=idep(i)
      py1=x(2,id)*p*fact
      py2=x(4,id)*fact/fu(id)
      sig=ro(id)*vp2(id)
      py3=-(x(3,id)/sig+p*(1.d0-2.d0*fu(id)/sig)*x(2,id))*fact
      if(id.le.ls) goto 60
      py1=0.d0
      py2=0.d0
      py3=0.d0
   60 write(iouf1) py1,py2,py3
      per=tpi/w
      if(q.ne.0.d0)q=1.d0/q
      write(iouf2,900) nord,per,cc,u,q,flan
  900 format(1x,i5,5g15.7)
      return
  999 write(*,950) flan
  950 format(' problem with eigenfunction : flan =',g15.7)
      stop
      end

      subroutine detlov(cc,w,de,ifeif)
c  computes the love stress-displacement vector and propagates it
c  upwards. if ifeif=0 only the determinant at the surface
c  is returned. if ifeif=1 derivatives,excitations etc. are computed
      implicit real*8(a-h,o-z)
      include '../include/sizes.inc'
      include '../include/units.inc'
      common/x/x(4,lyrs),scale(lyrs),der(3,lyrs),dummy(lyrs*5)
      common/m/d(lyrs),ro(lyrs),vp2(lyrs),vs2(lyrs),fu(lyrs),n,noc,ist
      common/bits/u,nsrce,idep(lsd),nord,tpi,sdep(lsd),ig,idisc,irdep
      common/q/qb(lyrs),qa(lyrs)
      p=1.d0/cc
      psq=p*p
      x(1,ist)=1.d0
      x(2,ist)=-fu(ist)*dsqrt(psq-1.d0/vs2(ist))
      i=ist
      scale(i)=0.d0
c  propagate the solution up
  100 i1=i
      i=i-1
      f1=1.d0/fu(i)
      hb=psq-1.d0/vs2(i)
      wd=w*d(i)
      call arg(wd,hb,cb,sb,fac)
      scale(i)=scale(i1)+fac
      fsb=f1*sb
      hsb=hb*sb*fu(i)
      x(1,i)=cb*x(1,i1)+fsb*x(2,i1)
      x(2,i)=hsb*x(1,i1)+cb*x(2,i1)
      if(i.gt.noc) goto 100
      de=x(2,noc)/(dabs(x(2,noc))+dabs(x(1,noc)))
      if(ifeif.eq.0) return
c  compute derivatives,excitations etc.
      xnorm=1.d0/x(1,noc)
      knt=noc-1
      do 25 i=noc,ist
      knt=knt+1
      xfac=1.d0
      if(scale(i)-scale(noc).ne.0.d0) xfac=dexp(scale(i)-scale(noc))
      ls=i
      x(1,i)=x(1,i)*xnorm*xfac
      x(2,i)=x(2,i)*xnorm*xfac
      if(dabs(de).gt.1.d-4) goto 25
      if(i.lt.noc+1) goto 25
c  if stress-displacement vector is small enough and solution is no
c  longer oscillatory reduce the model size
      pbsq=1.d0/vs2(i)
      if(xfac.lt.1.d-15.and.psq.ge.pbsq) goto 30
   25 continue
   30 if(dabs(de).ge.1.d-4)call fixlov(w,p)
      hb=w*dsqrt(psq-1.d0/vs2(ls))
      fi1=0.5d0*x(1,ls)*x(1,ls)/hb
      fi2=0.5d0*x(2,ls)*x(2,ls)/(hb*fu(ls))
c  si1,si2,si3 are the energy integrals
c  der(1,i),der(2,i)are the log. derivatives of phase velocity wrt rho and vs
      si1=ro(ls)*fi1
      si2=psq*fu(ls)*fi1+fi2
      si3=fu(ls)*fi1
      der(1,ls)=0.5d0*(si2-si1)
      der(2,ls)=si2
      i=knt
   35 i1=i
      i=i-1
      c0=psq*fu(i)-ro(i)
      c1=0.5d0*(x(1,i)*x(2,i)-x(1,i1)*x(2,i1))/w
      c2=0.5d0*d(i)*(c0*x(1,i1)*x(1,i1)-x(2,i1)*x(2,i1)/fu(i))
      fi1=(c2-c1)/c0
      fi2=-c2-c1
      sj1=ro(i)*fi1
      sj2=psq*fu(i)*fi1+fi2
      si1=si1+sj1
      si2=si2+sj2
      si3=si3+fu(i)*fi1
      der(1,i)=0.5d0*(sj2-sj1)
      der(2,i)=sj2
      if(i.gt.noc) goto 35
      u=p*si3/si1
c  flan should be close to zero if the eigenfunction is accurate
      flan=si2/si1-1.d0
      if(dabs(flan).ge.1.e-8) goto 999
      fnorm=1.d0/(psq*si3)
      q=0.d0
      do 40 i=noc,ls
      der(1,i)=der(1,i)*fnorm
      der(2,i)=der(2,i)*fnorm
   40 q=q+der(2,i)*qb(i)
      gam=0.5d0*w*q/cc
c  compute 'excitation' functions.py1,py2 correspond to a,b in mendiguren
      gb1=w*x(2,irdep)/fu(irdep)      
      write(iouf1) w,cc,gam,x(1,irdep),gb1
      fact=1.d0/(dsqrt(p*w)*si3*15.853309d-6)
      do 45 i=1,nsrce
      id=idep(i)
      py1=x(1,id)*p*fact
      py2=x(2,id)*fact/fu(id)
      if(id.le.ls) goto 45
      py1=0.d0
      py2=0.d0
   45 write(iouf1) py1,py2
      per=tpi/w
      if(q.ne.0.d0)q=1.d0/q
      write(iouf2,900) nord,per,cc,u,q,flan
  900 format(1x,i5,5g15.7)
      return
  999 write(*,950) flan
  950 format(' problem with eigenfunction : flan =',g15.7)
      stop
      end

      subroutine intrp(om,dom,jcom)
c  interpolates between bracketing c's to find the root.
c  uses a bisection scheme.
      implicit real*8(a-h,o-z)
      include '../include/sizes.inc'
      common/m/d(lyrs),ro(lyrs),vp2(lyrs),vs2(lyrs),fu(lyrs),n,noc,ist
      common/bits/u,nsrce,idep(lsd),nord,tpi,sdep(lsd),ig,idisc,irdep
      common/bran/ce(2),ke(2),de(2),ctry,ceps,um,cm,cmx,cmn
      data tol,dlta/1.d-11,1.d-5/

      fc=de(1)
      fb=de(2)
      if(fc*fb.ge.0.d0) return
      nord=ke(1)
      c=ce(1)
      b=ce(2)
      psq=1.d0/(b*b)
      call strtdp(psq,om)
      s=c
      fs=fc
c*** bisect ***
    5 h=0.5d0*(b+c)
      t=h*tol
c*** check for convergence ***
      db=dabs(fb)
      dc=dabs(fc)
      if(dabs(h-b).lt.t) goto 35
      if(db.le.dc) goto 10
      y=b
      fy=fb
      gg=b
      fg=fb
      s=c
      fs=fc
      goto 15
   10 y=s
      fy=fs
      gg=c
      fg=fc
      s=b
      fs=fb
   15 if(fy.eq.fs) goto 20
      b=(s*fy-y*fs)/(fy-fs)
      if(dabs(b-s).lt.t) b=s+dsign(t,gg-s)
      if((b-h)*(s-b).lt.0.d0) b=h
      goto 25
   20 b=h
   25 call detray(b,om,fb,0,jcom)
      if(fg*fb.lt.0.d0) goto 30
      c=s
      fc=fs
      goto 5
   30 c=gg
      fc=fg
      goto 5
   35 if(dc.lt.db) b=c
      call detray(b,om,fb,1,jcom)
      cp=b*(1.d0-b/u)/om
      ctry=0.d0
      if(cm.le.0.d0) goto 60
      clin=b-cp*dom
      ctry=5.d0*cm-4.d0*b-2.d0*dom*(um+2.d0*cp)
      ceps=dmax1(dabs(clin-ctry),ctry*tol)
   60 um=cp
      cm=b
      if(ctry.ne.0.d0) return
      ctry=b-cp*dom
      ceps=dmax1(dabs(cp*dom),ctry*tol)
      return
      end

      subroutine detray(cc,w,de,ifeif,jcom)
c  routine computes minor vector in each layer for
c  raleigh waves by progagation upwards. if ifeif=0
c  only the determinant at the surface is returned.
c  if ifeif=1 the stress-displacement vector is
c   computed in each layer and deriv is called.
      implicit real*8(a-h,o-z)
      include '../include/sizes.inc'
      common/x/x(4,lyrs),scale(lyrs),der(3,lyrs),dummy(lyrs*5)
      common/bits/u,nsrce,idep(lsd),nord,tpi,sdep(lsd),ig,idisc,irdep
      common/m/d(lyrs),ro(lyrs),vp2(lyrs),vs2(lyrs),fu(lyrs),n,noc,ist
      common/y/y1,y2,y3,y4,y5
      dimension ym(5,lyrs),y(5)
      equivalence (ym,x),(y1,y)
      if(jcom.eq.1) goto 1
      call detlov(cc,w,de,ifeif)
      return
c  compute minor vector y at bottom
    1 p=1.d0/cc
      ysav=0.d0
      psq=p*p
      r2=2.d0*fu(ist)*p
      y3=dsqrt(psq-1.d0/vp2(ist))
      y4=-dsqrt(psq-1.d0/vs2(ist))
      y1=-(y3*y4+psq)/ro(ist)
      y2=r2*y1+p
      y5=ro(ist)-r2*(p+y2)
      i=ist
      scale(i)=0.d0
      do 5 j=1,5
    5 ym(j,i)=y(j)
c**** propagate up layers ****
  100 i=i-1
      wd=w*d(i)
      ha=psq-1.d0/vp2(i)
      call arg(wd,ha,ca,sa,faca)
      hb=psq-1.d0/vs2(i)
      call arg(wd,hb,cb,sb,facb)
      scale(i)=scale(i+1)+faca+facb
      hbs=hb*sb
      has=ha*sa
      r1=1.d0/ro(i)
      r2=2.d0*fu(i)*p
      b1=r2*y1-y2
      g3=(y5+r2*(y2-b1))*r1
      g1=b1+p*g3
      g2=ro(i)*y1-p*(g1+b1)
      g1=g1*dexp(-faca-facb)
      e1=cb*g2-hbs*y3
      e2=-sb*g2+cb*y3
      e3=cb*y4+hbs*g3
      e4=sb*y4+cb*g3
      y3=ca*e2-has*e4
      y4=sa*e1+ca*e3
      g3=ca*e4-sa*e2
      b1=g1-p*g3
      y1=(ca*e1+has*e3+p*(g1+b1))*r1
      y2=r2*y1-b1
      y5=ro(i)*g3-r2*(y2-b1)
      do 10 j=1,5
   10 ym(j,i)=y(j)
      if(i.gt.noc) goto 100
      if(noc.eq.1) goto 15
c  *** propagate through ocean layer ***
      ha=psq-1.d0/vp2(1)
      wd=w*d(1)
      call arg(wd,ha,ca,sa,faca)
      y1=ca*y3-ha*sa*y5/ro(1)
      y5=ca*y5-ro(1)*sa*y3
      y2=y5
   15 de=y5/dsqrt(y1*y1+y2*y2)
      if(ifeif.eq.0)return
c  compute stress-displacement vector x=n*y
      ynorm=1.d0/ym(3,noc)
      if(noc.eq.1) goto 20
c  cope with possibility of stoneley mode on ocean floor by integrating
c  from surface
      y1=ca
      y2=ro(1)*sa
      ysav=y2/y1
      de1=de
      de=dmin1(dabs(de1),dabs(ynorm*ym(5,noc)/ysav-1.d0))
c  specify arbitrary solution y at the surface
   20 y1=0.d0
      y2=-ynorm
      y3=0.d0
      y4=0.d0
      xfac=1.d0
      sum=0.d0
      i=noc
c  minor  elements compose matrix n so that x=n*b.(x is the stress-displacement
c  vector) compute x in each layer. this happens to be numerically stable.
c  b is a solution to the equations of motion db/dz=ab
   25 xx1=-ym(2,i)*y1-ym(3,i)*y2+ym(1,i)*y4
      xx2=-ym(4,i)*y1+ym(2,i)*y2-ym(1,i)*y3
      xx3=-ym(5,i)*y2-ym(2,i)*y3-ym(4,i)*y4
      xx4= ym(5,i)*y1-ym(3,i)*y3+ym(2,i)*y4
      x(1,i)=xx1*xfac
      x(2,i)=xx2*xfac
      x(3,i)=xx3*xfac
      x(4,i)=xx4*xfac
      ls=i
      if(i.eq.ist) goto 30
      if(i.lt.noc+1) goto 35
      if(dabs(de).gt.1.d-4) goto 35
c  if x becomes  small and the solution is no longer oscillatory
c  reduce the model size
      pbsq=1.d0/vs2(i)
      if(xfac.lt.1.d-15.and.psq.ge.pbsq) go to 30
   35 wd=w*d(i)
      ha=psq-1.d0/vp2(i)
      call arg(wd,ha,ca,sa,faca)
      hb=psq-1.d0/vs2(i)
      call arg(wd,hb,cb,sb,facb)
      dfac=dexp(facb-faca)
      cb=dfac*cb
      sb=dfac*sb
      hbs=hb*sb
      has=ha*sa
      r2=2.d0*p*fu(i)
      e2=r2*y2-y3
      e3=ro(i)*y2-p*e2
      e4=r2*y1-y4
      e1=ro(i)*y1-p*e4
      e6=ca*e2-sa*e1
      e8=cb*e4-sb*e3
      y1=(ca*e1-has*e2+p*e8)/ro(i)
      y2=(cb*e3-hbs*e4+p*e6)/ro(i)
      y3=r2*y2-e6
      y4=r2*y1-e8
      i=i+1
      sum=sum+faca
      xfac=dexp(scale(i)-scale(noc)+sum)
      goto 25
   30 if(dabs(de).gt.1.d-4) call fixray(ysav,w,p)
c  compute excitations or partial derivatives
   50 call deriv(cc,w,ls)
      return
      end

      subroutine qcor(om,omref)
      implicit real*8(a-h,o-z)
      include '../include/sizes.inc'
      common/mq/vps(lyrs),vss(lyrs)
      common/m/d(lyrs),ro(lyrs),vp(lyrs),vs(lyrs),fu(lyrs),n,noc,ist
      common/q/qb(lyrs),qa(lyrs)
      data pii/.3183098d0/
      omscl=0.d0
      if(omref.ne.0.d0)omscl=pii*dlog(om/omref)
      do 1 i=1,n
      vs(i)=vss(i)*(1.d0+qb(i)*omscl)**2
      vp(i)=vps(i)*(1.d0+qa(i)*omscl)**2
    1 fu(i)=ro(i)*vs(i)
      return
      end
