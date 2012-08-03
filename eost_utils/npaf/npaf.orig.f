c---programme de filtrage  variable -version sun
c--  lecture des donnees GEOSCOPE en format ah
c---------------------------------------------------------------------
      character cod*6,chan*6,infile*50,outfile*50,outfile1*50,filevit*50
        include '/donnee/user0/ipg/stutzman/XDR/common_xdrhead'
      common/c1/y(33000),yy(33000),yf(33000),yref(33000)
      common /sleep/ u(3,54),t(54),x(54)
      common/tire/l1(30),l2(30),is(30),il(30)
      common/cogen/in,iout,iperf,dnu,pas
      common/codisp/ nq1,nq2
      common/ond/a(3,55),b(3,55),c(3,55),e(3,55)
      dimension aa(55),bb(55),cc(55),ee(55),uu(55),dy(55)
      dimension ipr(20)
	dimension tt1(30,55),uu1(30,55),ta1(55),ua1(55)
	integer nbper(30),nutrain(30),nmode(30)
	character*2 b1,b2
      real*8 r(440)
      pi=3.14159265
      ndim=32768
      luin=1
      luout=10
 1      continue
        infile='      '  
	outfile='     '
	npo=0
	n=0
	ntrain=0
	nk=0
	tdeb=0
	tfin=0.
	ipi=0
	ips=0
	do i=1,54
		t(i)=0.
		do k=1,3
			u(k,i)=0.
		enddo
	enddo
	
c----------------------------------------------------------------
c     lecture d un fichier ah
c---------------------------------------------------------------
      write(*,*)'name of the infile'
      read(*,'(a50)')infile
      write(6,'(a50)')infile
	if(infile.eq.'fin') goto 7777
      write(*,*)'nombre de points a lire'
      read(*,*)npo
      write(*,*)npo
      write(*,*)'name of the outfile'
      read(*,'(a50)')outfile1
      write(*,'(a50)')outfile1
c ntot= nb de mode   ntrain=nbre de train a traiter
      read(*,*)ntot,ntrain
      write(6,*)'ntot,ntrain',ntot,ntrain
      read (5,*) tdeb,tfin,ipi,ips
      write (6,*) tdeb,tfin,ipi,ips
	write(6,*) 'nom du fichier de vitesse'
	read(*,'(a20)')filevit
	write(6,'(a20)')filevit
	open(19,file=filevit,status='old')
	do  k=1,ntot
 	 	read(19,*) b1,b2,nbper(k)
		print*,b1,b2,nbper(k)
		read(b1(2:2),'(i1)') nutrain(k)
		read(b2(2:2),'(i1)') nmode(k)
      		read (19,*) (tt1(k,i),i=1,nbper(k))
      		write(6,*) (tt1(k,i),i=1,nbper(k))
		read (19,*) (uu1(k,i),i=1,nbper(k))
		write(6,*) (uu1(k,i),i=1,nbper(k))
	enddo
	close (19)
	do 109 kkk=0,4
		nk=0
		n=0
		do i=1,55
			ta1(i)=0.
			ua1(i)=0.
		enddo
		do i=1,ntot
			if (nmode(i).eq.kkk) then
				nk=nk+1
				if (nk.eq.1) then
					n=nbper(i)
					do j=1,n
						ta1(j)=tt1(i,j)
						ua1(j)=uu1(i,j)
					enddo
				else
				    k=1
			  	    do j=1,n
				    do jj=1,nbper(i)
					if(ta1(j).eq.tt1(i,jj)) then
				 	    t(k)=ta1(j)
					    u(1,k)=ua1(j)
					    u(2,k)=uu1(i,jj)
					    k=k+1
					endif
				    enddo
				    enddo	
				    n=k-1
				endif
				
			endif

		enddo
		if (nk.eq.1) then
			do j=1,n
				t(j)=ta1(j)
				u(1,j)=ua1(j)
			enddo
		endif
		if (nk.eq.0) then
			print*,'pas de courbe de dispersion pour n=',kkk
			goto 109 
		endif
		write(outfile,'(a1,i1,a28)') 'n',kkk,outfile1
        cod='  '
        chan='  '
        alat=0.
        alon=0.
        sdep=0.
        iy=0
        imo=0
        ida=0
        ihd=0
        imd=0
        sd=0.
        elate=0.
        elonge=0.
        edep=0.
        la=0
        lmo=0
        lj=0
        ih0=0
        im0=0
        s0=0.
        pas=0.
        npoint=0

	do  i=1,ndim
		y(i)=0.
		yy(i)=0.
		yref(i)=0.
		yf(i)=0.
	enddo

        luin=1
        mode=0
        call xdrf_ahopen(luin,infile,mode,ierr)
       if (ierr.ne.0) then
                write(6,*)'pb ouverture fichier de lecture'
     *                      ,luin,infile,mode,ierr
                goto 99
        endif
        call xdrf_rdahrecor(luin,cod,chan,alat,alon,sdep,la,
     *       lmo,lj,ihd,imd,sd,elate,elonge,edep,iy,imo,ida,ih0,
     *       im0,s0,pas,npoint,ierr,y)
        if(ierr.le.0) then
                write(6,*) 'pb lecture des donnee ',ierr
                goto 99

        endif
        call xdrf_ahclose(luin,ierr)
        if(ierr.ne.0) then
                write(*,*) 'code fermeture du fichier ',ierr
                goto 99
        endif

        print*,'fichier ',infile,' code et canal de la station ',cod,chan
        print*,'coordonnees de la station',alat,alon,sdep
        print*,'coordonnees de l epicentre',elate,elonge,edep
        print*,'date debut enregistrement',la,lmo,lj,ihd,imd,sd
        print*,'date seisme',iy,imo,ida,ih0,im0,s0
        print*,'nb points',npoint,' dt',pas,' npo=',npo
	write(*,*)'n  number of periods= ', n
	write(*,*)'nb de trains,nb de courbes de dispersion=' ,ntrain,nk
	write(*,*)'tdeb,tfin,ipi,ips(ipi,ips plateau)=' ,tdeb,tfin,ipi,ips
	write(*,*)' u(k,i):valeurs des points des courbes de dispersion'
	write (*,8888) ((u(k,i),i=1,n),k=1,nk)
	write(6,*) 'fichier de sortie:',outfile
	if (npo.gt.npoint) npo=npoint
	if (npo.lt.npoint) then
		do ik=npo+1,npoint
			y(ik)=0.
		enddo
	endif
	npoint=npo
		
      call julien1(ida,imo,iy,iy0,id,1)
      write(*,*)'day,month,year,year,julian day=',ida,imo,iy,iy0,id
c----------------------------------------------------------------
c     the epicentral distance is calculated (ellipticity correction)
      call great(elate,elonge,alat,alon,dist1,disd,gc,fai1,fai2)
      dist=dist1
      write(*,*) 'epicentral distance=',dist
c------------------------------------------------------------------
      s=(lj-ida)*24.*3600.+(ihd-ih0)*3600.+(imd-im0)*60.+sd-s0
      if (s) 40,40,41
   40 nst = s/pas-0.5
      go to 43
   41 nst= s/pas + 0.5
   43 write (*,*)'trecord -tevent,number of points  ',s,nst
      npoint=nst+npoint
      nst=nst+1

      npr=n-1
 8888 format (15f8.3)
      s=n
      sn=s/5.
      do 32 i=1,npr
   32 dy(i)=0.0001
      npr=npr+1
      do 36 i=npr,n
   36 dy(i)=0.0005
c  boucle sur les nk trains
      do 30 k=1,nk
c  boucle sur les n periodes
      do 31 i=1,n
      j=n-i+1
      x(i)=1./t(j)
   31 uu(i)=1./u(k,j)
   50 call lisse(n,x,uu,dy,s,aa,bb,cc,ee,r)
      print *,'la apres lisse'
      if(aa(1).ne.0.) go to 51
      write (1,161) s
  161 format ('  coefficients nuls  s= ',f10.2)
      s=s+sn
      isi=isi+1
      if(isi-10) 50,50,99
   51 s=n
      isi=0
      do 34 i=1,n
      a(k,i)=aa(i)
      b(k,i)=bb(i)
      c(k,i)=cc(i)
   34 e(k,i)=ee(i)
      print*,' resultat du lissage'
      write (*,160) (a(k,i),i=1,n)
      write (*,160) (b(k,i),i=1,n)
      write (*,160) (c(k,i),i=1,n)
      write (*,160) (e(k,i),i=1,n)
  160 format (8e15.6)
   30 continue
      write(*,*)'epicentral distance=',dist
      ip=1
    2 np=2**ip
      if(np-npoint) 3,4,4
    3 ip=ip+1
      go to 2
    4 dnu=1./(np*pas)
      nfrq=np/2.+1.
      write(*,102) npoint,np,ip,nfrq
      write (*,8889) pas,dnu
 8889 format (' pas en secondes=',f15.5,2x,' pas en frequence=',e15.8)
  102 format(8i10)
  103 format (8f10.5)
      print*,'enter szper'
      print*,'tdeb,tfin,ipi,ips,n',tdeb,tfin,ipi,ips,n
      call szper (tdeb,tfin,ipi,ips,n)
      print*,'szper passed'
      ideb =1
      if(nst.le.0) go to 208
      ideb=nst
      ymoy=0.
      do 211 i=1,100
  211 ymoy=ymoy + y(i)
      ymoy= ymoy/100.
      do 210 i=1,ideb
  210 yref(i) = ymoy
  208 do 201 i=ideb,npoint
      k=i-nst+1
  201 yref(i)=y(k)
      do 207 i=1,npoint
  207 y(i)=yref(i)
      np1=npoint+1
      do 209 i=np1,ndim
  209 y(i)=0.
      nst=1
      call sztr(y,nst,npoint)
      write (*,5313) nst,npoint
 5313 format(2x,'nst=',i6,'npoint=',i6)
      do 8 i=1,np
    8 yref(i)=y(i)
      call nlogn(ip,y,yy,+1.)
  777 format (2f10.3)
      j1=nq1
      j2=nq2
      write (1,5312) j1,j2
 5312 format (10i10)
c     apodisation du spectre a b.f.
      jm=2*j1
      fact=pi/2./j1
      do 90 j=j1,jm
      pond = sin(fact*(j-j1))
      y(j)=y(j)*pond
      yy(j) = yy(j)*pond
   90 continue
c     apodisation du spectre a hf
      jm=0.9*j2
      fact=5.*pi/j2
      do 92 j=jm,j2
      pond=sin(fact*(j2-j))
      y(j)=y(j)*pond
   92 yy(j)=yy(j)*pond
      jp=0
      print*,'nq1,nq2',nq1,nq2
      do 6 j=j1,j2
c     print*,'j1 j2 n ntrain',j1,j2,n,ntrain
      call szond (j,dist,n,ntrain,i1,i2,irec,nk)
c     print*,'szond passed'
      if(i1.gt.0)go to 223
      write (1,1000) j,i1,i2
      i1=1
  223 if(i2.le.npoint)go to 222
      write (*,1000) j,i1,i2
      i2=npoint
  222 continue

      jp=jp+1
      ipr(jp) = i1
      if (jp-20) 70,71,71
c     write (*,310) (ipr(jt),jt=1,20)
   71 continue
  310 format (1h ,20i5)
      jp=0
   70 continue
      pul=2*pi*(j-1)*dnu
      dphi=pul*pas
      dc=cos(dphi)
      ds=sin(dphi)
      if(irec) 99,13,12
   13 do 26 k=1,ntrain
      i1=l1(k)
      i2=l2(k)
      if(i1.le.0)i1=1
      if(i2.gt.npoint)i2=npoint
      phi=(i1-2)*dphi
      cosi=cos(phi)
      sini=sin(phi)
      do 25 i=i1,i2
      cosr=cosi
      cosi=cosi*dc-sini*ds
      sini=sini*dc+cosr*ds
      at=i-is(k)
      at=at/il(k)
      at=at**4
      apo=1.-at
   25 yf(i)=yf(i)+apo*y(j)*cosi+apo*yy(j)*sini
   26 continue
      go to 6
   12 k=1
      kk=2
      if(i1.le.0)i1=1
      if(i2.gt.npoint)i2=npoint
      phi=(i1-2)*dphi
      cosi=cos(phi)
      sini=sin(phi)
      do 5 i=i1,i2
      cosr=cosi
      cosi=cosi*dc-sini*ds
      sini=sini*dc+cosr*ds
      at=i-is(k)
      att=i-is(kk)
      at=at*at/(il(k)*il(k))
      att=att*att/(il(kk)*il(kk))
      if(att-at) 15,16,16
   15 k=kk
      kk=kk+1
      at=att
   16 at=at*at
      apo=1.-at
      if (apo) 5,5,11
   11 continue
      yf(i)=yf(i)+apo*y(j)*cosi+apo*yy(j)*sini
    5 continue
    6 continue
      do 7 i=1,ndim
   7  yf(i)=2*yf(i)
  77  format (f10.3)
  200 format (1x,15f6.0)
      do 10 i=1,np
   10 y(i)=yref(i)-yf(i)
c     desapodisation des donnees du pt nst au pt npoint
c     npoint nb de points de sortie le premier pont correspond
c     a l heure origine
      m=npoint
      write (*,*) m
      mi=(npoint-nst+1)/2
      dmi=mi
      p=1.e-07
      write (*,1000) m,nst,npoint
      do 300 i=nst,npoint
      di=i-nst
      ap=(1.-((di-dmi)**4)/(dmi**4))
      if(ap.gt.p)go to 301
      ap=1.
      write (*,1000) i
  301 yref(i)=yref(i)/ap
      yf(i)=yf(i)/ap
  300 y(i)=y(i)/ap
 1000 format(12i6)
c ecriture d un fichier ah
	modout=1
        call xdrf_ahopen(luout,outfile,modout,ierr)
      if(ierr.ne.0)then
		write(*,*)'erreur ouverture sortie ah'
		goto 99
	endif
        call xdrf_wrahrecord0(luout,ida,ih0,im0,s0,m,ierr,yref)
      	if(ierr.le.0)then
 	write(*,*)'erreur ecriture yref',ierr
 	goto 99
	endif
        call xdrf_wrahrecord0(luout,ida,ih0,im0,s0,m,ierr,yf)
      	if(ierr.le.0)then
 	write(*,*)'erreur ecriture yf',ierr
 	goto 99
	endif
        call xdrf_wrahrecord0(luout,ida,ih0,im0,s0,m,ierr,y)
      	if(ierr.le.0)then
 	write(*,*)'erreur ecriture y',ierr
 	goto 99
	endif
	call xdrf_ahclose(luout,ierr)
	if(ierr.ne.0) then
		write(*,*)'erreur fermeture fichier de sortie',ierr
		goto 99
	endif

      nst1=nst+20
c     write (*,101) (yref(i),i=nst,nst1)
 101  format (10f6.0)
	print*,'toto'
  109 continue
 99   continue
	goto 1
 7777 continue
      end
c------------------------------------------------------------------
      subroutine szak (r1,r2,ba,dtmin)
      pi=3.14159265389
      cdb=(6.*alog(10.))**(1./2.)/(2.*pi)
      a1=r1/(4.*pi)
       a2=r2/(24.*pi*pi)
      a2=abs(a2)
      a1=abs(a1)
      ai=(a2**(2./3.))*(3./2.)
      da=ai/100.
      dt=a1*2+(27./16.)*(a2**2)/(a1**2)*2
      alpha=a1-da
      do 1 i=1,102
      alpha=alpha+da
      dtt=(alpha+(a1**2)/alpha+(27./16.)*(a2**2)/(alpha**2))
      if(dtt.ge.dt) go to 2
    1 dt=dtt
      write (1,5) dtmin,r1,r2
      write (1,3)
    3 format (' recherche du filtre optimum impossible')
    2 alpha=alpha-da
      alpha=abs(alpha)
      ba=((alpha)**(-1./2.))*cdb
      dtmin=dt**(1./2.)
    5 format (e15.7)
      return
      end
c------------------------------------------------
      subroutine szper (tdeb,tfin,ipi,ips,n)
      common /sleep/ u(3,54),t(54),x(54)
      common /coond/ frq1,frq2,nf,dflat
      common /codisp/ nq1,nq2
      common /cogen/ in,iout,iperf,dnu,pas
      write (*,*) (t(i),i=1,n)
      nq1=1./(t(n)*dnu)+2.
      nq2=1./(t(2)*dnu)+1.
      nf=nq2-nq1+1
      frq1=(nq1-1)*dnu
      frq2=(nq2-1)*dnu
      if(ipi.eq.0) go to 30
      if(ips.eq.0) go to 30
      write (1,10) ipi,ips
   10 format (' periode debut et fin du plateau',2i10)
      dflat=1./ipi-1./ips
      go to 31
   30 dflat=0.
   31 return
      end
c-------------------------------------------------------
      subroutine nlogn(n,xr,xi,sign)
c      sign=-1.0  exp (-2*i*pi*...)
c      sign =1.0 (1/q) * exp(2*i*pi*...)
c      nmax=plus grande valeur de n avec
c      dimension x(2**n) ,m(nmax)
c
      dimension xr(2),xi(2),m(20)
      lx=2**n
      do 1 i=1,n
    1 m(i)=2**(n-i)
      do 4 l=1,n
      nbloc=2**(l-1)
      lbloc=lx/nbloc
      lbhaf=lbloc/2
      k=0
      do 4 ibloc=1,nbloc
      fk=k
      flx=lx
      v=sign*6.2831853*fk/flx
      wkr=cos(v)
      wki=sin(v)
      istat=lbloc*(ibloc-1)
      do 2 i=1,lbhaf
      j=istat+i
      jh=j+lbhaf
      qr=xr(jh)*wkr-xi(jh)*wki
      qi=xi(jh)*wkr+xr(jh)*wki
      xr(jh)=xr(j)-qr
      xi(jh)=xi(j)-qi
      xr(j)=xr(j)+qr
      xi(j)=xi(j)+qi
    2 continue
      do 3 i=2,n
      ii=i
      if (k-m(i)) 4,3,3
    3 k=k-m(i)
    4 k=k+m(ii)
      k=0
      do 8 j=1,lx
      if (k-j) 5,6,6
    6 holdr=xr(j)
      holdi=xi(j)
      xr(j)=xr(k+1)
      xi(j)=xi(k+1)
      xr(k+1)=holdr
      xi(k+1)=holdi
    5 do 7 i=1,n
      ii=i
      if (k-m(i)) 8,7,7
    7 k=k-m(i)
    8 k=k+m(ii)
      if (sign) 11,9,9
    9 do 10 i=1,lx
      xr(i)=xr(i)/flx
   10 xi(i)=xi(i)/flx
   11 return
      end
c-------------------------------------------------------------------
      subroutine sztr(y,nst,npoint)
      dimension y(1)
      som=0.
      do 1 i=nst,npoint
    1 som=y(i)+som
      som=som/(npoint-nst+1)
      do 2 i=nst,npoint
    2 y(i)=y(i)-som
      nst1=nst+20
      write (1,101) (y(i),i=nst,nst1)
  101 format(10f6.0)
      mi=(npoint-nst+1)/2
      do 3 i=nst,npoint
      di=i-nst
      dmi=mi
    3 y(i)=y(i)*(1.-((di-dmi)**4)/(dmi**4))
      return
      end
c--------------------------------------------------------------------
      subroutine great(alat1,alon1,alat2,alon2,dist,disd,gc,
     *  az12,az21)
c calculation of distance in km and in deg, length of great circle path,
c and azimuth
c  input (lat1,lon1) (lat2,lon2)
c  output	odist : distance in km
c		odisd : distance in degree
c		ogc   : length of great circle path
c		oz12  : azimuth of 2 at 1 measured clockwise from north
c		oz21  : azimuth of 1 at 2
       ath=6378.140
      bth=6356.755
      pi = 3.14159265
      rad = pi/180.
      h = 1. - bth*bth/(ath*ath)
      p = h/(1. - h)
      gr = alon1*rad
      tr = alat1*rad
      sintr =sin(tr)
      costr =cos(tr)
      if (sintr .eq. 0.) sintr = .00000100
      if (costr .eq. 0.) costr = .00000100
      r1 = ath/sqrt(1. - h*sintr*sintr)
      z1 = r1*(1. - h)*sintr
      g = alon2*rad
      t = alat2*rad
      if (t .eq. 0.) t = .0000100
      sint =sin(t)
      cost =cos(t)
      r2 = ath/sqrt(1. - h*sint*sint)
      dg = g - gr
      cosdg =cos(dg)
      sindg =sin(dg)
      dgr = gr - g
      dt = t - tr
      q = sint*costr/((1. + p)*cost*sintr) + h*r1*costr/(r2*cost)
      x = r2*cost*cosdg
      y = r2*cost*sindg
      z = r2*(1. - h)*sint
      az12 =atan2(sindg,(q - cosdg)*sintr)
      q = sintr*cost/(costr*sint*(1. + p)) + h*r2*cost/(r1*costr)
      az21 = atan2(sin(dgr),sint*(q-cos(dgr)))
      cos12 =cos(az12)
      cta2 = costr*costr*cos12*cos12
      p0 = p*(cta2 + sintr*sintr)
      b0 = (r1/(1. + p0))*sqrt(1. + p*cta2)
      e0 = p0/(1. + p0)
      gc = 2.*pi*b0*sqrt(1. + p0)*(1. - e0*(.25 + e0*(3./64.
     *                                          + 5.*e0/256.)))
      c0 = 1. + p0*(.25 - p0*(3./64. - 5.*p0/256.))
      c2 = p0*(-.125 + p0*(1./32. - 15.*p0/1024.))
      c4 = (-1./256. + 3.*p0/1024.)*p0*p0
      u0 =atan2(sintr,costr*cos12*sqrt(1. + p0))
      u =atan2(r1*sintr + (1. + p0)*(z - z1),(x*cos12 - y*sintr*
     *                                       sin(az12))*sqrt(1. + p0))
      disd = u - u0
      if (u .lt. u0) disd = pi + pi + disd
      dist = b0*(c0*( disd ) +c2*(sin(u + u) -sin(u0 + u0))
     *                       +c4*(sin(4.*u) -sin(4.*u0)))
      disd = disd/rad
      az12 = az12/rad
      az21 = az21/rad
      if (az12 .lt. 0.) az12 = 360. + az12
      if (az21 .lt. 0.) az21 = 360. + az21
      return
      end
c_____________________________________________________________________
      subroutine julien1(ij,im,ia,iy,id,isi)
c
c     converti ij,im,ia -> iy,id       pour isi=+1
c     converti iy,id    -> ij,im,ia    pour isi=-1
c
c     ex:  ij,im,ia = 21 02 80    iy,id = 80 52
c
c     integer*2 ij,im,ia,iy,id,isi
      dimension imo(12)
      data imo/31,28,31,30,31,30,31,31,30,31,30,31/
      if(isi.eq.-1)go to 10
      if(ia.eq.0)go to 20
      imo(2)=28
      if(mod(ia,4).eq.0)imo(2)=29
      id=0
      do 1 i=1,12
      if(i.eq.im)go to 2
1     id=id+imo(i)
2     id=id+ij
      iy=ia
      return
10    if(iy.eq.0)go to 20
      imo(2)=28
      if(mod(iy,4).eq.0)imo(2)=29
      idl=id
      do 3 i=1,12
      if(idl.le.imo(i))go to 4
3     idl=idl-imo(i)
4     ij=idl
      im=i
      ia=iy+1900
      return
20    ij=0
      im=0
      ia=0
      iy=0
      id=0
      return
      end
c-----------------------------------------------------------
      subroutine  lisse(n,x,y,dy,s,a,b,c,d,r)
c     cette subroutine permet le calcul de la fonction spline d'ajustement
c        d'ordre 2 sur les points x(i),y(i),i=1,n .
c        les x(i) sont donnes dans l'ordre croissant  . f est choisie telle
c        que sigma(((f(x(i))-y(i))/dy(i))**2) <s .
c     valeurs de sorties : a,b,c,d (vecteurs de dimension n+1)
c        f(t)=a(i)+b(i)*h+c(i)*h**2+d(i)*h**3  avec h=t-x(i-1) , x(i-1)<t<x(i)
c        f(t)=a(1)+b(1)*(t-x(1)) si t<x(1)
c        f(t)=a(n+1)+b(n+1)*(t-x(n)) si t>x(n)
c     r : vecteur de travail en double precision de dimension 8*(n+1)
c
      dimension x(1),y(1),dy(1),a(1),b(1),c(1),d(1),r(1)
      double precision bi,ci,di,ci1,di1,di2,r,dh,dg,ka,kb,
     1dyj,dp,ai,dv
      if (s.lt.1.e-6) go to 999
c
c     initialisations
c
      ss=s+1.e-06*s
      ka=2.d0/3.d0
      kb=1.d0/3.d0
      p=0.
      nit=0
      na=n + 1
      nb=na + na
      nc=nb+na
      nd=nc+na
      ne=nd+na
      nf=ne+na
      ng=nf+na
      nh=ng+na
      do 10 j=1,nh
   10 r(j)=0.d0
      do 15 i=1,na
      a(i)=0.
      b(i)=0.
      c(i)=0.
   15 d(i)=0.
c
c     calcul de c=q*ddq et de q*d=y
c
      h=x(2) -x(1)
      print 300,h
  300 format(5f12.7)
      dh=dble(h)
      f=(y(2)-y(1))/h
      do 20 i=3,n
      j=i-1
      g=h
      h=x(i)-x(j)
      dg=dh
      dh=dble(h)
      e=f
      f=(y(i)-y(j))/h
      a(i)=f-e
      r(nc+i)=ka*(dg+dh)
      r(nd+i)=kb*dh
      dyj=dble(dy(j))
      r(i)=-dyj*(1.d0/dg+1.d0/dh)
      r(na+i)=dble(dy(j-1))/dg
   20 r(nb+i)=dble(dy(i))/dh
      ia=na +2
      ib=nb+2
      do 30 i=3,n
      ia=ia +1
      ib=ib +1
      r(ne+i)=r(ia)*r(ia) +r(i)*r(i) +r(ib)*r(ib)
      r(nf+i)=r(i)*r(ia+1)+r(i+1)*r(ib)
   30 r(ib)=r(ib)*r(ia+2)
   35 if(nit.gt.200) go to 999
      nit=nit +1
      dp=dble(p)
c
c     decomposition choleski rr*=c
c
      do 40 i=3,n
      i1=i-1
      i2=i-2
      ai=dble(a(i))
      bi=r(ne+i)+r(nc+i)*dp
      ci=r(nf+i)+r(nd+i)*dp
      di=r(nb+i)
      tol=1.e-16*abs(sngl(bi))
      di1=dble(d(i1))
      ci1=dble(c(i1))
      di2=dble(d(i2))
      bi=bi-di2*di2-ci1*ci1
      if(sngl(bi).lt.tol) go to 999
      bi=1.d0/dsqrt(bi)
      r(i)=bi
      c(i)=sngl(bi*(ci-ci1*di1))
      d(i)=sngl(bi*di)
   40 r(ng+i)=(ai-ci1*r(ng+i1)-di2*r(ng+i2))*bi
c
c     resolution cu=y
c
      r(nh)=0.d0
      ii=nh-1
      ij=n
      r(ii)=r(ii)*r(ij)
      do 50 i=4,n
      ii=ii-1
      ij=ij-1
      ci=dble(c(ij))
      di=dble(d(ij))
   50 r(ii)=(r(ii)-ci*r(ii+1)-di*r(ii+2))*r(ij)
c
c     calcul de v=dqu et e=v*v
c
      res=0.
      h=0.
      f=0.
      ig=ng+1
      do 60 i=2,n
      ig=ig+1
      i1=i-1
      g=h
      h=x(i)-x(i1)
      e=f
      f=(sngl(r(ig+1)-r(ig)))/h
      b(i)=(f-e)*dy(i1)*dy(i1)
   60 res=res+b(i)*(f-e)
      b(na)=-f*dy(n)*dy(n)
      res=res-b(na)*f
c
c     test res>s
c
      if(res.lt.ss) go to 80
c
c     calcul de g=w*w et f=u*tu
c
      g=0.
      f=0.
      ia =na+2
      ic=nc+2
      id=nd+2
      ig=ng+2
      do 70 i=3,n
      ia=ia+1
      ic=ic+1
      id=id+1
      ig=ig+1
      dv=r(id-1)*r(ig-1)+r(ic)*r(ig)+r(id)*r(ig+1)
      ci1=dble(c(i-1))
      di2=dble(d(i-2))
      r(ia)=(dv-ci1*r(ia-1)-di2*r(ia-2))*r(i)
      g=g+sngl(r(ia)*r(ia))
   70 f=f+sngl(r(ig)*dv)
c
c     nouvelle valeur de p
c
      p=p+(res-sqrt(s*res))/(f-p*g)
      go to 35
c
c     calcul de a,b,c,d
c
   80 do 90 i=2,na
      c(i)=p*sngl(r(ng+i))
   90 a(i)=y(i-1)-b(i)
      do 100 i=2,n
      h=x(i)-x(i-1)
      d(i)=(c(i+1)-c(i))/(3.*h)
  100 b(i)=(a(i+1)-a(i))/h-h*(c(i)+h*d(i))
      b(1)=b(2)
      a(1)=a(2)
      b(na)=b(n)+(2.*c(n)+3.*d(n)*h)*h
      return
  999 do 2000 i=1,na
      a(i)=0.
      b(i)=0.
      c(i)=0.
 2000 d(i)=0.
      return
      end
c-----------------------------------------------------------
      function smoo(kt,om,a,b,c,d,omm)
      dimension om(1)
      dimension a(1),b(1),c(1),d(1)
c   test sur l'intervalle
      ktt=kt-1
      do 10 k=1,ktt
      test=om(k)
      if (omm.ge.test) go to 10
      if (k.eq.1) go to 8
      var=om(k-1)
      var=omm-var
      smoo=a(k)+b(k)*var+c(k)*var**2+d(k)*var**3
      go to 20
    8 smoo=a(1)+b(1)*(omm-test)
      go to 20
   10 continue
      smoo=a(kt)+b(kt)*(omm-test)
   20 continue
      return
      end
c-------------------------------------------------------------------------------
      subroutine szond (i,dist,nt,ntrain,i1,i2,irec,nk)
c     sousprogramme utilisant les derivees du lissage des vitesses de groupe
      dimension aa(55),bb(55),cc(55),ee(55)
      dimension d(3),dd(3),ud2(3)
      common /sleep/ u(3,54),t(54),x(54)
      common/tire/l1(30),l2(30),is(30),il(30)
      common/coond/ frq1,frq2,nf,dflat
      common/cogen/ in,iout,iperf,dnu,pas
      common/ond/a(3,55),b(3,55),c(3,55),e(3,55)
      er=3.5
      ef=1.3
      iflat=dflat
      flat=iflat
      if(iflat.eq.0) go to 5
      flat=ef/dflat
    5 continue
      frq=(i-1)*dnu
      n=nt
   24 dt=frq-x(n)
      if(dt) 21,22,22
   21 n=n-1
      go to 24
   22 continue
      do 3 k=1,nk
      do 33 m=1,nt
      aa(m)=a(k,m)
      bb(m)=b(k,m)
      cc(m)=c(k,m)
   33 ee(m)=e(k,m)
      ud2(k)= smoo(nt,x,aa,bb,cc,ee,frq)
      nn = n+1
      d(k) = b(k,nn) + 2.*c(k,nn)*dt + 3.*e(k,nn)*dt**2
      dd(k)= 2.*c(k,nn) + 6.*e(k,nn)*dt
      per = 1./frq
      ul = 1./ud2(k)
      df = x(nn)-frq
      if((dt.gt.dnu).and.(df.gt.dnu)) go to 3
c     write (1,100) frq,per,ud2(k),d(k),dd(k),ul
  100 format (1h ,6e15.6)
    3 continue
      do 1 j=1,ntrain
      nte=(j-1)/2.
      dt=nte*40030.
      npa=2*nte
      if(npa.ne.j-1) go to 2
      is(j)=(dist*ud2(1)+dt*ud2(3))/pas+1.5
      r1=dist*d(1)+dt*d(3)
      r2=dist*dd(1)+dt*dd(3)
      go to 10
    2 di=40030.-dist
      is(j)=(di*ud2(2)+dt*ud2(3))/pas+1.5
      r1=di*d(2)+dt*d(3)
      r2=di*dd(2)+dt*dd(3)
   10 call szak(r1,r2,ba,dtmin)
      all=er*dtmin
      if (all.ge.flat) go to 1
      all=flat
    1 il(j)=all/pas+0.5
      i1=is(1)-il(1)
      i2=is(ntrain)+il(ntrain)
      ntp=ntrain+1
      is(ntp)=is(ntrain)
      il(ntp)=il(ntrain)
      irec =0
      ntm=ntrain-1
      do 6 j=1,ntm
      jj=j+1
      ifin=is(j)+il(j)
      ideb=is(jj)-il(jj)
      if(ifin.le.ideb) go to 6
      irec=1
    6 continue
      if (irec.ne.0) go to 7
      do 11 j=1,ntrain
      l1(j)=is(j)-il(j)
   11 l2(j)=is(j)+il(j)
    7 continue
      return
      end

