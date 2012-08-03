c Version faite pour tourner avec la bibliotheque graphique Vogle

c  Version synt5:	tentative pour ameliorer le lissage de la
c			source. Dans les versions anterieures, on
c			lisse l'amplitude et la phase. Dans certains
c			cas, il y a un saut de 2*pi dans la phase qui
c			amene a des resultats aberrants.

c     parameter(ndim=4096, ndim2=ndim/2+1)
      parameter(ndim=8192, ndim2=ndim/2+1)
      dimension xr(ndim),xi(ndim),xrs(ndim),ainstr(ndim2),pinstr(ndim2),
     1ak(ndim2),aq(ndim2)
C
C Dimensions des tableaux:
C xr,xi,xrs .......................... 2**ip
C ainstr,pinstr,ak ................... 2**(ip-1)+1
C tb1,tb2,tb3,tb5 .................... nt
C c0,c1,c2,c3 ........................ nt+1
C r .................................. 8*(nt+1)
C
      character fi*80
      character bcd(4)*80, ligne*100
      real*8 pi, pi2
      character comp*12
      dimension amo0(10),amo(10)

      call getarg(1,comp)
c    dimension limitant dans ce pgm le nombre de seismes
c    inversibles.
      dimnbs=10
	print *,'composante: ', comp
      pi=3.1415926536
      pi2=2*pi
      r0=6371.
      io=10
      in=11
      is=12
      iq=13
      idessin=0
      itf=0
      write(*,'("voulez-vous dessiner ? ",$)')
      if(nonoui().eq.1)then
      xg=0.
      xd=60.
      yt=100.
      dh=5.
      y0=dh
      ey=0.      
      idessin=1
      write(*,'("veuillez PATIENTER...")')
      call initdes(xg,xd,y0,yt,7,0,0,0)
      write(*,'("veuillez patienter...")')
	call CreateRetainSeg(1)
      endif
c*
c* entree des parametres temporels et frequentiels du sismogramme
c*
21    write(*,'("entrer la distance epicentrale en km: ",$)')
      read(*,*) dist
      write(*,'("entrer les vitesses de groupe extremes en km/s: ",$)')
      read(*,*) v1,v2
      if(v1.gt.v2) go to 1
      t1=dist/v2
      t2=dist/v1
      go to 2
1     t1=dist/v1
      t2=dist/v2
2     dlt=t2-t1
      write(*,'("la duree du signal est ",f10.4," sec")')dlt
      write(*,'("voulez-vous imposer une autre duree ? ",$)')
      if(nonoui().eq.1)then
      write(*,'("entrer la duree en sec: ",$)')
      read(*,*)dlt
      endif
      dnu=1./dlt
      write(*,'("entrer la periode minimum en sec.: ",$)')
      read(*,*)f2
4     f2=1./f2
      nn=2*f2/dnu+1.5
      call szpu2(nn,np,ip)
      if(np.gt.ndim)then
      write(*,'("nombre de points trop grand (periode mini trop petite o
     *u duree trop grande)")')
      stop
      endif
C
C* Calcul de la frequence de Nyquist
      f2=np*dnu/2
C
      dt=dlt/np
      write(*,'("le pas en temps naturel du sismogramme est",/
     1 f10.4,"sec.")') dt
      write(*,'("voulez-vous imposer un pas different a la sortie de la 
     1T.F. ? ",$)') 
      if(nonoui().eq.0) go to 5
7     write(*,'("entrer le nouveau pas en seconde: ",$)')
      read(*,*) dtt
      if(dtt.lt.dt)go to 6
      write(*,'("il ne sert a rien d imposer un pas plus lache")')
      go to 7
6     n=dlt/dtt+ .5
      call szpu2(n,np,ip)
      if(np.gt.ndim)then
      write(*,'("nombre de points trop grand")')
      stop
      endif
      dlt=np*dtt
      dnu=1./dlt
      t2=t1+dlt
      write(*,'("la nouvelle duree du signal est ",f6.1," sec")')dlt
C
C* Calcul de la frequence maximum
C
      f2=np*dnu/2
      dt=dtt
5     continue
      npp=np/2+1

c*
c* entree de la T.F. du sismogramme entre la frequence zero et f2
c             convention de signe: la T.F. inverse est en
c                                  somme de exp(+iwt)dnu
c*
      call szopen(in,'("entrer le nom du fichier contenant la dispersion
     1: ",$)')
      read(in,'(a)')bcd(1)
      write(*,'(a)')bcd(1)
      read(in,*)km
      call szopen(is,'("entrer le nom du fichier contenant la source: ",
     1$)')
      read(is,'(a)')bcd(2)
      write(*,'(a)')bcd(2)
      read(is,*)kmm,dist0
      if(kmm.lt.km)km=kmm
      egeom=sqrt(sin(dist0/r0)/sin(dist/r0))
      print*,egeom
      call szopen(iq,'("entrer le nom du fichier contenant l attenuation
     1: ",$)')
      read(iq,'(a)')bcd(3)
      write(*,'(a)')bcd(3)
c--------------------------------------------------------------------
c  modification de la lecture des moments sismiques compatible avec la
c  lecture de plusieurs trajets(derniere modif,11/12/94).
c ancienne version (plus complexe) dans synt5.lr.F.old.
      write(*,*)'entrer le nombre de seismes inverses:'
      read (*,*) nbs
      if (nbs.gt.dimnbs) stop
     *'erreur de dimension dans SYNT3'
      write(*,*)'entrer le numero du seisme sur lequel vous
     *etes en train de travailler:'
      read (*,*) nos
c     lecture du nb de modes, des moments sismiques initial et inverse
      read(iq,*) kmm,(amo0(isei),amo(isei),isei=1,nbs)
      do 16 isei=1,nbs
      if (isei.eq.nos)then
      am0=amo0(isei)
      am=amo(isei)
      write(*,*)'no seisme',isei,'moment initial',am0
      write(*,*)'no seisme',isei,'moment inverse',am
      endif
16    continue
c--------------------------------------------------------------------
      if(kmm.lt.km)km=kmm
      write(*,'("les synthetiques sont calcules pour ",i2,"modes, voulez
     1-vous les conserver tous ? ",$)')km
      if(nonoui().ne.1)then
      write(*,'("entrer le nombre de modes a conserver: ",$)')
      read(*,*)km
      endif
      if(idessin.eq.1)then
      h=(yt-y0)/(km+1)-dh
      endif
      call szinst(npp,dnu,f2,ainstr,pinstr,bcd(4))
****   	write(22,*) 'Fonction de transfert extraite de synt2'//bcd(4)
****   	write(22,*) 1
****   	write(22,*) npp-1,dnu
****   	do 222 i=2,npp
****   	write(22,*)(i-1)*dnu,ainstr(i),pinstr(i)
****  222	continue
      do 10 i=1,np
10    xrs(i)=0.
      if(idessin.eq.1)yb=y0
      write(*,'("voulez vous stocker la T.F. ? ",$)')
      if(nonoui().eq.1)then
9     write(*,'("entrer le nom de fichier: ",$)')
      read(*,'(a)')fi
      open(io,status='new',form='unformatted',err=9,file=fi)
      write(io)npp,dnu,t1,dt,dlt,km,ip,dist
      itf=1
      endif
C
C* Boucle sur le rang k des modes
C
      write(*,'("voulez vous integrer la vitesse de groupe ? ",$)')
      ivg=nonoui()
      imark=0
      do 19 k=1,km
      if(ivg.eq.1)then
      call szvg(ak,dnu,npp,in,if11,if21)
      else
      call szvp(ak,dnu,npp,in,if11,if21)
      endif 
      call szs(xr,xi,dnu,is,if12,if22,npp,comp)
      call szq(aq,dnu,iq,if1,if2,npp)
      if1=max(if11,if12,if1)
      if2=min(if21,if22,if2)
      write(*,'("km,np,npp,ip,if1,if2",6i4)')km,np,npp,ip,if1,if2
C
C* Boucle sur les frequences
C

*      do 18 i=if1,if2

*c     Ecriture du spectre a la source incluant l'exp. geom.
*c     i.e. spectre final - instrument - phase propa - att. propa.
*c     pour comparaison avec Jeannot et Hugues
*
c     write(24,*) "Spectre a la source incluant l'exp. geom."
c     write(24,*) "1"
c     write(24,*) if2-if1+1
      do 18 i=if1,if2
      amsou=xr(i)*egeom*am/am0
      phsou=xi(i)
c     write(24,*) 1./((i-1)*dnu),phsou,phsou,amsou,amsou,phsou*180./pi

      f=(i-1)*dnu
      a=exp(-aq(i)*dist)*xr(i)*ainstr(i)*egeom*am/am0
      ph=-ak(i)*dist+xi(i)+pinstr(i)+pi2*f*t1
      xr(i)=a*cos(ph)/dt
18    xi(i)=a*sin(ph)/dt

      ido=if1-1
      do 11 i=1,ido
      xr(i)=0.
11    xi(i)=0.
      ido=if2+1
      do 12 i=ido,np
      xr(i)=0.
12    xi(i)=0.
C
C* Fin du calcul de la TF du mode k
C
      if(itf.eq.1) then
      write(io) (xr(i),xi(i),i=1,npp)
      endif
*      do 221 i=2,npp
*      perio=1./(i-1)/dnu
*      if(perio.ge.25..and.perio.le.300.) then
*        amp=sqrt(xr(i)*xr(i)+xi(i)*xi(i))
*        write(21,*) dist, amp, perio
*        endif
*221   continue
c
      write(*,'("voulez-vous calculer le signal temporel ? ",$)')
      if(nonoui().eq.0)go to 19
      imark=1
      call nlogn(ip,xr,xi,+1.)
      do 17 i=1,np
17    xr(i)=2*xr(i)
      do 20 i=1,np
20    xrs(i)=xrs(i)+xr(i)
      if(idessin.eq.1)then
      write(*,'("voulez vous dessiner le mode",i2," ? ",$)')k-1
      if(nonoui().eq.1)then
      yb=yb+h
      yh=yb+h
      call des(xi,xr,np,xg,xd,yb,yh,ey,dt)
      write(*,'("ey=",e10.3)')ey
      endif
      endif
19    continue
      close(in)
  191 read(is,'(a)') ligne
      if(ligne(:8).ne.'  moment') goto 191
      read(ligne,'(20x,e15.5)') amoment
      write(*,*) ' moments =',amoment, am0, am
      close(is)
      close(iq)
      if(itf.eq.1) then
      write(io) amoment
      write(*,*) ' amoment =',amoment
      close(io)
      endif
      if (imark.eq.0) goto 99
      write(*,'("voulez vous stocker le signal somme ? ",$)')
      if(nonoui().eq.1)then
      write(*,'("entrer le nom de fichier: ",$)')
      read(*,'(a)')fi
      open(io,file=fi)
*	bcd(4) = 1ere ligne du fichier-instrument
*	qui doit contenir: nom-station (a8) et nom-composante (a8)
	ligne=bcd(4)
	ligne='Synthetic       '//ligne(1:16)
      write(io,'(a)') ligne(1:lnblnk(ligne))
*	bcd(1) = 1ere ligne du fichier-dispersion
	ligne=bcd(1)
      write(io,'(a)') ligne(1:lnblnk(ligne))
	write(io,888) evla, evlo, evdp
	write(io,889) stla, stlo, stel
	write(*,*) "Entrer l'heure origine yyyy mm dd hh mm ss.ss"
	read(*,*) iye,ime,ide,ihe,imie,se
      write(io,'(2i3,f6.2,a)')ihe,imie,se,'               heure-origine'

	t1=t1+ihe*3600+imie*60+se
	ih1=int(t1/3600)
	s1=t1-3600*ih1
	im1=int(s1/60)
	s1=s1-60*im1
      write(io,'(2i3,f6.2,a,2i3.2,i5.4,a)') ih1,im1,s1,
     *				'   le',ide,ime,iye,'   (heure debut)'
      write(io,*)np,dt,dist
      write(io,'(10e14.6)') (xrs(i),i=1,np)
      close(io)
  888	format ("  event:  lat=",f11.5,"    lon=",f11.4,
     *		"    prof=",f11.5," km")
  889	format ("  statn:  lat=",f11.5,"    lon=",f11.4,
     *		"     alt=",f11.5," m")
      endif 
      if(idessin.eq.1)then
      write(*,'("voulez vous dessiner le signal somme ? ",$)')
      if(nonoui().eq.1)then
      yb=y0
      yh=h+y0
      ey=0.
      write(*,'("voulez-vous un menu de dessin special ? ",$)')
      if(nonoui().eq.1)then
      write(*,'("entrer yb et yh (entre 0 et 100): ",$)')
      read(*,*)yb,yh
      endif
      v1=dist/t1
      write(*,'(" dist,km  v1,km/s   t1,sec  dlt,sec"/4(f8.2,1x))')
     1dist,v1,t1,dlt
      call des(xi,xrs,np,xg,xd,yb,yh,ey,dt)
      write(*,'("ey=",e10.3)')ey
      endif
      endif
   99 write(*,'("voulez-vous calculer un autre sismogramme ? ",$)')
      if(nonoui().eq.1)then
      if(idessin.eq.1)then
      write(*,'("voulez-vous effacer l\'ecran ? ",$)')
      if(nonoui().eq.1)then
      call DelAllRetainSegs()
      call NewFrame()
      call CreateRetainSeg(1)
      endif
      ey=0.
      endif
      itf=0
      go to 21
      endif
      if(idessin.eq.1) then
	call CloseRetainSeg(1)
	call findes
	endif
      stop
      end
c---------------------------------------------------------------
      SUBROUTINE szinst(npp,dnu,f2,amp,phi,bcd)
      dimension AMP(1),phi(1)
      DIMENSION T(8),B(8),Tstd(8),Bstd(8)
      COMPLEX CS,CS2,CFAC,CDEN,CFTR,C1,CNUM,cfil
      DIMENSION CDEN(8),CFAC(8),OM(8) 
      character*80 bcd
      DATA Tstd/20.0,3635.,30.14,30.10,30.13,150.6,10.,10./
      DATA Bstd/0.707,0.711,0.966,0.2601,0.712,0.712,0.9238,0.3826/
c
C **************** Calcul de la fonction de transfert du sismo
C **************** Wielandt-Streickeisen;
c
c                  Ce sous-programme restitue une fonction de
c                  transfert par rapport au mouvement du sol
c
C **************** programme ecrit par D. Rouland, 1985, modifie
c                  Cara avril 1986, ****attention signe phase non
c                  encore verifie au 14 06 86.
c    signe de la phase HGLP verifie en aout 1986 (J.J. Leveque)
c    signe de la phase VBB-Echery verifie en Sept 1986 (J.J. Leveque)
c    signe de la phase WWSSN verifie le 4 mars 87 (J.J. Leveque)
c    signe de la phase SRO verifie le ?? ??? 88 (J.J. Leveque)
c    signe de la phase SRO reverifie le 26 avril 89 (J.J. Leveque)
c    signe de la phase Geoscope verifie le 26 avril 89 (J.J. Leveque)
c    compatible avec une TF(temporel --> frequentiel)
c                     en  NLOGN(X,Y,-1)
c    (la phase diminue quand la frequence augmente)
c
c	partie "Geoscope" modifiee le 5/11/91 pour tenir compte du
c     fait que les poles-zeros Geoscope sont maintenant exprimes
c     en omega et non plus en frequence.
  100 FORMAT(I3)
C
C AMPL est la SENSIBILITE MAX EN MV/MICRON/SEC**2 POUR LP.. 
C                             EN MV/MICRON/SEC POUR BRB 
C T(i) i=1,ki sont les periodes des poles
C B(i) i=1,ki sont les amortisements des poles
C INSTR = 1 POUR BRB,2 POUR POS,5 POUR HGLP,6 POUR VLP, 8 POUR VLP AVEC
C           FILTRE BUTTERWORTH, 10 POUR VBB ECHERY
C 
c      INSTR = -1 pour instrument type WWSSN (sismo+galva)
c      INSTR = -2 pour instrument type "SRO" defini par les
c		poles et les zeros de sa fonction de transfert
c      INSTR = -3 pour instrument "Geoscope" defini par les
c		poles et les zeros de sa fonction de transfert
c
      in=15
      call szopen(in,'("entrer le nom du fichier instrument: ",$)')
      pi=3.14159265
      PI2=2.*pi
      C1=(1.,0.)
      n1=2
      n2=f2/dnu+2.0001
      read(in,'(a)') BCD
      write(*,'(a)') BCD
      READ(IN,*)INSTR 
c     instr=-1 ==> instrument wwssn
      if(instr.eq.-1) then
		call instww(in,npp,dnu,f2,amp,phi,bcd)
c	        Impression de la fonction de transfert
c	        do 654 i=n2,n1,-1
c        write(*,*) log(1./(i-1)/dnu),amp(i),phi(i), 1./(i-1)/dnu
c 654 		continue
		close(in)
		return
		endif

	if(instr.eq.-2.or.instr.eq.-3) then
		call instsro(instr,in,npp,dnu,f2,amp,phi,bcd)
		close(in)
		return
		endif
c
c
c     les 3 lectures suivantes sont bidon si instr=10
c      et ne sont la que pour garder la meme structure de fichier.
      READ(IN,*) AMPL 
      ki=instr
      if(instr.eq.2) ki=1
      READ(IN,*)(T(I),I=1,KI) 
      READ(IN,*)(B(I),I=1,KI) 
      IF(INSTR.EQ.10) GO TO 2 
      do 321 i=1,ki
      if(t(i).ne.tstd(i).or.b(i).ne.bstd(i))then
        write(*,*)' Valeurs instrumentales non standard:'
        write(*,*)'  t=',t(i),' tstd=',tstd(i)
        write(*,*)'  b=',b(i),' bstd=',bstd(i)
      endif
  321 continue
      DO 1 I=1,KI 
    1 OM(I)=PI2/T(I)
      OMP2=OM(1)*OM(1)
C 
C  CALCUL DE LA FONCTION DE TRANSFERT SUR NFR  FREQUENCES
C
    2 DO 20 J=N1,N2 
      OMEG=PI2*DNU*(J-1)
      CS=CMPLX(0.,OMEG) 
      CS2=CS*CS 
      IF(INSTR.EQ.10) THEN
         CNUM=CS*CS2*(CS+3.5511) 
         CDEN(1)=CS2+2.8595*CS*0.01+7.948*0.0001 
         CDEN(2)=CS+3.431
         CDEN(3)=CS+2.1276*0.01    
         ampl=4390.
         cfil=(7.09/(7.09+cs))**4
         CFTR=cfil*CNUM*AMPL/(CDEN(1)*CDEN(2)*CDEN(3))
         GO TO 6 
         END IF
      CDEN(1)=CS2+2.*B(1)*OM(1)*CS+OM(1)*OM(1)*C1 
      IF(KI.EQ.1) THEN
         CNUM=CS2
         IF(INSTR.EQ.2) CNUM=CMPLX(OMP2,0.)
         CFTR=CNUM*AMPL/CDEN(1)
         GO TO 6 
         END IF
      CFAC(1)=OM(1)*OM(1)/CDEN(1) 
      DO 3 I=2,INSTR
    3 CDEN(I)=CS2+2.*B(I)*OM(I)*CS+OM(I)*OM(I)*C1 
      CFAC(2)=CS2/CDEN(2) 
      DO 4 I=3,INSTR
    4 CFAC(I)=OM(I)*OM(I)/CDEN(I) 
      CFTR=CMPLX(AMPL,0.) 
      DO 5 I=1,INSTR
    5 CFTR=CFTR*CFAC(I) 
    6 AMP(J)=CABS(CFTR)           
      phi(j)=atan2(aimag(cftr),real(cftr))
   20 CONTINUE
c  INSTR=1 ou INSTR=10 : Broadband ou Geotech ZLB
c  dans ce cas, la fonction de transfert est en vitesse
c  et non pas en acceleration.
      IF(INSTR.EQ.1.OR.INSTR.EQ.10)then  
         do 55 i=n1,n2
         phi(i)=phi(i)+pi/2
   55    amp(i)=amp(i)*pi2*(i-1)*dnu
      else
         do 66 i=n1,n2
         phi(i)=phi(i)+pi
         OMEG=PI2*DNU*(i-1)
   66    amp(i)=amp(i)*omeg*omeg
         endif
      phi(1)=0.
      amp(1)=0.
      nn=n2+1
      do 777 i=nn,npp
      phi(i)=0.
777   amp(i)=0.
c
      close(in)
      RETURN
      END
c---------------------------------------------------------------
      subroutine instww(in,npp,dnu,f2,amp,phi,bcd)
      dimension AMP(1),phi(1)
      character*80 bcd
C******************CALCUL DE LA FONCTION DE TRANSFERT D UN SISMOGRAPHE
C                  ELECTROMAGNETIQUE COUPLE A UN GALVANOMETRE 
c
c utilisation simple: prendre IK=1
c
c     DIMENSION AMP(150),PHI(150)
c     DIMENSION DP(150),PER(150)
c
      READ(in,*) IK 
      READ(in,*) CG,CS
      READ (in,*) TG,TS
      READ(in,*) AOG,AOS
      READ(in,*) RIS,RES,RIG,REG,S
c     write(*,*) IK 
c     write(*,*) CG,CS
c     write(*,*) TG,TS
c     write(*,*) AOG,AOS
c     write(*,*) RIS,RES,RIG,REG,S
C****************** 
C  CG  CS          CONSTANTES ELECTRIQUES DU GALVA ET DU SISMO
C  TG  TS          PERIODES PROPRES 
C  AOG AOS         AMORTISSEMENTS MECANIQUES
C  RIS,RES,RIG,RIS,S    RESISTANCES INTERNES,EXTERNES ET SHUNT
C  TDEB,TFIN,TPAS      PERIODES DU DEBUT ET DE LA FIN, PAS EN PERIODE 
C*****CORRECTION POUR COMPOSANTE NORD 
      SM=10.655 
      SLOG=0.345
C*************
      SM=7.2
      SLOG=0.34 
      SLR=0.42
      GMI=4.2E-07 
      D2=2. 
C                  SLOG DISTANCE AU CENTRE DE GRAVITE, TOUJOURS 
C                  INFERIEURE A SLR LONGUEUR REDUITE
C                  SM MASSE MOBILE
C                  GMI MOMENT D INERTIE 
C                  D2 LEVIER OPTIQUE
C   COEFFICIENT DE COUPLAGE STANDARD AU CARRE  PRIS LORSQUE IK=1
      SI2S=0.04 
C 
      PI=3.1415926
      PI2=PI*2
      S2=S*S
      RS=RES+RIS
      RG=REG+RIG
      Q2=(RS+S)*(RG+S)-S2 
      ALPHA1=CG*(RS+S)/Q2 
      BETA1=CS*(RG+S)/Q2
      ALPHA=AOG+ALPHA1
      BETA=AOS+BETA1
      IF(IK.NE.0)  ALPHA=1. 
      IF(IK.NE.0)  BETA=1.
      OMG=PI2/TG
      OMS=PI2/TS
      RO=SQRT(OMS/OMG)
      SI2=(ALPHA1*BETA1*S2)/(ALPHA*BETA*(Q2+S2))
      IF(IK.NE.0)  SI2=SI2S 
      A=2*(ALPHA/RO+BETA*RO)
      B=RO*RO+1./RO/RO+4*(1.-SI2)*ALPHA*BETA
      C=2*(ALPHA*RO+BETA/RO)
      V=SQRT(SLR*SM*4.*SI2*ALPHA*BETA/SLOG/GMI) 
c     N=(TFIN-TDEB)/TPAS+1.5
      ONOR=SQRT(OMS*OMG)
c
      n1=2
      n2=f2/dnu+2.0001
c
      DO 1 I=n1,n2
c     T=(I-1)*TPAS+TDEB 
c     OB=PI2/T/ONOR 
      omega=pi2*dnu*(i-1)
      ob=omega/onor
      OB2=OB*OB 
      OB3=OB2*OB
      AMP(I)=D2*V/SQRT((OB-B/OB+1./OB3)**2+(C/OB2-A)**2)
      X=OB3-OB*B+1./OB
      XX=OB2*A-C
      PHI(I)=-ATAN2(X,XX)
c     DP = derivee de la phase
c     R1=-(3*OB2-B-1./OB2)*XX 
c     R2=2*X*OB*A 
c     R3=XX*XX+X*X
c     DP(I)=(R1+R2)/R3/ONOR 
    1 continue
c
      phi(1)=0.
      amp(1)=0.
      nn=n2+1
      do 777 i=nn,npp
      phi(i)=0.
777   amp(i)=0.
      END 
c---------------------------------------------------------------
      subroutine instsro(instr,in,npp,dnu,f2,amp,phi,bcd)
      dimension amp(1),phi(1), zr(30),zi(30), polr(30),poli(30)
      character*80 bcd

c   Calcul de la fonction de transfert a partir des poles et des zeros

c	irep=1  ==> deplacement
c	irep=2  ==> vitesse
c	irep=3  ==> acceleration

      complex poles(30),zeros(30),ZZ,rep
	character*10 keyw

      pi2=2.*3.1415927

	if(instr.eq.-2) then

	 write(0,*)' Instrument SRO'
    1 	read(in,*)keyw,a
	write(0,'(3a)')'--',keyw,'--'
	if(keyw(:5).eq.'ZEROS') then
		nz=a+.1
		do 10 i=1,nz
		read(in,*) ar,ai
		zeros(i)=cmplx(ar,ai)
   10		continue
		iz=1
		endif
	if(keyw(:5).eq.'POLES') then
		np=a+.1
		do 11 i=1,np
		read(in,*) ar,ai
		poles(i)=cmplx(ar,ai)
   11		continue
		ip=1
		endif
	if(keyw(:8).eq.'CONSTANT') then
		a0=a
		ic=1
		endif
	if(iz.ne.1.or.ip.ne.1.or.ic.ne.1) goto 1

	endif

	if(instr.eq.-3) then
c	nouvelle convention Geoscope = convention GDSN/SAC
c         ( C * ( i.w - z ) / ( i.w - p )  )

c	On suppose ici que la fonction de transfert Geoscope
c	est donnee en deplacement (si ce n'est pas le cas,
c     c'est facile d'y remedier en ajoutant 1 ou 2 zeros).
c	Les amplitudes obtenues sont en digit/micron

	write(0,*)' Instrument Geoscope'
	read(in,*) a,nz,np
	read(in,*) (zr(i), i=1,nz)
	read(in,*) (zi(i), i=1,nz)
	read(in,*) (polr(i), i=1,np)
	read(in,*) (poli(i), i=1,np)
	k=0
	do 21 i=1,nz
	k=k+1
	zeros(k)=cmplx(zr(i),zi(i))
   21	continue
	nz=k

	k=0
	do 22 i=1,np
	k=k+1
	poles(k)=cmplx(polr(i),poli(i))
   22	continue
	np=k
	a0=a
c	passage digit/m --> digit/micron
	a0=a0*1.e-6

		write (*,*) '-2'
		write(*,'(3x,a,i6)') 'ZEROS  ',nz
		do 23 i=1,nz
		write(*,*) real(zeros(i)), aimag(zeros(i))
   23		continue
		write(*,'(3x,a,i6)') 'POLES  ',np
		do 24 i=1,np
		write(*,*) real(poles(i)), aimag(poles(i))
   24		continue
		write(*,*) 'CONSTANT  ',a0,'   DU/micron (pour HGLP)'

	endif

c	ici, calcul de la fonction de transfert en deplacement
	irep=1

	n1=2
	n2=f2/dnu+2.0001

	do 5 i=n1,n2
	omega=pi2*dnu*(i-1)
	zz=rep(a0,np,poles,nz,zeros,omega,irep)
	amp(i)=cabs(zz)
	phi(i)=atan2(aimag(zz),real(zz))
    5	continue

	amp(1)=0.
	phi(1)=0.
	do 6 i=n2+1,npp
	amp(i)=0.
	phi(i)=0.
    6	continue

c     signe de la phase: il doit etre compatible avec
c     les autres instruments (WWSSN, HGLP, ...) qui supposent
c     que la TF(temps-->frequence) est en -iwt  et donc que
c     la phase instrumentale decroit quand la frequence augmente.
c     Il semble qu'ici, il faille: ph = + atan2(... , ...)

      end
c---------------------------------------------------------------------
      complex function rep(a0,np,poles,nz,zeroes,omega,ifl)
C
C Resp returns T(s)=a0*(s-z(1))*...*(s-z(m))/(s-p(1))*...*(s-p(n)),
C where p,z,and s are complex.  s = cmplx(0.,omega).
C ifl flags the callers desire for
C displacemt (ifl=1),velocity (ifl=2),
C or acceleration (ifl=3).
C
C
C Calls no other routine.
C
C      Programmed by R. Buland
C         August 23,1979
C      Updated to double precision January 1986
C                    M. Zirbes
C
C ANGULAR FREQUENCY - INPUT
      REAL OMEGA
C
      COMPLEX*16ZZ, S, ZD
      DOUBLE PRECISION OO
C
C NORMALIZATION CONSTANT
      REAL A0
C NUMBER OF POLES
      INTEGER NP
C POLES
      COMPLEX POLES(30)
C NUMBER OF ZEROES
      INTEGER NZ
C ZEROES
      COMPLEX ZEROES(20)
C
      IF (.NOT.(np .LE. 0)) GOTO 2000
        rep = (0.0)
        RETURN
2000  CONTINUE
      OO = OMEGA*1D0
      ZD = DCMPLX(1.D5, 0.D0)
      S = DCMPLX(0.D0, OO)
C
C  COMPUTE THE CONTRIBUTION OF THE ZEROES IF ANY
C
      iz = nz
      DO 2020 J = 1, IZ
        ZZ = ZEROES(J)
        ZD = ZD*(S - ZZ)
C
C  COMPUTE THE CONTRIBUTION OF THE POLES
C
2020  CONTINUE
	ip=np
      DO 2040 J = 1, IP
        ZZ = POLES(J)
        ZD = ZD/(S - ZZ)
2040  CONTINUE
      ZD = ZD*1D - 5
      OO = A0*1.D0
C
C  DISPLACEMENT
C
      IF (.NOT.(IFL .EQ. 1)) GOTO 2060
        ZD = OO*ZD
C
C  VELOCITY
C
        GOTO 2070
2060  CONTINUE
      IF (.NOT.(IFL .EQ. 2)) GOTO 2080
        ZD = OO*ZD/S
C
C  ACCELERATION
C
        GOTO 2070
2080  CONTINUE
      IF (.NOT.(IFL .EQ. 3)) GOTO 2100
        ZD = OO*ZD/(S*S)
        GOTO 2070
2100  CONTINUE
        ZD = (0.D0, 0.D0)
2070  CONTINUE
      rep = ZD
      RETURN
      END
c---------------------------------------------------------------
      subroutine szs(as,ps,dnu,is,if1,if2,npp,comp)
	parameter(NPERMAX=30)
	real*8 pi
      DOUBLE PRECISION R
      dimension R(8*(NPERMAX+1)),c0(NPERMAX+1),c1(NPERMAX+1),
     *  c2(NPERMAX+1), c3(NPERMAX+1),ph(NPERMAX),u(NPERMAX),
     *  q(NPERMAX),f(NPERMAX),t(NPERMAX)
      dimension as(1),ps(1)
      character*12 comp
      ndim=NPERMAX
      pi=3.1415926536
      DC=1000.
      read(is,*)nt
      if(nt.gt.ndim)then
      write(*,'("nt=",i10," trop grand, valeur forcee nt=",i2)')nt,ndim
      nt=ndim
      endif
      if(comp.eq.'longitudinal'.or.comp.eq.'L'.or.comp.eq.'l')then
      		do 109 i=1,nt
      		READ (is,*) t(I),bidon,ph(I),bidon, u(i)
109   		continue
		else
      		do 110 i=1,nt
      		READ (is,*) t(I),ph(I),bidon, u(i),bidon
110   		continue
		endif
      tsup=t(1)
      tmin=t(1)
      DO 6 I=1,NT
      IF (t(I).GT.TSUP) TSUP=t(I)
6     IF (t(I).LT.TMIN) TMIN=t(I)
      if1=1./tsup/dnu+1.5
      if2=1./tmin/dnu+1.5
      if(if2.gt.npp)then
      write(*,'("szs: if2.gt.npp, correction faite en posant if2=npp")')
      if2=npp
      endif
      if(t(nt).gt.t(1))then
      call szrenv(t,nt)
      call szrenv(u,nt)
      call szrenv(ph,nt)
      endif

c	Lissage de la phase

c		elimination des sauts de 2*pi
	phref=ph(1)
	phdec=0.
      DO 4 I=1,NT
	test=ph(i)+phdec-phref
	if(abs(test).gt.pi) then
		phdec=phdec-2*pi*test/abs(test)
		endif
c           write(25,*) t(i),ph(i),ph(i)+phdec,phref
	ph(i)=ph(i)+phdec
	phref=ph(i)
c	lissage
      F(I)=1./T(I)
      Q(I)=ph(I)/dc
    4	continue
      S=NT
      SO=S/5
      ILI=0
   18 CONTINUE
      CALL LISSE (NT,F,ph,Q,S,C0,C1,C2,C3,R)
      IF (C0(1).NE.0.) GO TO 17
      write(*,213)
  213 FORMAT (' INCIDENT DANS LISSE, routine szs()')
      ILI=ILI+1
      S=SO*ILI
      IF (ILI.LT.10) GO TO 18
17    continue
      DO 10 I=IF1,IF2
      FRQ=(I-1)*DNU
      Vp=SMOO(NT,F,C0,C1,C2,C3,FRQ)
      ps(i)=vp
   10 continue

c	Lissage de l'amplitude
      DO 104 I=1,NT
      F(I)=1./T(I)
104   Q(I)=u(i)/dc
      S=NT
      SO=S/5
      ILI=0
  108 CONTINUE
      CALL LISSE (NT,F,U,Q,S,C0,C1,C2,C3,R)
      IF (C0(1).NE.0.) GO TO 107
      write(*,213)
      ILI=ILI+1
      S=SO*ILI
      IF (ILI.LT.10) GO TO 108
107   continue
      DO 100 I=IF1,IF2
      FRQ=(I-1)*DNU
      as(i)=SMOO(NT,F,C0,C1,C2,C3,FRQ)
  100 continue
      return
      end
c---------------------------------------------------------------
      subroutine szq(aq,dnu,iq,if1,if2,npp)
	parameter(NPERMAX=30)
      DOUBLE PRECISION R
      dimension R(8*(NPERMAX+1)),c0(NPERMAX+1),c1(NPERMAX+1),
     *  c2(NPERMAX+1), c3(NPERMAX+1),aqq(NPERMAX),q(NPERMAX),
     *  f(NPERMAX),t(NPERMAX)
      dimension aq(1)
      ndim=NPERMAX
      DC=1000.   
      read(iq,*)nt
      if(nt.gt.ndim)then
      write(*,'("nt=",i10," trop grand, valeur forcee nt=",i2)')nt,ndim
      nt=ndim
      endif
      do 109 i=1,nt
      READ (iq,*) t(I),aqq(I)
109   continue
      tsup=t(1)  
      tmin=t(1)
      DO 6 I=1,NT
      IF (t(I).GT.TSUP) TSUP=t(I)
6     IF (t(I).LT.TMIN) TMIN=t(I)
      if1=1./tsup/dnu+1.5
      if2=1./tmin/dnu+1.5
      if(if2.gt.npp)then
      write(*,'("szq: if2.gt.npp, correction faite en posant if2=npp")')
      if2=npp
      endif
      if(t(nt).gt.t(1))then
      call szrenv(t,nt)
      call szrenv(aqq,nt)
      endif
      DO 4 I=1,NT
      F(I)=1./T(I)
    4 Q(I)=aqq(i)/dc 
      S=NT
      SO=S/5
      ILI=0
   18 CONTINUE
      CALL LISSE (NT,F,aqq,Q,S,C0,C1,C2,C3,R)
      IF (C0(1).NE.0.) GO TO 17
      write(*,213)
  213 FORMAT (' INCIDENT DANS LISSE, routine szq()')
      ILI=ILI+1
      S=SO*ILI
      IF (ILI.LT.10) GO TO 18
17    continue
      DO 10 I=IF1,IF2
      FRQ=(I-1)*DNU
      aq(i)=SMOO(NT,F,C0,C1,C2,C3,FRQ)
   10 continue
      return
      end
c---------------------------------------------------------------
      subroutine szvg(ak,dnu,npp,in,if1,if2)
	parameter(NPERMAX=30)
      DOUBLE PRECISION R
      dimension R(8*(NPERMAX+1)),c0(NPERMAX+1),c1(NPERMAX+1),
     *  c2(NPERMAX+1), c3(NPERMAX+1),vgg(NPERMAX),u(NPERMAX),
     *  q(NPERMAX),f(NPERMAX),t(NPERMAX)
      dimension ak(1)
C
C        DC PRECISION SUR LA VITESSE DE PHASE EN KM/SEC
C
      ndim=NPERMAX
      DC=1.E-03
      PI=3.14159265
      read(in,*)nt
      if(nt.gt.ndim)then
      write(*,'("nt=",i10," trop grand, valeur forcee nt=",i2)')nt,ndim
      nt=ndim
      endif
      do 257 i=1,nt
257   READ (in,*) t(I),U(i),vgg(I)
      READ (in,*) VC
      vggMAX=vgg(1)
      vggMIN=vgg(1)
      TSUP=t(1)
      TMIN=t(1)
      if(t(nt).gt.t(1))then
      call szrenv(t,nt)
      call szrenv(vgg,nt)
      endif
      DO 6 I=1,NT
      IF (t(I).GT.TSUP) TSUP=t(I)
      IF (t(I).LT.TMIN) TMIN=t(I)
      IF (vgg(I).LT.vggMIN) vggMIN=vgg(I)
    6 IF (vgg(I).GT.vggMAX) vggMAX=vgg(I)
      if1=1./tsup/dnu+1.5
      if2=1./tmin/dnu+1.5
      if(if2.gt.npp)then
      write(*,'("szvg:if2.gt.npp, correction faite en posant if2=npp")')
      if2=npp
      endif
      UL1=1./(VC*T(1))
      DO 4 I=1,NT
      CT=vgg(I)
      F(I)=1./T(I)
      Q(I)=DC/(CT**2)
    4 U(I)=1./CT
      S=NT
      SO=S/5
      ILI=0
   18 CONTINUE
      CALL LISSE (NT,F,U,Q,S,C0,C1,C2,C3,R)
      IF (C0(1).NE.0.) GO TO 17
      write(*,213)
  213 FORMAT (' INCIDENT DANS LISSE, routine szvg()')
      ILI=ILI+1
      S=SO*ILI
      IF (ILI.LT.10) GO TO 18
      DO 19 I=1,NT
   19 U(I)=0.
      GO TO 20
   17 continue
      DO 5 I=1,NT
      FRQ=F(I)
      CALL SZSOML (C0,C1,C2,C3,FRQ,F,NT,UL)
      UL=UL+UL1
    5 U(I)=FRQ/UL
   20 CONTINUE
      DO 10 I=IF1,IF2
      FRQ=(I-1)*DNU
      CALL SZSOML (C0,C1,C2,C3,FRQ,F,NT,UL)
      CC=UL+UL1
      ak(i)=2*PI*CC
   10 continue
      return
      end
c---------------------------------------------------------------
      subroutine szvp(ak,dnu,npp,in,if1,if2)
	parameter(NPERMAX=30)
      DOUBLE PRECISION R
      dimension R(8*(NPERMAX+1)),c0(NPERMAX+1),c1(NPERMAX+1),
     *  c2(NPERMAX+1), c3(NPERMAX+1),vph(NPERMAX),u(NPERMAX),
     *  q(NPERMAX),f(NPERMAX),t(NPERMAX)
      dimension ak(1)
C
C        DC PRECISION SUR LA VITESSE DE PHASE EN KM/SEC
C
      ndim=NPERMAX
      DC=1.E-03
      PI=3.14159265
      read(in,*)nt
      if(nt.gt.ndim)then
      write(*,'("nt=",i10," trop grand, valeur forcee nt=",i2)')nt,ndim
      nt=ndim
      endif
      do 257 i=1,nt
257   READ (in,*)t(I),vph(I),U(I)
      tsup=t(1)
      tmin=t(1)
      DO 6 I=1,NT
      IF (t(I).GT.TSUP) TSUP=t(I)
6     IF (t(I).LT.TMIN) TMIN=t(I)
      if1=1./tsup/dnu+1.5
      if2=1./tmin/dnu+1.5
      if(if2.gt.npp)then
      write(*,'("szvp:if2.gt.npp, correction faite en posant if2=npp")')
      if2=npp
      endif
      if(t(nt).gt.t(1))then
      call szrenv(t,nt)
      call szrenv(vph,nt)
      endif
      DO 4 I=1,NT
      CT=vph(I)
      F(I)=1./T(I)
      Q(I)=DC/(CT**2)
    4 U(I)=CT
      S=NT
      SO=S/5
      ILI=0
   18 CONTINUE
      CALL LISSE (NT,F,U,Q,S,C0,C1,C2,C3,R)
      IF (C0(1).NE.0.) GO TO 17
      write(*,213)
  213 FORMAT (' INCIDENT DANS LISSE, routine szvp()')
      ILI=ILI+1
      S=SO*ILI
      IF (ILI.LT.10) GO TO 18
17    continue
      DO 10 I=IF1,IF2
      FRQ=(I-1)*DNU
      Vp=SMOO(NT,F,C0,C1,C2,C3,FRQ)
      ak(i)=2*PI*FRQ/Vp
   10 continue
      return
      end
      function nonoui()
      character*10 yo(14),yn(8)
      character c*10
      data yo/'oui','ouai','da','ya','yes','o','y','si','yo','ye','yeeh'
     1,'ok','positif','affirmatif'/
      data yn/'non','n','no','niet','nicht','neni','nein','pas'/
C
C Cette fonction retourne 1 si reponse positive
C                         0 si reponse negative
C                        -1 si reponse incorrecte apres 3 essais
C
      nonoui=-1
      k=0
3       read(*,'(a)')c
      do 1 i=1,14
      if(c.eq.yo(i))nonoui=1
1      continue
      do 2 i=1,8
      if(c.eq.yn(i))nonoui=0
2      continue
      if(nonoui.eq.-1)then
      if(k.ge.3)then
      write(*,'("advienne ce qu il adviendra")')
      return
      else
      k=k+1
      write(*,'("repondez a la question, rappel no ",i2," ",$)')k
      go to 3
      endif
      else
      return
      endif
      end
c---------------------------------------------------------------
      SUBROUTINE SZSOML (C0,C1,C2,C3,FRQ,FR,NT,SS)
C******************C0 C1 C2 C3    TABLEAU DES COEFFICIENTS DE LA
C                                 SPLINE DE LISSAGE *****DIMENSION  NT+1
C******************FR  TABLEAU DES ABCISSES DIMENSION  NT
C******************FRQ  ABCISSE OU L ON CALCULE    LA VALEUR DE L INTEGR
C                       FRQDOIT ETRE COMPRIS ENTRE  FR(1)ET FR(NT)
      DIMENSION C0(1),C1(1),C2(1),C3(1),FR(1)
      SS=0.
      S=FR(1)-FRQ
      IF (S.EQ.0.) RETURN
      DO 1 I=1,NT
      ST=FR(I)-FRQ
      IF (ST*S.LT.0.) GO TO 2
    1 CONTINUE
      write(*,3) FRQ
    3 FORMAT (' INCIDENT      FRQ',F10.5)
    2 II=I-1
      DF=FRQ-FR(II)
      S=C0(I)*DF+C1(I)*(DF**2)/2+C2(I)*(DF**3)/3+C3(I)*(DF**4)/4
      IF (II.EQ.1) GO TO 5
      DO 4 K=2,II
      KK=K-1
      DF=FR(K)-FR(KK)
    4 SS=SS+C0(K)*DF+C1(K)*(DF**2)/2+C2(K)*(DF**3)/3+C3(K)*(DF**4)/4
    5 SS=SS+S
      RETURN
      END
c---------------------------------------------------------------
      SUBROUTINE  LISSE(N,X,Y,DY,S,A,B,C,D,R)
C-----------------------------------------------------------------------
C     CETTE SUBROUTINE PERMET LE CALCUL DE LA FONCTION SPLINE D'AJUSTEME
C        D'ORDRE 2 SUR LES POINTS X(I),Y(I),I=1,N .
C        LES X(I) SONT DONNES DANS L'ORDRE CROISSANT  . F EST CHOISIE TE
C        LES X(I) SONT DONNES DANS L'ORDRE CROISSANT  . F EST CHOISIE TE
C        QUE SIGMA(((F(X(I))-Y(I))/DY(I))**2) )S .
C     VALEURS DE SORTIES : A,B,C,D (VECTEURS DE DIMENSION N+1)
C     VALEURS DE SORTIES : A,B,C,D (VECTEURS DE DIMENSION N+1)
C        F(T)=A(I)+B(I)*H+C(I)*H**2+D(I)*H**3  AVEC H=T-X(I-1) , X(I-1))
C        F(T)=A(1)+B(1)*(T-X(1)) SI T)X(1)
C        F(T)=A(N+1)+B(N+1)*(T-X(N)) SI T>X(N)
C     R : VECTEUR DE TRAVAIL EN DOUBLE PRECISION DE DIMENSION 8*(N+1)
C
C-----------------------------------------------------------------------
c Modif eric 6/02/2002 compilation sous linux
c     DIMENSION X(1),Y(1),DY(1),A(1),B(1),C(1),D(1),R(1)
C-----------------------------------------------------------------------
      DIMENSION X(*),Y(*),DY(*),A(*),B(*),C(*),D(*),R(*)
      DOUBLE PRECISION BI,CI,DI,CI1,DI1,DI2,R,DH,DG,KA,KB,
     1DYJ,DP,AI,DV
      IF (S.LT.1.E-6) GO TO 999
C
C     INITIALISATIONS
C
      SS=S+1.E-06*S
      KA=2.D0/3.D0
      KB=1.D0/3.D0
      P=0.
      NIT=0
      NA=N + 1
      NB=NA + NA
      NC=NB+NA
      ND=NC+NA
      NE=ND+NA
      NF=NE+NA
      NG=NF+NA
      NH=NG+NA
      DO 10 J=1,NH
   10 R(J)=0.D0
      DO 15 I=1,NA
      A(I)=0.
      B(I)=0.
      C(I)=0.
   15 D(I)=0.
C
C     CALCUL DE C=Q*DDQ ET DE Q*D=Y
C
      H=X(2) -X(1)
      DH=DBLE(H)
      F=(Y(2)-Y(1))/H
      DO 20 I=3,N
      J=I-1
      G=H
      H=X(I)-X(J)
      DG=DH
      DH=DBLE(H)
      E=F
      F=(Y(I)-Y(J))/H
      A(I)=F-E
      R(NC+I)=KA*(DG+DH)
      R(ND+I)=KB*DH
      DYJ=DBLE(DY(J))
      R(I)=-DYJ*(1.D0/DG+1.D0/DH)
      R(NA+I)=DBLE(DY(J-1))/DG
   20 R(NB+I)=DBLE(DY(I))/DH
      IA=NA +2
      IB=NB+2
      DO 30 I=3,N
      IA=IA +1
      IB=IB +1
      R(NE+I)=R(IA)*R(IA) +R(I)*R(I) +R(IB)*R(IB)
      R(NF+I)=R(I)*R(IA+1)+R(I+1)*R(IB)
   30 R(IB)=R(IB)*R(IA+2)
   35 IF(NIT.GT.200) GO TO 999
      NIT=NIT +1
      DP=DBLE(P)
C
C     DECOMPOSITION CHOLESKI RR*=C
C
      DO 40 I=3,N
      I1=I-1
      I2=I-2
      AI=DBLE(A(I))
      BI=R(NE+I)+R(NC+I)*DP
      CI=R(NF+I)+R(ND+I)*DP
      DI=R(NB+I)
      TOL=1.E-16*ABS(SNGL(BI))
      DI1=DBLE(D(I1))
      CI1=DBLE(C(I1))
      DI2=DBLE(D(I2))
      BI=BI-DI2*DI2-CI1*CI1
      IF(SNGL(BI).LT.TOL) GO TO 999
      BI=1.D0/DSQRT(BI)
      R(I)=BI
      C(I)=SNGL(BI*(CI-CI1*DI1))
      D(I)=SNGL(BI*DI)
   40 R(NG+I)=(AI-CI1*R(NG+I1)-DI2*R(NG+I2))*BI
C
C     RESOLUTION CU=Y
C
      R(NH)=0.D0
      II=NH-1
      IJ=N
      R(II)=R(II)*R(IJ)
      DO 50 I=4,N
      II=II-1
      IJ=IJ-1
      CI=DBLE(C(IJ))
      DI=DBLE(D(IJ))
   50 R(II)=(R(II)-CI*R(II+1)-DI*R(II+2))*R(IJ)
C
C     CALCUL DE V=DQU ET E=V*V
C
      RES=0.
      H=0.
      F=0.
      IG=NG+1
      DO 60 I=2,N
      IG=IG+1
      I1=I-1
      G=H
      H=X(I)-X(I1)
      E=F
      F=(SNGL(R(IG+1)-R(IG)))/H
      B(I)=(F-E)*DY(I1)*DY(I1)
   60 RES=RES+B(I)*(F-E)
      B(NA)=-F*DY(N)*DY(N)
      RES=RES-B(NA)*F
C
C     TEST RES>S
C
      IF(RES.LT.SS) GO TO 80
C
C     CALCUL DE G=W*W ET F=U*TU
C
      G=0.
      F=0.
      IA =NA+2
      IC=NC+2
      ID=ND+2
      IG=NG+2
      DO 70 I=3,N
      IA=IA+1
      IC=IC+1
      ID=ID+1
      IG=IG+1
      DV=R(ID-1)*R(IG-1)+R(IC)*R(IG)+R(ID)*R(IG+1)
      CI1=DBLE(C(I-1))
      DI2=DBLE(D(I-2))
      R(IA)=(DV-CI1*R(IA-1)-DI2*R(IA-2))*R(I)
      G=G+SNGL(R(IA)*R(IA))
   70 F=F+SNGL(R(IG)*DV)
C
C     NOUVELLE VALEUR DE P
C
      P=P+(RES-SQRT(S*RES))/(F-P*G)
      GO TO 35
C
C     CALCUL DE A,B,C,D
C
   80 DO 90 I=2,NA
      C(I)=P*SNGL(R(NG+I))
   90 A(I)=Y(I-1)-B(I)
      DO 100 I=2,N
      H=X(I)-X(I-1)
      D(I)=(C(I+1)-C(I))/(3.*H)
  100 B(I)=(A(I+1)-A(I))/H-H*(C(I)+H*D(I))
      B(1)=B(2)
      A(1)=A(2)
      B(NA)=B(N)+(2.*C(N)+3.*D(N)*H)*H
      RETURN
  999 DO 2000 I=1,NA
      A(I)=0.
      B(I)=0.
      C(I)=0.
 2000 D(I)=0.
      RETURN
      END
c---------------------------------------------------------------
        FUNCTION SMOO(N,X,C0,C1,C2,C3,XX)
      DIMENSION X(1),C0(1),C1(1),C2(1),C3(1)
      SS=XX-X(N)
      S=XX-X(1)
      IF(S.LE.0.) GOTO 5
      IF(SS.GE.0.) GO TO 6
      DO 1 I=1,N
      T=XX-X(I)
      IF (T.LE.0.)GO TO 2
    1 CONTINUE
    2 D=XX-X(I-1)
      SMOO=C0(I)+D*C1(I)+D*D*C2(I)+D*D*D*C3(I)
      RETURN
    5 SMOO=C0(1)+C1(1)*S
      RETURN
    6 SMOO=C0(N+1)+C1(N+1)*SS
      RETURN
      END
c---------------------------------------------------------------
      SUBROUTINE SZPU2 (N,NP,I)
      I=1
    1 NP=2**I
      IF(NP-N) 2,3,3
    2 I=I+1
      GO TO 1
    3 RETURN
      END
c---------------------------------------------------------------
      SUBROUTINE NLOGN(N,xr,xi,SIGN)
C ALGORITHME COOLEY TUKEY.
C      SIGN=-1.0  EXP (-2*I*PI*...)
C      SIGN =1.0 (1/Q) * EXP(2*I*PI*...)
C      NMAX=PLUS GRANDE VALEUR DE N AVEC
C      DIMENSION X(2**N) ,M(NMAX)
      DIMENSION XR(2),XI(2),M(20)
      LX=2**N
      DO 1 I=1,N
    1 M(I)=2**(N-I)
      DO 4 L=1,N
      NBLOC=2**(L-1)
      LBLOC=LX/NBLOC
      LBHAF=LBLOC/2
      K=0
      DO 4 IBLOC=1,NBLOC
      FK=K
      FLX=LX
      V=SIGN*6.2831853*FK/FLX
      WKR=COS(V)
      WKI=SIN(V)
      ISTAT=LBLOC*(IBLOC-1)
      DO 2 I=1,LBHAF
      J=ISTAT+I
      JH=J+LBHAF
      QR=XR(JH)*WKR-XI(JH)*WKI
      QI=XI(JH)*WKR+XR(JH)*WKI
      XR(JH)=XR(J)-QR
      XI(JH)=XI(J)-QI
      XR(J)=XR(J)+QR
      XI(J)=XI(J)+QI
    2 CONTINUE
      DO 3 I=2,N
      II=I
      IF (K-M(I)) 4,3,3
    3 K=K-M(I)
    4 K=K+M(II)
      K=0
      DO 8 J=1,LX
      IF (K-J) 5,6,6
    6 HOLDR=XR(J)
      HOLDI=XI(J)
      XR(J)=XR(K+1)
      XI(J)=XI(K+1)
      XR(K+1)=HOLDR
      XI(K+1)=HOLDI
    5 DO 7 I=1,N
      II=I
      IF (K-M(I)) 8,7,7
    7 K=K-M(I)
    8 K=K+M(II)
      IF (SIGN) 11,9,9
    9 DO 10 I=1,LX
      XR(I)=XR(I)/FLX
   10 XI(I)=XI(I)/FLX
   11 RETURN
      END
c---------------------------------------------------------------
      subroutine des(x,y,n,xg,xd,yb,yh,ey,pas)
      dimension x(1),y(1)
      ymin=y(1)
      ymax=ymin
      dh=yh-yb
      do 1 i=1,n
      if(y(i).lt.ymin) ymin=y(i)
1     if(y(i).gt.ymax) ymax=y(i)
      if(ey.eq.0.) ey=dh/(ymax-ymin)
      ex=(xd-xg)/(n-1)
      yc=(yb+yh)/2
      ys=0.
      do 3 i=1,n
3     ys=ys+y(i)
      ys=ys/n
      do 2 i=1,n
      y(i)=(y(i)-ys)*ey+yc
2     x(i)=(i-1)*ex
      ycc=yc+(dh)/4
      top=(dh)/20
      xmn=60*ex/pas
      call MoveAbs2(xg,ycc+top)
      call LineAbs2(xg,ycc)
      call LineAbs2(xmn,ycc)
      call LineAbs2(xmn,ycc+top)
      call MoveAbs2(xd,yc-dh)
      call LineAbs2(xd,yc+dh)
      call MoveAbs2(xg,yc)
      call PolyLineAbs2(x,y,n)
      call MoveAbs2(xg,yc)
      call LineAbs2(xd,yc)
      return
      end
c---------------------------------------------------------------
      subroutine szopen(l,c)
      character*80c
c
c Ce sous-programme ouvre l unite logique "l" en lecture
c apres reponse a la question posee dans la chaine de caracteres "c"
c 
      character*80f
1     write(*,c)
      read(*,'(a)')f
      open(l,status='old',err=2,file=f)
      return
2     write(*,'("fichier inexistant")')
      go to 1
      end
c---------------------------------------------------------------
      subroutine szrenv(x,n)
      dimension x(1)
c
c Renversement de l ordre des indices dans le tableau x
c
      m=n/2
      do 1 i=1,m
      is=n-i+1
      xs=x(is)
      x(is)=x(i)
1     x(i)=xs
      return
      end
