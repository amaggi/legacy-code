C  CALCUL DU SPECTRE AU FOYER DES MODES DE RAYLEIGH
c  adapte pour calculer l'excitation avec un pas de 10deg (=360/36)
c  et creer un fichier dessinable par destout
C UTILISE UN FICHIER CREE PAR LE PROGRAMME Takeuchi & Saito (UNITE LOGIQUE 10)
      DIMENSION BCD(20)
      COMMONRAD
      COMMON/C/PROF,T,REP,VPH,VGR,vrupt,icentre,azim
C* AZIMUT COMPTE POSITIVEMENT DU NORD VERS L EST
C* ENTREE
C* PPL PLONGEMENT AXE P ENTRE 0 ET 90 DEGRES
C* PAZ AZIMUT AXE P ENTRE 0 ET 360 DEG
C* TPL PLONGEMENT AXE T
C* TAZ AZIMUT AXE T
C* PROF     PRONFONDEUR DU FOYER
C* ALPHA VITESSE DES ONDES P AU FOYER
C* BETA VITESSE DES ONDES S AU FOYER
C* RO  DENSITE AU FOYER
      character*100 entree, source,outamp,outpha,carte
c
c     Vitesse de rupture sur la faille en km/s
      vrupt=3.
c
312   write(*,'(" nom du fichier de sortie Saito?")')
      read(*,'(a)')entree
      open(10,status='old',form='formatted',err=312,file=entree)
      isort=0
      isort=1
      write(*,'("nom fichier de sortie pour le spectre-source?")')
      read(*,'(a)')source
	outamp=source(:lnblnk(source))//".amp"
	outpha=source(:lnblnk(source))//".pha"
      open(12,file=outamp)
      open(13,file=outpha)
	write(0,'(10a)')" Fichiers de sortie: ",outamp(:lnblnk(outamp))
	write(0,'(10a)')"                     ",outpha(:lnblnk(outpha))
      PI=3.14159265
      RAD=PI/180.
    1 continue
      write(*,'("entrer plongement et azimuth de l axe P en deg")')
      read(*,*) PPL,PAZ
      plp=PPL
      azp=PAZ
      write(*,'("entrer plongement et azimuth de l axe T en deg")')
      read(*,*) TPL,TAZ
      plt=TPL
      azt=TAZ
      write(*,'("entrer la profondeur du foyer (en km)")')
      read(*,*) PROF
      write(*,
     *'("entrer le moment sismique en cgs (NB: 1 Nm=1.e+7 cgs)")')
      read(*,*) amoment
      write(*,'("entrer la distance epicentrale (en km)")')
      read(*,*) dist
c Distance de reference pour le calcul du spectre en km:
      if(dist.eq.0.) dist=10000.
      write(*,'("entrer l azimuth de la station (en degres)")')
	read(*,'(a)') carte
      read(carte,*) AZ
      read(carte,*,end=72) bid, az1, az2, deltaz
	az1=az+az1
	az2=az+az2
	goto 73
   72		continue
		az1=az-180.
		az2=az+180.
		deltaz=10.
   73	continue
      azim=AZ
      write(*,'("entrer la demie duree de la source en sec.")')
      read(*,*)htd
      alf=htd*vrupt*2 
      write(*,'("la longueur de la faille est ",f5.0," km")') alf
      if(alf.ne.0.)then
      write(*,'("temps de montee de la dislocation en sec.?")')
      read(*,*)riset
      write(*,'("voulez-vous calculer le spectre a la source pour le cen
     1troide?")')
      if(nonoui().eq.1)then
      icentre=1
      else
      icentre=0
      write(*,'("vous calculez le spectre a un epicentre suppose etre en
     1 bout de faille...")')
      endif
      endif
      READ (10,'(20a4)')(BCD(I),I=1,20)
	write(carte,'(20a4)') (bcd(i),i=1,20)
      PRINT 102,BCD
      CALL ABR(ALPHA,BETA,RO,PROF)
      PRINT 206,PROF
  206 FORMAT (' PROFONDEUR SEISME',F10.4)
  205 FORMAT(1X,3F10.4)
      PRINT 200,PPL,PAZ
      PRINT 201,TPL,TAZ
      PRINT 202,AZ
  200 FORMAT(' AXE DE COMPRESSION P   PLONGEMENT',F10.4,' AZIMUT',F10.4)
  201 FORMAT (' AXE DE TENSION T   PLONGEMENT',F10.4,' DEG   AZIMUT',F10
     1.4,' DEG')
  202 FORMAT (' AZIMUT STATION',F10.4)
  102 FORMAT(20A4)

      READ(10,*) NMODE
      if(isort.eq.1) then
		write(12,'(a)')carte(:lnblnk(carte))
		write(13,'(a)')carte(:lnblnk(carte))
		endif

	naz=(az2-az1)/deltaz + .5
	daz=deltaz
	nadd=3
	if(az+180..ge.az1.and.az+180..le.az2) nadd=4

      DO 3 KM=1,NMODE
      READ(10,*) NPER
      if(isort.eq.1) then
		write(12,'(i2)')NPER
		write(13,'(i2)')NPER
		endif
      DO 2 K=1,NPER
  204 FORMAT (/' PERIODE',F10.4)

	write(12,*) naz+3,'     mode',km,'     per',k
	write(13,*) naz+nadd,'     mode',km,'     per',k
	write( *,*) 'restent:',nmode-km,' modes',nper-k,' periodes'

	do 1422 iaz=0,naz
	azs=az1+iaz*daz
	az=azs
      AZ=PI/2.-AZ*RAD
      AZO=AZ
c	az0 = angle trigo entre l'Est et la station
      CAZ=COS(AZ)
      SAZ=SIN(AZ)
	if(ideb.eq.0) then
      PPL=PPL*RAD
      PAZ=PI/2.-PAZ*RAD
      TPL=TPL*RAD
      TAZ=PI/2.-TAZ*RAD
      CALL GEOF(PPL,PAZ,TPL,TAZ,DELTA,AMBDA,CK,SK)
	str=90-atan2(sk,ck)/rad
	print *,'strike (+ du N vers l''E)         :',str,' degres'
	print *,'dip    (from horizontal)         :',delta/rad,' degres'
	print *,'slip (counterclockw. from strike):',ambda/rad,' degres'
	ideb=1
	endif
c	K = angle trigo Est-strike
      AZ=ATAN2((-SK*CAZ+CK*SAZ),(CK*CAZ+SK*SAZ))
      STR=(PI/2.-(AZO-AZ))/RAD
cc	print 203, str
  203 FORMAT(' AZIMUT STRIKE + DU N VERS L E',F10.4,' DEG')
      if(iaz.eq.0) CALL ABC(A,B,C,ALPHA,BETA,RO,EO,y1solid,ierr,iend)
      if(ierr.eq.1.or.iend.eq.1)go to 44
c	dans RAY(), AZ doit etre l'angle strike-station
c	compte + dans le sens clockwise (?)
      CALL RAY(DELTA,AMBDA,AZ,A,B,C,EO,y1solid,isort,dist,alf,riset,
     *                              ambdru,amoment,azs)
	if(abs(azs-azim).lt.1e-3) then
		write(12,*) '0.      0.     0.'
		write(13,*) azim,'   0.        0.'
      CALL RAY(DELTA,AMBDA,AZ,A,B,C,EO,y1solid,isort,dist,alf,riset,
     *                              ambdru,amoment,azs)
		endif

 1422	continue

    2 CONTINUE
	write(12,*)'----------- Fin mode',km,'  -----------------------'
    3 CONTINUE
44    continue
  101 FORMAT (I2)
      close(10)
      if(isort.eq.1)then
      write(12,200)plp,azp
      write(12,201)plt,azt
      write(12,206)PROF
      write(12,*)'moment sismique :',amoment,' cgs'
      write(12,*)'distance epicentrale :',dist,' km'
      write(12,202)azim
      write(12,*)'demi-duree de la source : ',htd
      if(icentre.eq.0) then
		write(12,*) 'non-centroide'
		else
		write(12,*) 'centroide'
      endif
      write(12,*)'rise time : ',riset
      write(12,*)'angle(strike,propag. rupt.) : ',ambdru
      close(12)
      endif
      STOP
      END
c------------------------------------------------------------
      SUBROUTINE GEOF(PPL,PAZ,TPL,TAZ,DELTA,AMBDA,CK,SK)
      COMMON RAD
      write(*,'("preferez-vous entrer les vecteurs n et nu?")')
      if(nonoui().eq.1)go to 1
      CPHI=COS(PAZ)
      CKHI=COS(PPL)
      SPHI=SIN(PAZ)
      SKHI=SIN(PPL)
      PX=CKHI*CPHI
      PY=CKHI*SPHI
      PZ=-SKHI
      CPHI=COS(TAZ)
      CKHI=COS(TPL)
      SPHI=SIN(TAZ)
      SKHI=SIN(TPL)
      TX=CKHI*CPHI
      TY=CKHI*SPHI
      TZ=-SKHI
      SQ=SQRT(2.)
      XN=-(PX+TX)/SQ
      YN=-(PY+TY)/SQ
      ZN=-(PZ+TZ)/SQ
      XU=(PX-TX)/SQ
      YU=(PY-TY)/SQ
      ZU=(PZ-TZ)/SQ
      GO TO 2
    1 write(*,'("entrer XN,YN,ZN,XU,YU,ZU")')
      READ(*,*)XN,YN,ZN,XU,YU,ZU
    2 CONTINUE
      PRINT 500,XN,YN,ZN
      PRINT 500,XU,YU,ZU
  500 FORMAT (1X,3F10.4)
  501 FORMAT (6F10.4)
      AKSI=ATAN2(YN,XN)+RAD*90.
      DEL=SQRT(XN*XN+YN*YN)/sqrt(xn*xn+yn*yn+zn*zn)
      DELTA=ASIN(DEL)
      CK=COS(AKSI)
      SK=SIN(AKSI)
      CD=COS(DELTA)
      SD=SIN(DELTA)
      XFU=XU*CK+YU*SK
      YFU=-XU*SK*CD+YU*CK*CD+ZU*SD
      ZFU=XU*SK*SD-YU*CK*SD+ZU*CD
      AMBDA=ATAN2(YFU,XFU)
      DELTAD=DELTA/RAD
      AMBDAD=AMBDA/RAD
      PRINT 203,DELTAD,AMBDAD,ZFU
  203 FORMAT(' INCLINAISON DU PLAN DE FAILLE',F10.4,' DEG',/' ANGLE VECT
     1EUR  BURGER-STRIKE',F10.4,' DEG'/' COMPOSANTE RESIDUELLE VECTEUR B
     1URGER SUR  LA NORMALE AU PLAN DE FAILLE',F10.4)
      RETURN
      END
c------------------------------------------------------------
      SUBROUTINE ABC(A,B,C,ALPHA,BETA,RO,EO,y1solid,ierr,iend)
      COMMON/C/PROF,T,REP,VPH,VGR,vrupt,icentre,azim
      real*8 pi
c
c     Les A, B, C  calcules ici tendent vers (1/2/6371.)* les A,B,C
c     de Harkrider (1970) calcules dans le prog. FORAY plan.
c
      pi=3.14159265359
      ao=6371.
      r=ao-prof
      CALL input(y1solid,y1,y2,y3,y4,eo,ierr,iend)
c
c y2 et y4 doivent etre en 10**10 u. cgs
c
      ak=2*pi/vph/T
      an=ao*ak-0.5
      RIG=RO*BETA*BETA
      ambda=alpha*alpha*ro-2*rig

      a=-an*y3/(2*r)
      b=((3*ambda+2*rig)*y1/r -y2 -an*(an+1)*(3*ambda+2*rig)*y3/(2*r))
     *                  /(an*(ambda+2*rig))
      c=y4/(2*rig)
      RETURN
      END
c-------------------------------------------------------------
      SUBROUTINE RAY(DELTA,AMBDA,AZ,A,B,cc,EO,y1solid,isort,dist,
     1				alf,riset,ambdru,amoment,azss)
      COMMON/C/PROF,T,REP,VPH,VGR,vrupt,icentre,azim
c
c    Les formules (These M. Cara, p. 31 corrigee) donnent, pour une source
c    double couple, impulsion temporelle et moment sismique unite,
c    le deplacement vertical en une station situee a la distance dist
c    et a la profondeur 0. Si on entre le modele en km/s pour les
c    vitesses, en g/cm**3 pour les densites et en km pour les
c    profondeurs, le deplacement obtenu est en km, et pour un moment
c    sismique de 10**25 CGS (dyne.cm). Dans ce programme, on multiplie
c    par le facteur ECH=10**9, qui peut s'interpreter comme donnant un
c    deplacement en microns pour un moment sismique de 10**25 CGS.
c    NB: le 10**25 vient du rapport (km/cm) qui intervient a la
c        puissance 5 dans l'expression du deplacement.
c
c 17-09-1987: les expressions programmees convergent bien vers les
c             formules Rayleigh de Harkrider(1970) et Aki(1980)
c
	c=CC
      r0=6371.
      PI=3.14159265
	rad=pi/180.
      distrad=dist/r0
      ech=1.e+09
c
      SL=SIN(AMBDA)
      CL=COS(AMBDA)
      SD=SIN(DELTA)
      CD=COS(DELTA)
      S2D=2*SD*CD
      C2D=CD*CD-SD*SD
      D0=SL*S2D*B/2
      D1=-SL*C2D*C
      D2=-CL*CD*C
      D3=CL*SD*A
      D4=-SL*S2D*A/2
      AZZ=2*AZ
      S=D1*SIN(AZ)+D2*COS(AZ)
      C=D0+D3*SIN(AZZ)+D4*COS(AZZ)
      AMP=SQRT(S*S+C*C)
      OME=2.*PI/T
      ak=ome/vph
      an=40030./vph/t -0.5

c     amplitude du mouvement vertical pour une source impulsive
c     NB: rep=1/(C*U*(I1+n(n+1)I2)) = 1/I3(Saito)
c	Pour un modele avec couche d'eau, on suppose que le sismogramme
c	est calcule pour une station sur une ile (et non pas a la surface
c	de l'eau), et que le deplacement sur l'ile est identique au
c	deplacement de la 1ere couche solide ==> on multiplie par
c	y1solid au lieu de multiplier par 1 (valeur normalisee en surface).
c	Un calcul sur le modele P7 a montre que cette hypothese est
c	valable a 5-6% pres en deplacement.

      ampz=ech*r0*rep*amp*an*(an+0.5) / (ak*sqrt(2.*an*pi*sin(distrad)))
	ampz=ampz*y1solid
      ampz=ampz*amoment*1.e-25

C** PHIZ PHASE AU FOYER  COMP Z VERS LE HAUT
C** TF EN DT*EXP(-IWT)
      phiz=atan2(s,c)+3.*pi/4.

      if(abs(alf).lt.1e-3) then
c      Pour un echelon (sans rampe): ampz=ampz/ome
c      Pour un echelon (sans rampe): phiz=phiz-pi/2.
         ampz=ampz/ome
         phiz=phiz-pi/2.
      else
c     source rampe temporelle --> (-1+exp(-2iw*hd))/(2w^2*hd)
c     pour TF en exp(-iwt).dt
         hd=alf/vrupt/2.
         spr=-1.+cos(-2.*ome*hd)
         spi=    sin(-2.*ome*hd)
         ampz=ampz*sqrt(spr*spr+spi*spi)/(2.*ome*ome*hd)
         phiz=phiz+atan2(spi,spr)
      endif

C** AMPL AMPLITUDE COMP L
      AMPL=AMPZ*EO
C** PHIL PHASE AU FOYER COMP L  MEMES CONVENTIONS
      PHIL=PHIZ+PI-SIGN(PI/2.,EO)
    1 FORMAT (' PHASE LONGITUDINALE',F10.4,' RD   PHASE VERTICALE',F10.4
     1,' RD')
    4 FORMAT (' AMPLITUDE A',f11.4,' KM POUR UN MOMENT SISMIQUE DE'
     *                    ,e10.3,' DYNES CM ')
	azs=azss
	if(abs(azs-azim-180).lt.1.)
     *	write(13,'(3e14.6,3f10.4)') azs,PHIZ,PHIL,azim,azss,azs
	if(azs.ge.(azim+180.)) azs=azs-360.
      if(isort.eq.1) then
		 write(12,'(4e14.6)')
     *		ampz*sin(azss*rad),ampz*cos(azss*rad),azss,T
		 write(13,'(4e14.6)')
     *		azs,PHIZ,PHIL,T
		endif
    5 FORMAT(' AMPLITUDE VERTICALE ',E10.3,' micron'/' AMPLITUDE LONGITU
     *DINALE ',E10.3,' micron')
      RETURN
      END
c------------------------------------------------------------
      SUBROUTINE ABR (ALPHA,BETA,RO,PROF)
      parameter (ndim=200)
      DIMENSION A(ndim),B(ndim),R(ndim),H(ndim)
      READ (10,*) N
      READ (10,*)(h(I),B(I),R(I),a(I),I=1,N)
      IH=0
    1 IH=IH+1
      IF(PROF.GT.H(IH)) GO TO 1
      ALPHA=A(IH)
      BETA=B(IH)
      RO=R(IH)
      RETURN
      END
c------------------------------------------------------------
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
3     read(*,'(a)')c
      do 1 i=1,14
      if(c.eq.yo(i))nonoui=1
1     continue
      do 2 i=1,8
      if(c.eq.yn(i))nonoui=0
2     continue
      if(nonoui.eq.-1)then
      if(k.ge.3)then
      write(*,'("advienne ce qu il adviendra")')
      return
      else
      k=k+1
      write(*,'("repondez a la question, rappel no",i2)')k
      go to 3
      endif
      else
      return
      endif
      end
c------------------------------------------------------------
      subroutine input(y1solid,yy1,yy2,yy3,yy4,eo,ierr,iend)
      COMMON/C/PROF,T,REP,VPH,VGR,vrupt,icentre,azim
      parameter(NCOUCHMAX=200)
      dimension h(NCOUCHMAX), y1(NCOUCHMAX), y2(NCOUCHMAX),
     1                        y3(NCOUCHMAX), y4(NCOUCHMAX)
      character*80 titre
c
      ierr=0
      iend=0
      lu=10
      nysave=0

      read(lu,'(a)',err=88,end=99) titre
      read(lu,*,err=88,end=99) bid
      read(lu,*,err=88,end=99) bid,bid,t,vph,vgr,rep
      read(lu,*,err=88,end=99) ncouch
c      le deplacement horizontal est: y3*(n+1/2)  = y3*an2   
c      la contrainte horizontale est: y4*(n+1/2)  = y4*an2   
      an2=40030./(vph*t)
c

      do 1 i=1,ncouch
      read(lu,'(i5,5e15.7)',err=88,end=99) ny,h(i),
     *                                     y1(i),y2(i),y3(i),y4(i)
c
c     calcul de l'ellipticite a la premiere couche solide
      if(ny.eq.4.and.nysave.ne.4)then
		eo=-y3(i)*an2/y1(i)
		y1solid=y1(i)
		endif
c
      nysave=ny
    1 continue

      if(prof.lt.h(1))goto 2

      do 3 i=2,ncouch
      if(prof.lt.h(i)) then
          yy1=y1(i-1)+(prof-h(i-1))*(y1(i)-y1(i-1))/(h(i)-h(i-1))
          yy2=y2(i-1)+(prof-h(i-1))*(y2(i)-y2(i-1))/(h(i)-h(i-1))
          yy3=y3(i-1)+(prof-h(i-1))*(y3(i)-y3(i-1))/(h(i)-h(i-1))
          yy4=y4(i-1)+(prof-h(i-1))*(y4(i)-y4(i-1))/(h(i)-h(i-1))
          return
      endif

    3 continue
    2 continue
      print *, 'prof=', prof,'   profmin=',h(1),'  profmax=', h(ncouch)
      stop 'interpolation des fonctions propres impossible'


   88 ierr=1
      return

   99 iend=1
      return

      end
