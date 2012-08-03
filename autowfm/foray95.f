C  CALCUL DU SPECTRE AU FOYER DES MODES DE RAYLEIGH
C UTILISE UN FICHIER CREE PAR LE PROGRAMME Takeuchi & Saito (UNITE LOGIQUE 10)
      DIMENSION BCD(20)
      COMMONRAD
      COMMON/C/PROF,T,REP,VPH,VGR,vrupt,icentre
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
      character*80 entree, disp, source
	character nomprog*69

      call getarg(0,nomprog)

c
c Vitesse de rupture sur la faille en km/s
      vrupt=3.
c
312   write(*,'(" nom du fichier de sortie Saito?",$)')
      read(*,'(a)')entree
      write(*,'(a)')entree
      open(10,status='old',form='formatted',err=312,file=entree)
      isort=0
      write(*,'("voulez-vous creer un fichier pour calcul de synthetique
     1?")')
      if(nonoui().eq.1)then
      isort=1
      write(*,'("nom fichier de sortie pour le spectre-source?",$)')
      read(*,'(a)')source
      write(*,'(a)')source
      write(*,'("nom fichier de sortie pour la dispersion?",$)')
      read(*,'(a)')disp
      write(*,'(a)')disp
      open(11,file=disp)
      open(12,file=source)
      endif
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
      read(*,*) AZ
      azim=AZ
      write(*,'("entrer la demie duree de la source en sec.")')
      read(*,*)htd
      alf=htd*vrupt*2 
      write(*,'("la longueur de la faille est ",f5.0," km")') alf
**      if(alf.ne.0.)then
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
**      endif
      READ (10,'(20a4)')(BCD(I),I=1,20)
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

c	Dans Saito 67 (applique ici) et Harkrider 70, les
c	azimuts croissent dans le sens "trigo" et non pas
c	du Nord vers l'Est ==> passage en repere "trigo"

      AZ=PI/2.-AZ*RAD
      AZO=AZ
      PPL=PPL*RAD
      PAZ=PI/2.-PAZ*RAD
      TPL=TPL*RAD
      TAZ=PI/2.-TAZ*RAD
      CALL GEOF(PPL,PAZ,TPL,TAZ,DELTA,AMBDA,CK,SK)
      CAZ=COS(AZ)
      SAZ=SIN(AZ)
      AZ=ATAN2((-SK*CAZ+CK*SAZ),(CK*CAZ+SK*SAZ))
c	Strike en repere "N >0 E"
      STR=(PI/2.-(AZO-AZ))/RAD
      PRINT 203,STR,delta/rad,ambda/rad
  203 FORMAT(' AZIMUT STRIKE + DU N VERS L E',F10.4,' DEG  dta=',f6.2,
     & '  lda=',f6.2)
**      if(alf.ne.0.)then
      write(*,'("voulez-vous entrer la direction de propagation de la ru
     1pture?")')
      if(nonoui().eq.1)then
      write(*,'("angle en degres entre le strike et la direction",/
     1,",dans le plan de faille, de propagation de la rupture?")') 
      read(*,*) ambdru
      ambdru=rad*ambdru
      else
      write(*,'("la direction de propagation de la rupture est prise 
     1dans la direction du vecteur dislocation")')
      endif
**      endif
      READ(10,*) NMODE
      if(isort.eq.1)write(11,102)BCD
      if(isort.eq.1)write(12,102)BCD
      if(isort.eq.1)write(11,'(i2,f10.1)')NMODE,dist
      if(isort.eq.1)write(12,'(i2,f10.1)')NMODE,dist
      DO 3 KM=1,NMODE
      READ(10,*) NPER
      if(isort.eq.1)write(11,'(i2)')NPER
      if(isort.eq.1)write(12,'(i2,a)')NPER,
     *       '     T    phiz-rad  phil-rad      ampz           ampl'
     *     //'         phiz-deg  phiz-0/99   dg-z      dg-l'
      DO 2 K=1,NPER
  204 FORMAT (/' PERIODE',F10.4)
      CALL ABC(A,B,C,ALPHA,BETA,RO,EO,y1solid,ierr,iend)
      diam=6371.*2.
cc      print *, 'A, B, C =', a*diam,b*diam,c*diam,t,rep,vph,vgr
      if(ierr.eq.1.or.iend.eq.1)go to 44
cc      PRINT 204,T
	print *, nomprog(:lnblnk(nomprog)),'    mode',km,'  periode:',t
      if(isort.eq.1)write(11,'(f10.4,2f10.7)')T,VPH,VGR
      CALL RAY(DELTA,AMBDA,AZ,A,B,C,EO,y1solid,isort,dist,alf,riset,
     *                                       ambdru,amoment)
    2 CONTINUE
    3 CONTINUE
44    continue
      if(isort.eq.1)write(11,'(f10.7)')VPH
  101 FORMAT (I2)
      close(10)
      if(isort.eq.1)then
      write(12,200)plp,azp
      write(12,201)plt,azt
      write(12,206)PROF
      write(12,'(a,e15.6,a)')'  moment sismique :',amoment,' cgs'
      write(12,*)'  distance epicentrale :',dist,' km'
      write(12,202)azim
      write(12,*)'  demi-duree de la source : ',htd
      if(icentre.eq.0) then
		write(12,*) '  non-centroide'
		else
		write(12,*) '  centroide'
      endif
      write(12,*)'  rise time : ',riset
      write(12,*)'  angle(strike,propag. rupt.) : ',ambdru
      close(11)
      close(12)
      endif

      END
c----------------------------------------------------------------------
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
cc cc     PRINT 500,XN,YN,ZN
cc      PRINT 500,XU,YU,ZU
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
cc      PRINT 203,DELTAD,AMBDAD,ZFU
  203 FORMAT(' INCLINAISON DU PLAN DE FAILLE',F10.4,' DEG',/' ANGLE VECT
     1EUR  BURGER-STRIKE',F10.4,' DEG'/' COMPOSANTE RESIDUELLE VECTEUR B
     1URGER SUR  LA NORMALE AU PLAN DE FAILLE',F10.4)
      RETURN
      END
c--------------------------------------------------------------------------
      SUBROUTINE ABC(A,B,C,ALPHA,BETA,RO,EO,y1solid,ierr,iend)
      COMMON/C/PROF,T,REP,VPH,VGR,vrupt,icentre
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
c--------------------------------------------------------------------------
      SUBROUTINE RAY(DELTA,AMBDA,AZ,A,B,C,EO,y1solid,isort,dist,
     1				alf,riset,ambdru,amoment)
      COMMON/C/PROF,T,REP,VPH,VGR,vrupt,icentre
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
	real*8 pi
	save omesave,phizsave,philsave

      r0=6371.
      PI=3.14159265
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
C** TF EN DT*EXP(-IWT) et pi/4 venant du passage Bessel -> cosinus
c   et autres (voir formules Saito 67, Rayleigh vertical des notes JJL)
      phiz=atan2(s,c)+3.*pi/4.
cccc	write(*,*)"azim.relat.=",az*180/pi,"   phases:",
cccc     &     atan2(s,c)*180/pi, phiz*180/pi,(phiz-pi/2)*180/pi

	if(abs(alf).lt.1e-3) then
c      Pour un echelon (sans rampe): ampz=ampz/ome
c      Pour un echelon (sans rampe): phiz=phiz-pi/2.
	   ampz=ampz/ome
	   phiz=phiz-pi/2.
	else
c     source rampe temporelle de duree 2*hd ==>
c       pour TF en exp(-iwt).dt, spectre a multiplier par:
c		(-1+exp(-2iw*hd))/(2w^2*hd) * exp(-iw*td)
c	  avec td= 0  pour un t0 "ondes de volume",
c            td=-hd pour un t0 "centroid"

	   hd=alf/vrupt/2.
	   td=0.
         if(icentre.eq.1) td=-hd
	   spr=-1.+cos(-2.*ome*hd)
	   spi=    sin(-2.*ome*hd)
	   ampz=ampz*sqrt(spr*spr+spi*spi)/(2.*ome*ome*hd)
	   phiz=phiz+atan2(spi,spr)-ome*td
	endif

cXX       if(tlf.ne.0.)then
cXX 
cXX CCC	Pour tenir compte de l'effet Doppler, il faudrait
cXX CCC	preciser quel est le plan de faille qui casse.
cXX CCC	Or, cela est impossible avec les seuls axes P et T.
cXX CCC      if(ambdru.eq.0)ambdru=ambda
cXX CCC      cpsi=cos(ambdru)*cos(az)+sin(ambdru)*sin(az)*cd
cXX CCC      tlf=alf*(1./vrupt-cpsi/vph)
cXX 
cXX C	On tient compte de l'extension temporelle de la source,
cXX C	mais sans effet Doppler
cXX 	tlf=alf/vrupt
cXX       ampz=ampz*(sin(ome*tlf/2.)/(ome*tlf/2.))
cXX       if(icentre.eq.0)phiz=phiz-ome*tlf/2.
cXX 	endif
cXX       if(riset.ne.0)then
cXX 	write(*,*) "Option non disponible dans cette version"
cXX 	write(*,*) "du programme. Contactez votre vendeur pour"
cXX 	write(*,*) "une mise a jour eventuelle (payante)."
cXX CCC      ampz=ampz*(sin(ome*riset/2.)/(ome*riset/2.))
cXX CCC      if(icentre.eq.0)phiz=phiz-ome*tlf/2.
cXX       endif

C** AMPL AMPLITUDE COMP L
      AMPL=AMPZ*EO
C** PHIL PHASE AU FOYER COMP L  MEMES CONVENTIONS
      PHIL=PHIZ+PI-SIGN(PI/2.,EO)
c	Les phases sont remises dans [0,2*pi]
	phiz=mod(phiz+10.*pi,2.*pi)
	phil=mod(phil+10.*pi,2.*pi)
cc      PRINT 1,PHIL,PHIZ
    1 FORMAT (' PHASE LONGITUDINALE',F10.4,' RD   PHASE VERTICALE',F10.4
     1,' RD')
cc      PRINT 4,dist,amoment
    4 FORMAT (' AMPLITUDE A',f11.4,' KM POUR UN MOMENT SISMIQUE DE'
     *                    ,e10.3,' DYNES CM ')
cc      PRINT 5,AMPZ,AMPL
      if(isort.eq.1) write(12,'(f10.4,2f10.6,2e15.7,4f10.2)')
     *                          T,PHIZ,PHIL,AMPZ,AMPL,phiz*180/pi
     *		,phiz*50/pi,-(phiz-phizsave)/(ome-omesave)
     *		,-(phil-philsave)/(ome-omesave)
    5 FORMAT(' AMPLITUDE VERTICALE ',E10.3,' micron'/' AMPLITUDE LONGITU
     *DINALE ',E10.3,' micron')
	phizsave=phiz
	philsave=phil
	omesave=ome
      RETURN
      END
c--------------------------------------------------------------------------
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
c--------------------------------------------------------------------------
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
c--------------------------------------------------------------------------
      subroutine input(y1solid,yy1,yy2,yy3,yy4,eo,ierr,iend)
      COMMON/C/PROF,T,REP,VPH,VGR,vrupt,icentre
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
      if(ny.eq.4.and.nysave.ne.4) then
		eo=-y3(i)*an2/y1(i)
		y1solid=y1(i)
		endif
c
      nysave=ny
    1 continue

c-Modif eric 05/02/2002 pour compilation Linux-----
c     if(prof.lt.h(1))goto 2
      if(prof.lt.h(1))goto 3
c-------------------------------------------------
      do 2 i=2,ncouch
      if(prof.lt.h(i)) then
          yy1=y1(i-1)+(prof-h(i-1))*(y1(i)-y1(i-1))/(h(i)-h(i-1))
          yy2=y2(i-1)+(prof-h(i-1))*(y2(i)-y2(i-1))/(h(i)-h(i-1))
          yy3=y3(i-1)+(prof-h(i-1))*(y3(i)-y3(i-1))/(h(i)-h(i-1))
          yy4=y4(i-1)+(prof-h(i-1))*(y4(i)-y4(i-1))/(h(i)-h(i-1))
          return
      endif

    2 continue
    3 continue
      print *, 'prof=', prof,'   profmin=',h(1),'  profmax=', h(ncouch)
      stop 'interpolation des fonctions propres impossible'


   88 ierr=1
      return

   99 iend=1
      return

      end
