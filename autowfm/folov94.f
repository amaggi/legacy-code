C  CALCUL DU SPECTRE AU FOYER DES MODES DE Love
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
      character*80entree
      character*80disp
      character*80source
	character nomprog*69

      call getarg(0,nomprog)

c
c Vitesse de rupture sur la faille en km/s
      vrupt=3.
c
312   write(*,'(" nom du fichier de sortie Saito?")')
      read(*,'(a)')entree
      open(10,status='old',form='formatted',err=312,file=entree)
      isort=0
      write(*,'("voulez-vous creer un fichier pour calcul de synthetique
     1?")')
      if(nonoui().eq.1)then
      isort=1
      write(*,'("nom fichier de sortie pour le spectre-source?")')
      read(*,'(a)')source
      write(*,'("nom fichier de sortie pour la dispersion?")')
      read(*,'(a)')disp
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
      write(*,'("entrer la profondeur du foyer")')
      read(*,*) PROF
      write(*,*)
     * 'entrer le moment sismique en cgs (NB: 1 Nm=1.e+7 cgs)' 
      read(*,*) amoment
      write(*,'("entrer la distance epicentrale (en km)")')
      read(*,*) dist
c Distance de reference pour le calcul du spectre en km:
      if(dist.eq.0.) dist=10000.
      write(*,'("entrer l azimuth de la station")')
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
      STR=(PI/2.-(AZO-AZ))/RAD
      PRINT 203,STR
  203 FORMAT(' AZIMUT STRIKE + DU N VERS L E',F10.4,' DEG')
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
      if(isort.eq.1)write(12,'(i2)')NPER
      DO 2 K=1,NPER
  204 FORMAT (/' PERIODE',F10.4)
      CALL vg(v,g,ALPHA,BETA,RO,ierr,iend)
      if(ierr.eq.1.or.iend.eq.1)go to 44
cc      PRINT 204,T
	print *, nomprog(:lnblnk(nomprog)),'    mode',km,'  periode:',t
      if(isort.eq.1)write(11,'(f10.4,2f10.7)')T,VPH,VGR
      CALL lov(DELTA,AMBDA,AZ,v,g,isort,dist,alf,riset,ambdru,amoment)
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
      STOP
      END

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
      SUBROUTINE vg(v,g,ALPHA,BETA,RO,ierr,iend)
      COMMON/C/PROF,T,REP,VPH,VGR,vrupt,icentre
      real*8 pi
c
c     Les v, g  calcules ici tendent vers ak* les v, g
c     de Harkrider (1970) calcules dans le prog. FOLOV plan.
c     Le signe de G dans Harkrider(1970) semble faux, bien
c     qu'il soit cense corriger un precedent "sign misprint"
c
      pi=3.14159265359
      ao=6371.
      r=ao-prof
      CALL input(y1,y2,ierr,iend)
c
c y2 doit etre en 10**10 u. cgs
c
      ak=2*pi/vph/T
      an=ao*ak-0.5
      RIG=RO*BETA*BETA

      v=an*y1/r
      g=-y2/rig
      RETURN
      END
      SUBROUTINE lov(DELTA,AMBDA,AZ,v,g,isort,dist,alf,riset
     1,ambdru,amoment)
      COMMON/C/PROF,T,REP,VPH,VGR,vrupt,icentre
c
c    Les formules donnent, pour une source
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
c             formules Love de Harkrider(1970) et Aki(1980)
c             si on change le signe de G dans Harkrider(1970)
c
	real*8 pi
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
      D0=0.
      D1=cl*cd*g
      D2=-sl*c2d*g
      D3=sl*s2d*v/2.
      D4=cl*sd*v
      AZZ=2*AZ
      S=D1*SIN(AZ)+D2*COS(AZ)
      C=D0+D3*SIN(AZZ)+D4*COS(AZZ)
      AMP=SQRT(S*S+C*C)
      OME=2.*PI/T
      ak=ome/vph
      an=40030./vph/t -0.5

c     amplitude du mouvement transversal pour une source impulsive
c     NB: rep=1/(C*U*I1) = 1/I3(Saito)
      ampt=ech*r0*rep*amp*(an+0.5)*(an+0.5) /
     *         (2.*(an+1.)*ak*sqrt(2.*an*pi*sin(distrad)))
      ampt=ampt*amoment*1.e-25

C** PHIt PHASE AU FOYER  COMP transversale, Z >0 VERS LE HAUT
c   L,T,Z triedre direct, L >0 de la source vers la station
C** TF EN DT*EXP(-IWT)
      phit=atan2(s,c)+5.*pi/4.

      if(abs(alf).lt.1e-3) then
c      Pour un echelon (sans rampe): ampt=ampt/ome
c      Pour un echelon (sans rampe): phit=phit-pi/2.
        ampt=ampt/ome
        phit=phit-pi/2
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
         ampt=ampt*sqrt(spr*spr+spi*spi)/(2.*ome*ome*hd)
         phit=phit+atan2(spi,spr)-ome*td
      endif

c     Les phases sont remises dans [0,2*pi]
      phit=mod(phit+10.*pi,2.*pi)

cc      PRINT 1,PHIt
    1 FORMAT (' PHASE transversale',F10.4,' RD')
cc      PRINT 4,dist,amoment
    4 FORMAT (' AMPLITUDE A',f11.4,' KM POUR UN MOMENT SISMIQUE DE'
     *                     ,e10.3,' DYNES CM ')
cc      PRINT 5,AMPt
      if(isort.eq.1)write(12,'(f10.4,2f10.6,2e15.7,f10.2)')
     &                         T,PHIt,phit,ampt,ampt,phit*180/pi
    5 FORMAT(' AMPLITUDE transversale ',E10.3,' micron')
      RETURN
      END
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
      subroutine input(yy1,yy2,ierr,iend)
      COMMON/C/PROF,T,REP,VPH,VGR,vrupt,icentre
      parameter(NCOUCHMAX=200)
      dimension h(NCOUCHMAX), y1(NCOUCHMAX), y2(NCOUCHMAX)
      character*80 titre
c
      ierr=0
      iend=0
      lu=10

      read(lu,'(a)',err=88,end=99) titre
      read(lu,*,err=88,end=99) bid
      read(lu,*,err=88,end=99) bid,bid,t,vph,vgr,rep
      read(lu,*,err=88,end=99) ncouch
      an2=40030./(vph*t)
c

      do 1 i=1,ncouch
      read(lu,'(i5,5e15.7)',err=88,end=99) ny,h(i), y1(i),y2(i)
c
    1 continue
c-modif eric 5/02/2002 pour compilation sous Linux------
c     if(prof.lt.h(1))goto 2
      if(prof.lt.h(1))goto 3
c------------------------------------------------------
      do 2 i=2,ncouch
      if(prof.lt.h(i)) then
          yy1=y1(i-1)+(prof-h(i-1))*(y1(i)-y1(i-1))/(h(i)-h(i-1))
          yy2=y2(i-1)+(prof-h(i-1))*(y2(i)-y2(i-1))/(h(i)-h(i-1))
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
