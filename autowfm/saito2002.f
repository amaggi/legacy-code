c            ()main3.f:1.2    3/7/86
c     saito93.f:  identique a saito90.f sauf realignement des
c                 variables dans les common (i4 apres les d8)
c                 et ajout de virgules dans les formats
c
c     Version 88: doit creer les fichiers deplacements-tensions et
c                 derivees partielles directement lisibles par les
c                 programmes de Waveform Modelling.
c
*     VERSION SUN 2/120  24 FEVRIER 1986
*     D'APRES VERSION UNIVAC-CRONENBOURG SUR BANDE.
*     FICHIERS SAITO.MAINJJL2 + SAITO.MANTEAU + SAITO.SSPM
 
*    MISE EN CONFORMITE DE TOUS LES COMMON /MODELE/
*    MODIF DES INSTRUCTIONS:  REAL FUNCTION TATA*8(A,B,...)
*    MODIF DES INSTRUCTIONS:  REAL*8 TOTO/3.14159D+0/
*    VERSION DOUBLE PRECISION: J'AI TOUT MIS EN REAL*8
 
C PROGRAM FOR DISPERSION PROBLEMS
C
C DOUBLE PRECISION
C
C SUBROUTINE REQUIRED - INIVAL, CONT, DIFCOE, INTGND, DENERG, DCDP
C                       THESE ARE PROBLEM DEPENDENT
C                       WNFRQ, CHAREQ ARE ALSO PROBLEM DEPENDENT,
C                       BUT INCLUDED IN MAIN ROUTINES
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/MODEL /
     1              RAD,TMASS,GS,RHO0,GA,
     2              HI(500),RHOI(500),VPI(500),VSI(500),NAME(20),
     3              XII(500),PHII(500),ETAI(500),qbetai(500),
     4              D(1000),H(1000),RHO(1000),ELAMB(1000),EMU(1000),
     5              GR(1000),XI(1000),PHI(1000),ETA(1000),qbeta(1000),
     6              DMAX,ICONT,NLAY,ISO,LAME,NDIV,ILAY,idpmi,idpmax
      COMMON/VALUE /MODE,IERROR,I,ISTEP,J,K,L,N,NSOL,ISUM,KMAX,
     1              NMAX,IBOTM,ITOP,LY,LD,ACR,ELLIP,
     2              RMAX,RC,HC,H1,H2,ALP,WN,WN2,C,C2,FRQ,FRQ2,T,
     3              U,ENGY,DELTA,
     4              A(6,6),F(20),Y(6),repon,njob,nper,jyy,jdd,npar,
     5			ncouc,hdeb,hmin,hmax
      COMMON/SOL   /YN(6,1000),YB(6,3,20),SUM(20),Q(3,21),
     1              WNB(100),TT(100),CC(100),UU(100),ENG(100),ELL(100),
     2              AC(100)
      DIMENSION  PERD(100),CMN(100),DC(100),CMX(100),DEP(100),EPSI(100),
     1           JYI(100),JDERI(100)
      DIMENSION  EP(100), APRM(100,7)
      character nomfich*80, carte*80
C
  33  write(*,*) ' Entrer le nom du fichier-modele:'
      read(*,'(a)') nomfich
      open(12,file=nomfich,err=33,status='old')
  34  write(*,*) ' Entrer le nom du fichier-listing:'
      read(*,'(a)') nomfich
      open(14,file=nomfich,err=34)
      READ(12,1)  ISET,JSET
    1 FORMAT(3I5)
C ISET=INPUT DATA SET,  JSET=OUTPUT DATA SET
C     IF(ISET.GT.7)  REWIND  ISET
      WRITE(14,601)
C     IF(JSET.GT.7)  REWIND  JSET
CCC   CALL  INTIME(.TRUE.,.TRUE.,.FALSE.,.FALSE.,'  DISP  ',ITIME)
C
C
      READ(12,*)  MJOB
   10 FORMAT(I5)
C
C MJOB
C
      DO  9000  MJ=1,MJOB
      WRITE(14,901)
  901 FORMAT(1H1)
C
      CALL  RDMDL(ISET,JSET)
      RMAX=RAD
C
      WRITE(14,910)  NAME,NDIV,DMAX,RAD,NLAY,LAME,ISO,TMASS,GS,RHO0
  910 FORMAT(1H ,20A4,/,
     1       1H ,5HDIV.=,I5,2X,5HUP TO,F10.3,/,
     2       1H ,5HR   =,F9.3,/,
     3       1H ,5HLAY.=,I5,/,1H ,5HLAME=,I5,/,1H ,5HISO =,I5,/,
     4       1H ,5HMASS=,D15.5,/,1H ,5HGS  =,F11.5,/,1H ,5HRHO =,F11.5)
C
      CALL  MLIST
C
      read(12,'(a)') carte
      read(carte,*) jyy, jdd

      if(jyy.ne.0) then
 333  write(*,*) ' Entrer le nom du fichier-deplacements:'
      read(*,'(a)') nomfich
      open(13,file=nomfich,err=333,status='new')
      endif

      if(jdd.ne.0) then
ccc	write(*,*) ' entrer les indices min et max des couches'
ccc	write(*,*) ' auxquelles on ecrit les deriv. part.  '
ccc	write(*,*) ' (la couche min doit etre plus bas que toutes les'
ccc	write(*,*) '  discontinuites du modele)'
	read(carte,*) bid, bid, idpmi, idpmax
 334  write(*,*) ' Entrer le nom du fichier-der.part.:'
      read(*,'(a)') nomfich
      open(15,file=nomfich,err=334,status='new')
	write(*,*) ' ( verifier que la couche min est plus bas que'
	write(*,*) '   toutes les discontinuites du modele)'
	write(*,*)
	
        npar=7
	endif

	write(*,'(2a)') '       perio  C-top   C-step  C-end  ',
     *			'zbtm         C      U    hmin hmax'
      READ(12,*)  NJOB
  920 FORMAT(I5)
      NY=0
      ND=0
      IERROR=0
C
C NJOB
C
      DO  8000  NJ=1,NJOB
      read(12,*)  nper
      if(nj.ne.1.and.jyy.ne.0) write(13,*) nper
      if(nj.ne.1.and.jdd.ne.0) write(15,*) nper
      do  8000  npm=1,nper
*     READ(12,800)  NP,DEPTHO,EPSO,MODE,JYO,JDERO
  800 FORMAT(I5,F10.0,E10.3,3I2)
C NP=NUMBER OF POINTS IN ONE NJOB
C DEPTHO,EPSO,JYO,JDERO=DEFAULT VALUES FOR DEPTH,EPSIL,JY,JD
*     READ(12,802)  FMT
  802 FORMAT(20A4)
C***** CHANGED INPUT FORMAT BY BOB NORTH , ELIMINATING NEED FOR CARDS DI
C***** ONLY IN PERIOD.
CX    READ(12,FMT)  (PERD(NT),CMN(NT),DC(NT),CMX(NT),DEP(NT),EPSI(NT),
CX   1              JYI(NT),JDERI(NT),NT=1,NP)
* 810 READ(12,FMT) TSTART,TSTEP,TEND,CTOP,CSTEP,CEND,DEPMX
  810 READ(12,*) NP,DEPTHO,EPSO,MODE,JYO,JDERO,
     *           TSTART,TSTEP,TEND,CTOP,CSTEP,CEND,DEPMX
      if(jyy.ne.0) jyo=2
      if(jdd.ne.0) jdero=2
      T1=CTOP
      CMN(1)=T1
	cmx(1)=cend
	if(cend.eq.0.) then
      IF(MODE.LE.4) THEN
      CMX(1)=T1*3.
      ELSE
      CMX(1)=T1/3.
      ENDIF
	endif
c on cherche une solution comprise entre CMN et CMX (=3*CMN pour MODE<5)
c    par increment initial de CMN/30 = DC
c  NT correspond aux modes, la periode etant fixe si TSTEP=0
c a l'etape NT+1, CMN est pris egal a la solution de l'etape NT
c                 augmente d'un pourcentage defini par TEND
c
      DO 830 NT=1,NP
      PERD(NT)=TSTART+(NT-1)*TSTEP
      DEP(NT)=DEPMX
      EPSI(NT)=EPSO
      JYI(NT)=JYO
  830 JDERI(NT)=JDERO
  840 CONTINUE
C***** END OF INSERTION BY BOB NORTH
      LY=0
      LD=0
      IF(IERROR.GT.1)  GO TO  8000
C
C PJOB
C
      DO  7000  NT=1,NP
	ntm=nt
	if(np.eq.1) ntm=nj
      IF(IERROR.LE.1)  IERROR=0
	dc(nt)=cstep
	if(cstep.eq.0.) DC(NT)=CMN(NT)/30.
      WRITE(14,187)ntm,'.',npm,': ',deptho,perd(nt),CMN(NT),CMX(NT),
     *              ' INC:', DC(NT)
  187 FORMAT(i3.2,A,I2.2,A,f6.0,f6.2,2F8.4,A,F7.4)
      WRITE(*,188)ntm,'.',npm,':',perd(nt),CMN(NT),dc(nt),CMX(NT)
  188 FORMAT(i2.2,A,I2.2,A,f6.1,3F8.4," ",$)
      WRITE(66,189)ntm,'.',npm,':',deptho,epso,mode,0,0,perd(nt),
     *			tstep,tend,CMN(NT),dc(nt),CMX(NT),depmx
  189 FORMAT(i2.2,A,I2.2,A,f6.0,e9.2,3i2,
     *		f6.1,f3.0,f6.3,3F8.4,f3.0," ",$)
      DEPTH=DEP(NT)
      EPSIL=EPSI(NT)
      JY   =JYI(NT)
      JD   =JDERI(NT)
C JY=PARAMETER FOR YLIST,  JD=PARAMETER FOR DLIST
      IF(DEPTH.EQ.0.0)  DEPTH=DEPTHO
      IF(EPSIL.EQ.0.0)  EPSIL=EPSO
      IF(JY.EQ.0)       JY=JYO
      IF(JD.EQ.0)       JD=JDERO
C
      CALL  TTL1(MODE)
C
      CALL  TABLE(IBOTM,DEPTH,ndiv,idisc)
ccc	print *
ccc	print *, 'prof. demandee:', depth,'   obtenue:',d(ibotm),'(',
ccc     *	ibotm,')'
	write(*,'(f5.0," ",$)') d(ibotm)
	call flush(6)
      IF(IBOTM.LE.1)  GO TO  7000
      DEPTH=D(IBOTM)
C
      IF(PERD(NT))  730,730,710
  710 CONTINUE
      IF(MODE-5)  720,730,730
  720 CONTINUE
      WRITE(14,721)  PERD(NT),DEPTH
  721 FORMAT(1H ,2HT=,F10.5,3X,6HDEPTH=,F8.3,/,1H ,5X,1HC,14X,5HDELTA)
      GO TO  750
  730 CONTINUE
      APERD=ABS(PERD(NT))
      WRITE(14,731)  APERD,DEPTH
  731 FORMAT(1H ,2HN=,E12.5,3X,6HDEPTH=,F8.3,/,1H ,5X,1HT,14X,5HDELTA)
  750 CONTINUE
C
C TO FIND A ROOT OF CHARACTERISTIC EQUATION
C
      ITOP=1
      NMAX=0
C
      CALL INTERP(XX,PERD(NT),CMN(NT),DC(NT),CMX(NT),EPSIL,JUMP)
C
      IF(JUMP)  6000,7100,6000
 7100 CONTINUE
      IF(IERROR-2)  7900,8500,9900
C
C TO COMPUTE EIGENFUNCTIONS
C
 6000 CONTINUE
      J=1
      IF((MODE.EQ.2).OR.(MODE.EQ.4))  CALL  INTEG
      J=2
      LY=LY+1
      DEP(LY)=DEPTH
C
      CALL  INTEG
C
C PRINT AND/OR PUNCH EIGENFUNCTIONS
C
      CALL  YLIST(JY,JSET,idisc)
      WRITE(14,601)
  601 FORMAT(1H )
C
C DERIVATIVES W.R.T. PARAMETERS
C
      IF(JD.EQ.0)  GO TO  7900
C
      CALL  DERIV
      LD=LD+1
C
C PRINT AND/OR PUNCH DERIVATIVES
C
      CALL  DLIST(JD,JSET,ep,aprm,incr)
      WRITE(14,601)
C
 7900 WRITE (14,700)
  700 FORMAT (/1X,70(1H*),/)
C
C END OF PJOB
C
	if(jump.eq.1) then
		write(*,'(a,2f7.4,2f6.0)') 'ok ', c, u, hmin, hmax
		write(66,'(a,2f7.4,2f6.0)') 
		else
		write(*,'(a,f7.4,3x,2f6.0)') 'PAS OK ',xx, hmin, hmax
		write(66,'(3f8.4,f3.0,a,f7.4)') xx-0.01,cstep,0.,0.,
     *                                          'PAS OK ',xx
		endif
C     XX = RESULTAT DU CALCUL EN COURS, SI OK (JUMP=1)
C        = DERNIERE VALEUR OBTENUE SI DEPASST TAILLE DANS INTERP
C        = 2EME BORNE DE L'INTERVALLE DE RECHERCHE SI 'NO ROOT'
          T1=XX
      IF(MODE.LE.4) THEN
C     T1=CC(NT)
C     IF(JUMP.EQ.0) T1=XX
      CMN(NT+1)=T1*(1.+TEND)
      CMX(NT+1)=CMN(NT+1)*3.
      ELSE
C     T1=TT(NT)
C     IF(JUMP.EQ.0) T1=XX
      CMN(NT+1)=T1*(1.-TEND)
      CMX(NT+1)=CMN(NT+1)/3.
      ENDIF
 7000 CONTINUE
C
 8500 CONTINUE
      IF(LY.EQ.0)  GO TO  8900
      NY=NY+LY
      ND=ND+LD
      WRITE(14,901)
C     WRITE(14,910)  NAME,NDIV,DMAX,RAD,NLAY,LAME,ISO,TMASS,GS,RHO0
      CALL  TTL1(MODE)
      WRITE(14,851)
  851 FORMAT(1X,4HK(N),5X,1HT,8X,1HC,6X,1HU,7X,3HK.E,4X,
     1 5HY3/Y1,3X,5HERROR,4X,5HDEPTH)
      WRITE(14,852) (WNB(LYY),TT(LYY),CC(LYY),UU(LYY),ENG(LYY),ELL(LYY),
     1               AC(LYY),DEP(LYY),LYY=1,LY)
  852 FORMAT(1X,F4.0,F11.5,2F7.4,E10.4,F6.4,E10.4,F7.1)
C
 8900 CONTINUE
      WRITE(14,601)
C
C END OF NJOB
C
 8000 CONTINUE
      if(jdd.ne.0) then
      WRITE(15,'(8F10.4)') hdeb, hmin, (EP(I),I=2,ncouc)
      DO 683 J=1,7
      WRITE(15,'(6F10.6)') (APRM(I,J),I=1,ncouc)
  683 CONTINUE
      endif
C
 9100 CONTINUE
C     IF(JSET.LE.7)  GO TO  9200
      GO TO 9200
 9201 WRITE(14,930)  NAME
  930 FORMAT(1H /1H ,20A4)
      WRITE(14,931)  NY,ND
  931 FORMAT(1H ,I3,2X,24HSET(S) OF EIGENFUNCTIONS/,
     1       1H ,I3,2X,21HSET(S) OF DERIVATIVES)
 9200 CONTINUE
C
C END OF MJOB
C
 9000 CONTINUE
 9900 CONTINUE
C     IF(ISET.GT.7)  REWIND  ISET
C     IF (JSET.GT.7) ENDFILE JSET
C     IF(JSET.GT.7)  REWIND  JSET
      WRITE(14,601)
C     CALL  INTIME(.FALSE.,.FALSE.,.TRUE.,.TRUE.,'  DISP  ',ITIME)
      write(0,*) ' Le fichier correspondant a l''unite 66 contient'
      write(0,*) ' les informations necessaires a un nouveau run'
      write(0,*) ' (trier par "sort" et changer les "pas OK").'
      STOP
      END
c-----------------------------------------------------------------------
      SUBROUTINE  listm(ndisc2)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/MODEL /
     1              RAD,TMASS,GS,RHO0,GA,
     2              HI(500),RHOI(500),VPI(500),VSI(500),NAME(20),
     3              XII(500),PHII(500),ETAI(500),qbetai(500),
     4              D(1000),H(1000),RHO(1000),ELAMB(1000),EMU(1000),
     5              GR(1000),XI(1000),PHI(1000),ETA(1000),qbeta(1000),
     6              DMAX,ICONT,NLAY,ISO,LAME,NDIV,ILAY,idpmi,idpmax
      COMMON/VALUE /MODE,IERROR,I,ISTEP,J,K,L,N,NSOL,ISUM,KMAX,
     1              NMAX,IBOTM,ITOP,LY,LD,ACR,ELLIP,
     2              RMAX,RC,HC,H1,H2,ALP,WN,WN2,C,C2,FRQ,FRQ2,T,
     3              U,ENGY,DELTA,
     4              A(6,6),F(20),Y(6),repon,njob,nper,jyy,jdd,npar,
     5			ncouc,hdeb,hmin,hmax
      COMMON/SOL   /YN(6,1000),YB(6,3,20),SUM(20),Q(3,21),
     1              WNB(100),TT(100),CC(100),UU(100),ENG(100),ELL(100),
     2              AC(100)
C
      if(jyy.ne.0) then
      WRITE(13,201)  NAME
  201 FORMAT(20A4)

c     ndisc = nb de discontinuites dans le modele
***      ndisc=2*nmax-ibotm-1
c  On a: nlay=ilay+(ndiv-1)*(ilay-1-ndisc) ==>
	 ndisc=ilay-1-(nlay-ilay)/(ndiv-1)
c     indice=(indice initial)*NDIV - (NDIV-1)*(1+nb de profondeurs doublees)
c	puis on prend une couche sur 2 ==>
	if((ndiv/2)*2.ne.ndiv) stop 'NDIV doit etre multiple de 2'
**	avant Love: ncou2=idpmax*ndiv/2-(ndiv/2-1)*(1+ndisc)
CJJ mai 2002 ncou2=(idpmax*ndiv-(ndiv-1)*(1+ndisc)-itop-ndisc2)/2+ndisc2+1
CJJ mai 2002 formule ci-dessus incompréhensible,
CJJ    on reprend l'ancienne formule (avant Love) et on tient compte par nwater
CJJ    de la suppression de la couche d'eau dans le cas Love (i.e. ndisc2.ne.ndisc)
CJJ    Pour Rayleigh, pas de correction, donc on met nwater=0
       nwater=0
       if(ndisc2.ne.ndisc) then
         do 2002 ibz=1,ilay
         if(abs(vsi(ibz)).gt. 1e-6) goto 2003
         nwater=ibz
 2002   continue
2003   continue
       endif
       ncou2=(idpmax-nwater)*ndiv/2-(ndiv/2-1)*(ndisc2+1)
*REM -------------------------------------------------------------------
*REM  ilay= nb de couches du modele fourni en entree
*REM  nlay= nb de couche du modele apres subdivision par ndiv
*REM  itop,ibotm= indices des couches sup et inf de l'integration
*REM  itop=1 sauf pour Love dans un modele avec eau
*REM  ibotm= point-modele (apres ndiv) de début d'integration,
*REM         qui depend donc du mode considere.
*REM  nmax= meme prof que ibotm, mais compte en ne gardant qu'un pt sur deux
*REM  ndep= nb de lignes de la fonction deplacement ecrites sur LU 13
*REM  ndisc= nb de discontinuites dans le modele fourni
*REM       = nb de couches d'épaisseur nulle, sauf la surface
*REM  idisc= nb de discontinuites jusqu'a la couche ibotm (du mode)
*REM  ndisc2= nb de discontinuites dans le modele apres suppression
*REM          eventuelle de la couche d'eau (pour Love)
*REM  ncou2=  nb de couches du modele liste en tete du fichier-deplacement
*REM       =  de "itop" a "idpmax" en prenant une couche sur deux 
*REM -------------------------------------------------------------------
*REM  Pour les cas autres que Love dans un modele avec eau,
*REM  ndisc2=ndisc et itop=1 donc l'ancien ncou2 vaut:
*REM    (idpmax*ndiv-(ndiv-1)*(1+ndisc)-itop-ndisc2)/2+ndisc2+1
*REM  = idpmax*ndiv/2 - ((ndiv-1)(ndisc+1)+itop+ndisc)/2 + ndisc+1
*REM  = idpmax*ndiv/2   - ( ndiv  * (ndisc+1) ) /2         + ndisc+1
*REM  = (idpmax-ndisc-1)*ndiv/2+ndisc+1,  qui est la bonne formule pour ncou2
*REM     et le modele mis dans le fichier deplacement etait donc correct.
*REM     Par contre, il était faux pour Love dans un modele avec eau.
*REM -------------------------------------------------------------------
*     write(13,*) "ilay,nlay,ndiv,ndisc,ndisc2,idpmi,idpmax,",
*    &            "itop,ibotm,nmax,ncou2,ncou2.old"
*     write(13,*)  ilay,nlay,ndiv,ndisc,ndisc2,idpmi,idpmax, 
*    &             itop,ibotm,nmax,ncou2,
*    &(idpmax*ndiv-(ndiv-1)*(1+ndisc)-itop-ndisc2)/2+ndisc2+1
*REM -------------------------------------------------------------------
*REM   Dans listm, on print le modele sur LU 13 en commençant à itop
*REM   et en comptant ncou2 couches, pour arriver a idpmax
*REM   Dans YLIST, on print sur LU 13 les déplacements de itop a ibotm,
*REM   mais le test est fait sur NMAX plutot que sur ibotm
*REM   Dans DLIST, on print sur LU 15 les der.part de idpmi a idpmax sans subdivision
*REM -------------------------------------------------------------------

c      write(13,*) NMAX  pas bon: NMAX depend des profs d'integration

      write(13,*) ncou2

      N=0
      I=ITOP-1
C
C
  500 CONTINUE
      I=I+1
C
  600 CONTINUE
      N=N+1
      DD=D(I)
      beta= emu(i)/rho(i)
      alpha= elamb(i)/rho(i) +2.*beta
      alpha=sqrt(alpha)
      beta=sqrt(beta)
      ro=rho(i)
      WRITE(13,801)  DD,beta,ro,alpha
  801 FORMAT(5E15.7)
C
c     NMAX  pas bon: NMAX depend des profs d'integration
**     IF(N.GE.NMAX)  GO TO  1000
      IF(N.GE.ncou2)  GO TO  1000
      IF(H(I+1).EQ.0.0)  GO TO  500
      I=I+2
      GO TO  600
C
 1000 CONTINUE

      write(13,*) njob
      write(13,*) nper
      endif

      RETURN
      END
c-----------------------------------------------------------------------
C   PROGRAMME LFOR.SAITOSSPM
C-----------------------------------------------------------------------
      SUBROUTINE  RDMDL(ISET,JSET)
      IMPLICIT REAL*8(A-H,O-Z)
C
C PURPOSE -   TO READ AN EARTH MODEL.
C
C             THIS PROGRAM READS AN EARTH MODEL, DIVIDES ONE STEP INTO
C     NDIV STEPS UP TO DMAX DEPTH, AND CONVERT VPI AND VSI TO LAMBDA
C     AND MU.
C
C MODEL -     AN EARTH MODEL IS SPECIFIED BY SIX PHYSICAL PARAMETERS,
C     RHOI, VPI, VSI, XII, PHII, AND ETAI AT SUCCESSIVE DEPTHS.
C     WHEN THE MODEL IS ISOTROPIC (ISO=0) XII, PHII, AND ETAI SHOULD NOT
C     BE GIVEN.
C             DEPTH INFORMATION IS GIVEN BY THE ARRAY HI.
C     THE FIRST HI CONTAINS THE DEPTH (USUALLY ZERO) THAT CORRESPONDS TO
C     THE PHYSICAL PARAMETERS GIVEN BY THE FIRST CARD.
C     I-TH ELEMENT OF HI,E.G. HI(I), SPECIFIES THE DISTANCE BETWEEN
C     (I-1)-TH AND I-TH CARD.     IF(HI(I).EQ.0) TWO ORDINATES ARE
C     INTERPRETED AS UPPER(I-1) AND LOWER(I) SIDES OF A DISCONTINUITY.
C
      COMMON/MODEL /
     1              RAD,TMASS,GS,RHO0,GA,
     2              HI(500),RHOI(500),VPI(500),VSI(500),NAME(20),
     3              XII(500),PHII(500),ETAI(500),qbetai(500),
     4              D(1000),H(1000),RHO(1000),ELAMB(1000),EMU(1000),
     5              GR(1000),XI(1000),PHI(1000),ETA(1000),qbeta(1000),
     6              DMAX,ICONT,NLAY,ISO,LAME,NDIV,ILAY,idpmi,idpmax
C**** FOLLOWING CARDS HAVE BEEN INSERTED OR CHANGED BY BOB NORTH********
C     DIMENSION DR(50),R(50),EL(50),EM(50)
C**** END OF INSERTION/REPLACEMENT *************************************
      DIMENSION  FMT(20)
C
C INPUT -     ILAY  =NO. OF ORDINATES(DEPTHS) AT WHICH PARAMETERS ARE
C             GIVEN
C             RAD   =RADIUS OF THE EARTH IN KM
C             TMASS =TOTAL MASS OF THE EARTH IN C.G.S.
C             GS    =SURFACE VALUE OF THE GRAVITY IN C.G.S.
C             RHO0  =DENSITY IN THE INNERMOST HOMOGENENEOUS CORE
C             ISO   =0 WHEN THE MODEL IS ISOTROPIC
C             ICONT=1 M0DELE CONTINU
C             ICONT=0 MODELE A COUCHES HOMOGENES
C             LAME  =0 WHEN VPI AND VSI ARE VELOCITIES.
C             IF(LAME.NE.0)  THEN VPI AND VSI ARE INTERPRETED AS LAME'S
C             CONSTANTS.
C             HI    =DISTANCES BETWEEN TWO ORDINATES IN KM
C             RHOI  =DENSITIES IN C.G.S.
C             VPI   =COMPRESSIONAL WAVE VELOCITIES IN KM/SEC OR
C                    LAMBDA                        IN (10**10) C.G.S.
C             VSI   =SHEAR WAVE VELOCITIES         IN KM/SEC  OR
C                    MU                            IN (10**10) C.G.S.
C             IF ISO=1 (ANISOTROPIC MODEL) VPI=VPH,VSI=VSV
C             XII   =(C11-C12)/2C44
C             PHII  =C33/C11
C             ETAI  =C13/(C11-2C44)
c             qbetai= inverse du facteur de qualite
C             JSET  =DATA SET REF. NO. FOR OUTPUT
C             ISET  =DATA SET REF. NO. FOR INPUT
C
C OUTPUT -    NLAY  =NO. OF ORDINATES IN DIVIDED MODEL
C             D     =DEPTHS IN KM
C             H     =DISTANCES IN KM     H(1)=D(1),  H(I)=D(I)-D(I-1)
C             RHO   =DENSITIES IN C.G.S.
C             ELAMB =LAMBDA IN (10**10) C.G.S.
C             EMU   =MU     IN (10**10) C.G.S.
C             GR    =GRAVITIES IN (10**5) C.G.S.
C             XI    =XI
C             PHI   =PHI
C             ETA   =ETA
c             qbeta = inverse du facteur de qualite
C SUBROUTINE REQUIRED - DIVIDE, GRAVTY
C
C     IF(ISET-7)  100,100,500
      GO TO 100
C
C CARD INPUT
C
  100 CONTINUE
      READ(12,110)  NAME
  110 FORMAT(20A4)
      READ(12,120)  ILAY,RAD,TMAS ,GS,RHO0,ISO,LAME
  120 FORMAT(I5,F10.2,E15.7,2F10.0,2I2)
      IF (RAD.EQ.0.) RAD=6371.0
      IF (GS.EQ.0.) GS=982.0
      TMASS=TMAS
      READ(12,130)  FMT
  130 FORMAT(20A4)
      READ(12,135) ICONT
  135 FORMAT(I5)
      IF(ISO.EQ.1) GO TO 160
  140 CONTINUE
      READ(12,FMT) (HI(I),RHOI(I),VPI(I),VSI(I),
     * qbetai(i),I=1,ILAY)
      DO 150 I=1,ILAY
      XII(I)=0.0
      PHII(I)=0.0
      ETAI(I)=0.0
      qbetai(i)=0.0
  150 CONTINUE
      GO TO 170
  160 CONTINUE
      READ(12,FMT) (HI(I),RHOI(I),VPI(I),VSI(I),XII(I),
     *PHII(I),ETAI(I),qbetai(i),I=1,ILAY)
  170 CONTINUE
      IF(ICONT.EQ.1) GO TO 200
  180 CONTINUE
      DO 190 I=1,ILAY
      J=ILAY+1-I
      K=2*J
      L=K-1
      HI(K)=HI(J)
      HI(L)=0.0
      RHOI(K)=RHOI(J)
      RHOI(L)=RHOI(J)
      VPI(K)=VPI(J)
      VPI(L)=VPI(J)
      VSI(K)=VSI(J)
      VSI(L)=VSI(J)
      XII(K)=XII(J)
      XII(L)=XII(J)
      PHII(K)=PHII(J)
      PHII(L)=PHII(J)
      ETAI(K)=ETAI(J)
      ETAI(L)=ETAI(J)
      qbetai(k)=qbetai(j)
      qbetai(l)=qbetai(j)
  190 CONTINUE
      ILAY=2*ILAY
  200 CONTINUE
      READ(12,210)  NDIV,DMAX,MPUNCH
  210 FORMAT(I5,F10.0,I2)
C
      CALL  DIVIDE
      CALL  GRAVTY
C
      IF(LAME)  260,240,260
  240 CONTINUE
      DO  250  I=1,NLAY
      EMU(I)=RHO(I)*EMU(I)*EMU(I)
      ELAMB(I)=RHO(I)*ELAMB(I)*ELAMB(I)-2.0D+0*EMU(I)
  250 CONTINUE
C
  260 CONTINUE
      LAME1=1
C     IF(JSET-7)  290,290,270
      GO TO 290
  270 CONTINUE
      WRITE(JSET,5)  NAME
    5 FORMAT (20A4)
      WRITE (JSET,10) NLAY
   10 FORMAT (I10)
      DO  280  I=1,NLAY
      WRITE(JSET,15)  H(I),RHO(I),ELAMB(I),EMU(I),GR(I),XI(I),PHI(I),ETA
     *(I),qbetai(i)
   15 FORMAT (9E13.5)
  280 CONTINUE
C
  290 CONTINUE
C
      IF(MPUNCH.EQ.0)  GO TO  1000
      WRITE(10,310)  NAME
  310 FORMAT(20A4)
      WRITE(10,312)  NLAY,RAD,TMAS ,GS,RHO0,ISO,LAME1
  312 FORMAT(I5,F10.5,E15.7,F10.5,F10.6,2I2)
      IF(ISO.NE.0)  GO TO  350
      WRITE(10,314)
  314 FORMAT(9H(4D20.12))
      WRITE(10,316)  (H(I),RHO(I),ELAMB(I),EMU(I),I=1,NLAY)
  316 FORMAT(4D20.12)
      GO TO  1000
C
  350 CONTINUE
      WRITE(10,354)
  354 FORMAT(17H(4D20.12/3D20.12))
      WRITE(10,356)
     1     (H(I),RHO(I),ELAMB(I),EMU(I),XI(I),PHI(I),ETA(I),
     2     qbeta(i),I=1,NLAY)
  356 FORMAT(4D20.12/4D20.12)
      GO TO  1000
C
C TAPE INPUT
C
  500 CONTINUE
      READ(ISET)  NAME
      READ(ISET)  NLAY,RAD,TMASS,GS,RHO0,ISO,LAME,NDIV,DMAX
      TMAS=TMASS
      DO  510  I=1,NLAY
      READ(ISET)  H(I),RHO(I),ELAMB(I),EMU(I),GR(I),XI(I),PHI(I),ETA(I)
  510 CONTINUE
C
      DD=H(1)
      HH=0.0
      DO  520  I=1,NLAY
      DD=DD+HH
      HH=H(I+1)
      D(I)=DD
  520 CONTINUE
      GO TO  260
C
C EXIT
C
 1000 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE WNFRQ(ARG1,ARG2)
      IMPLICIT REAL*8(A-H,O-Z)
C PURPOSE -   TO COMPUTE WAVE NUMBER(WN), PHASE VELOCITY(C),
C     FREQUENCY(FRQ) AND PERIOD(T) FROM ARG1 AND ARG2.
C     WN2=WN**2(FLAT MODEL),  WN2=WN*(WN+1.0) (SPHERICAL MODEL)
C     C2=C**2
C     FRQ2=FRQ**2
C
      COMMON/VALUE /MODE,IERROR,I,ISTEP,J,K,L,N,NSOL,ISUM,KMAX,
     1              NMAX,IBOTM,ITOP,LY,LD,ACR,ELLIP,
     2              RMAX,RC,HC,H1,H2,ALP,WN,WN2,C,C2,FRQ,FRQ2,T,
     3              U,ENGY,DELTA,
     4              A(6,6),F(20),Y(6),repon,njob,nper,jyy,jdd,npar,
     5			ncouc,hdeb,hmin,hmax
      DATA PEI/6.283185307179586D+00/
      GO TO  (100,100,300,300,500,500),MODE
  100 CONTINUE
C
      T=ARG1
      C=ARG2
      FRQ=PEI/T
      WN=FRQ/C
      WN2=WN*WN
      GO TO  1000
C
C SPHERICAL MODEL
C ARG1=PERIOD,  ARG2=PHASE VELOCITY
C IF(ARG1.LT.0)  ARG1=-(WAVE NUMBER),  ARG2=PERIOD
C
  300 CONTINUE
      IF(ARG1)  500,1000,310
  310 CONTINUE
      T=ARG1
      C=ARG2
      FRQ=PEI/T
      WN    =(FRQ*RMAX)/C-0.5D+0
      WN2   =WN*(WN+1.0D+0)
      GO TO  1000
C
C FREE OSCILLATION
C ARG1=WAVE NUMBER,  ARG2=PERIOD
C
  500 CONTINUE
      WN=ABS(ARG1)
      T=DABS(ARG2)
      FRQ=PEI/T
      C     =(FRQ*RMAX)/(WN+0.5D+0)
      WN2   =WN*(WN+1.0D+0)
      GO TO  1000
C
C EXIT
C
 1000 CONTINUE
      C2=C*C
      FRQ2=FRQ*FRQ
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE  FUNCT(YC,XC,P,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/VALUE /MODE,IERROR,I,ISTEP,J,K,L,N,NSOL,ISUM,KMAX,
     1              NMAX,IBOTM,ITOP,LY,LD,ACR,ELLIP,
     2              RMAX,RC,HC,H1,H2,ALP,WN,WN2,C,C2,FRQ,FRQ2,T,
     3              U,ENGY,DELTA,
     4              A(6,6),F(20),Y(6),repon,njob,nper,jyy,jdd,npar,
     5			ncouc,hdeb,hmin,hmax
C
C A DUMMY SUBPROGRAM
C
      CALL  WNFRQ(P,XC)
      J=-1
      CALL  INTEG
      IER=IERROR
      YC=DELTA
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE  INTEG
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/MODEL /
     1              RAD,TMASS,GS,RHO0,GA,
     2              HI(500),RHOI(500),VPI(500),VSI(500),NAME(20),
     3              XII(500),PHII(500),ETAI(500),qbetai(500),
     4              D(1000),H(1000),RHO(1000),ELAMB(1000),EMU(1000),
     5              GR(1000),XI(1000),PHI(1000),ETA(1000),qbeta(1000),
     6              DMAX,ICONT,NLAY,ISO,LAME,NDIV,ILAY,idpmi,idpmax
      COMMON/VALUE /MODE,IERROR,I,ISTEP,J,K,L,N,NSOL,ISUM,KMAX,
     1              NMAX,IBOTM,ITOP,LY,LD,ACR,ELLIP,
     2              RMAX,RC,HC,H1,H2,ALP,WN,WN2,C,C2,FRQ,FRQ2,T,
     3              U,ENGY,DELTA,
     4              A(6,6),F(20),Y(6),repon,njob,nper,jyy,jdd,npar,
     5			ncouc,hdeb,hmin,hmax
      COMMON/SOL   /YN(6,1000),YB(6,3,20),SUM(20),Q(3,21),
     1              WNB(100),TT(100),CC(100),UU(100),ENG(100),ELL(100),
     2              AC(100)
      COMMON/Y123  /Y1(6,3),Y2(6,3),Y3(6,3)
C
C PURPOSE - TO INTEGRATE DIFFERENTIAL EQUATION UPWARD
C INPUT - J=-1 WHEN SEARCHING FOR AN EIGENVALUE
C          = 2  WHEN CALCULATING EIGENFUNCTIONS
C         ITOP,IBOTM= UPPER AND LOWER LIMIT OF I
C
C INITIALIZATION
C
      I=IBOTM
C
      CALL  INIVAL
C
      IF(IERROR)  900,1,900
    1 CONTINUE
C
      DO  2  NS=1,NSOL
      DO  2  LL=1,L
      Y1(LL,NS)=YB(LL,NS,1)
    2 CONTINUE
C
      IBOTM=I
      NN=0
      N=NMAX+1
      K=1
      ISTEP=-1
C
C
  100 CONTINUE
      NN=NN+1
      N=N-1
      IF(J-1)  200,200,110
  110 CONTINUE
C
      CALL  EIGFUN(Y1)
C
  200 CONTINUE
      IF(I-ITOP)  1000,1000,210
  210 CONTINUE
      IF(H(I))  500,300,500
C
C DISCONTINUITY
C
  300 CONTINUE
      I=I-1
      IF(EMU(I+1))  310,320,310
  310 CONTINUE
      IF(EMU(I))  100,330,100
  320 CONTINUE
      IF(EMU(I))  350,100,350
  330 CONTINUE
      GO TO  (340,350,340,350,340,350),MODE
C
C LOVE WAVE TYPE
C
  340 CONTINUE
      I=I+1
      ITOP=I
      GO TO  1000
C
C DISCONTINUITY OF DIFFERENTIAL EQUATION
C
  350 CONTINUE
      K=K+1
C
      CALL  CONT
C
      GO TO  100
C
C NORMAL SEQUENCE - TO INTEGRATE 2-STEPS
C
  500 CONTINUE
C
      CALL  STEP
C
      IF(J-1)  700,700,600
C
C ENERGY INTEGRAL
C
  600 CONTINUE
      N=N-1
      CALL  EIGFUN(Y2)
      N=N-1
      CALL  EIGFUN(Y3)
      N=N+2
      CALL  ENERGY
      N=N-2
C
  700 CONTINUE
      I=I-4
      NN=NN+2
C
      DO  710  NS=1,NSOL
      DO  710  LL=1,L
      Y1(LL,NS)=Y3(LL,NS)
  710 CONTINUE
      GO TO  200
C
C
 1000 CONTINUE
      NMAX=NN
      K=K+1
      KMAX=K
      IF(J-1)  1100,1100,1200
C
C STORE SURFACE VALUES
C
 1100 CONTINUE
      DO  1110  NS=1,NSOL
      DO  1110  LL=1,L
      YB(LL,NS,K)=Y1(LL,NS)
 1110 CONTINUE
C
      IF(J)  1120,900,900
C
C COMPUTE CHARACTERISTIC EQUATION
C
 1120 CONTINUE
C
      CALL  CHAREQ
      GO TO  900
C
C NORMALIZATION OF EIGENFUNCTION
C
 1200 CONTINUE
      W=1.0D+0
      IF(MODE.GE.3)  W=RAD-D(ITOP)
      W=W/YN(1,1)
      WW=W*W
      DO  1210  NNN=1,NMAX
      DO  1210  LL=1,6
      YN(LL,NNN)=W*YN(LL,NNN)
 1210 CONTINUE
      DO  1220  IS=1,ISUM
      SUM(IS)=WW*SUM(IS)
 1220 CONTINUE
C
C ENERGY INTEGRALS IN THE HOMOGENEOUS SPHERE
C
      CALL  DENERG
C
      GO TO  900
C
C EXIT
C
  900 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE  STEP
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/MODEL /
     1              RAD,TMASS,GS,RHO0,GA,
     2              HI(500),RHOI(500),VPI(500),VSI(500),NAME(20),
     3              XII(500),PHII(500),ETAI(500),qbetai(500),
     4              D(1000),H(1000),RHO(1000),ELAMB(1000),EMU(1000),
     5              GR(1000),XI(1000),PHI(1000),ETA(1000),qbeta(1000),
     6              DMAX,ICONT,NLAY,ISO,LAME,NDIV,ILAY,idpmi,idpmax
      COMMON/VALUE /MODE,IERROR,I,ISTEP,J,K,L,N,NSOL,ISUM,KMAX,
     1              NMAX,IBOTM,ITOP,LY,LD,ACR,ELLIP,
     2              RMAX,RC,HC,H1,H2,ALP,WN,WN2,C,C2,FRQ,FRQ2,T,
     3              U,ENGY,DELTA,
     4              A(6,6),F(20),Y(6),repon,njob,nper,jyy,jdd,npar,
     5			ncouc,hdeb,hmin,hmax
      COMMON/Y123  /Y1(6,3),Y2(6,3),Y3(6,3)
C
C PURPOSE - TO INTEGRATE 2-STEPS BY A RUNGE-KUTTA METHOD
C INITIAL VALUE  Y1(RC,I)
C 1-ST STEP      Y2(RC+HC,I-2)
C S-ND STEP      Y3(RC+2*HC,I-4)
C IF(J.NE.0)  RC AND HC ARE COMPUTED FROM D(I)
C
C SUBROUTINE REQUIRED - RUNGE
C
      SAVE=RC
      III=I
      ALP=0.5D+0
C
C 1-ST STEP
C
      IF(J)  100,200,100
  100 CONTINUE
      RC=RAD-D(I)
      H1=H(I)
      H2=H(I-1)
      HC=H1+H2
      ALP=H1/HC
  200 CONTINUE
      CALL  RUNGE(Y1,Y2)
C
C 2-ND STEP
C
      RC=RC+HC
      I=I-2
      IF(J)  300,400,300
  300 CONTINUE
      RC=RAD-D(I)
      H1=H(I)
      H2=H(I-1)
      HC=H1+H2
      ALP=H1/HC
  400 CONTINUE
      CALL  RUNGE(Y2,Y3)
C
C EXIT
C
 1000 CONTINUE
      RC=SAVE
      I=III
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE  RUNGE(Y1,Y2)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/VALUE /MODE,IERROR,I,ISTEP,J,K,L,N,NSOL,ISUM,KMAX,
     1              NMAX,IBOTM,ITOP,LY,LD,ACR,ELLIP,
     2              RMAX,RC,HC,H1,H2,ALP,WN,WN2,C,C2,FRQ,FRQ2,T,
     3              U,ENGY,DELTA,
     4              A(6,6),F(20),Y(6),repon,njob,nper,jyy,jdd,npar,
     5			ncouc,hdeb,hmin,hmax
      DIMENSION  Y1(6,3),Y2(6,3),YY(6,3),DY(6)
C
C PURPOSE -   TO INTEGRATE DIFFERENTIAL EQUATION BY A RUNGE-KUTTA METHOD
C             INITIAL VALUE IS Y1(RC,I)
C             OUTPUT IS        Y2(RC+HC,I+2*ISTEP)
C FORMULA     WHEN ALP=1/2, THE SCHEME USED COINCIDES WITH THE FOURTH
C     ORDER RUNGE-KUTTA FORMULA.   OTHERWISE IT IS  A THIRD ORDER
C     FORMULA.
C             HC=H1+H2,  ALP=H1/HC  ASSUMED
C             ON EXIT BOTH RC AND I KEEP INPUT VALUES.
C SUBROUTINE REQUIRED - DIFCOE
C
      SAVE=RC
      II=I
      WT    =0.5D+0-1.0D+0/(6.0D+0*ALP)
C
      CALL  DIFCOE
C
      DO  100  NS=1,NSOL
      DO  100  LL=1,L
      W=0.0
      DO  120  LLL=1,L
      W=W+A(LL,LLL)*Y1(LLL,NS)
  120 CONTINUE
      W=HC*W
      Y2(LL,NS)=Y1(LL,NS)+WT*W
      YY(LL,NS)=Y1(LL,NS)+ALP*W
  100 CONTINUE
C
      RC=RC+H1
      I=I+ISTEP
      WT=WT/ALP
C
      CALL  DIFCOE
C
      DO  200  NS=1,NSOL
      DO  210  LL=1,L
      W=0.0
      DO  220  LLL=1,L
      W=W+A(LL,LLL)*YY(LLL,NS)
  220 CONTINUE
      W=HC*W
      DY(LL)=W
      Y2(LL,NS)=Y2(LL,NS)+WT*W
  210 CONTINUE
      DO  230  LLL=1,L
      YY(LLL,NS)=Y1(LLL,NS)+ALP*DY(LLL)
  230 CONTINUE
  200 CONTINUE
      WT    =1.0D+0/(6.0D+0*ALP*(1.0D+0-ALP))-WT
C
      DO  300  NS=1,NSOL
      DO  310  LL=1,L
      W=0.0
      DO  320  LLL=1,L
      W=W+A(LL,LLL)*YY(LLL,NS)
  320 CONTINUE
      W=HC*W
      DY(LL)=W
      Y2(LL,NS)=Y2(LL,NS)+WT*W
  310 CONTINUE
      DO  330  LLL=1,L
      YY(LLL,NS)=Y1(LLL,NS)+DY(LLL)
  330 CONTINUE
  300 CONTINUE
C
      RC=RC+H2
      I=I+ISTEP
      WT    =0.5D+0-1.0D+0/(6.0D+0*(1.0D+0-ALP))
C
      CALL  DIFCOE
C
      DO  400  NS=1,NSOL
      DO  400  LL=1,L
      W=0.0
      DO  420  LLL=1,L
      W=W+A(LL,LLL)*YY(LLL,NS)
  420 CONTINUE
      Y2(LL,NS)=Y2(LL,NS)+WT*HC*W
  400 CONTINUE
C
C EXIT
C
 1000 CONTINUE
      RC=SAVE
      I=II
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE  CHAREQ
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/VALUE /MODE,IERROR,I,ISTEP,J,K,L,N,NSOL,ISUM,KMAX,
     1              NMAX,IBOTM,ITOP,LY,LD,ACR,ELLIP,
     2              RMAX,RC,HC,H1,H2,ALP,WN,WN2,C,C2,FRQ,FRQ2,T,
     3              U,ENGY,DELTA,
     4              A(6,6),F(20),Y(6),repon,njob,nper,jyy,jdd,npar,
     5			ncouc,hdeb,hmin,hmax
      COMMON/SOL   /YN(6,1000),YB(6,3,20),SUM(20),Q(3,21),
     1              WNB(100),TT(100),CC(100),UU(100),ENG(100),ELL(100),
     2              AC(100)
C
C PURPOSE -   TO EVALUATE CHARACTERISTIC EQUATION
C INPUT -     K
C OUTPUT -    DELTA
C
      GO TO  (100,200,100,200,100,600),MODE
C
C LOVE WAVE TYPE
C
  100 CONTINUE
      DELTA=YB(2,1,K)
      GO TO  1000
C
C RAYLEIGH WAVE TYPE
C
  200 CONTINUE
      IF(NSOL-1)  100,100,210
C DETERMINANT OF A 2X2 MATRIX
  210 CONTINUE
      DELTA=YB(2,1,K)*YB(4,2,K)-YB(2,2,K)*YB(4,1,K)
      GO TO  1000
C
C SPHEROIDAL OSCILLATION
C
  600 CONTINUE
      IF(NSOL-2)  100,210,610
C DETERMINANT OF A 3X3 MATRIX
  610 CONTINUE
      DELTA=YB(2,1,K)*(YB(4,2,K)*YB(6,3,K)-YB(4,3,K)*YB(6,2,K))
     1     +YB(2,2,K)*(YB(4,3,K)*YB(6,1,K)-YB(4,1,K)*YB(6,3,K))
     2     +YB(2,3,K)*(YB(4,1,K)*YB(6,2,K)-YB(4,2,K)*YB(6,1,K))
      GO TO  1000
C
C EXIT
C
 1000 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE  EIGFUN(YY)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/VALUE /MODE,IERROR,I,ISTEP,J,K,L,N,NSOL,ISUM,KMAX,
     1              NMAX,IBOTM,ITOP,LY,LD,ACR,ELLIP,
     2              RMAX,RC,HC,H1,H2,ALP,WN,WN2,C,C2,FRQ,FRQ2,T,
     3              U,ENGY,DELTA,
     4              A(6,6),F(20),Y(6),repon,njob,nper,jyy,jdd,npar,
     5			ncouc,hdeb,hmin,hmax
      COMMON/SOL   /YN(6,1000),YB(6,3,20),SUM(20),Q(3,21),
     1              WNB(100),TT(100),CC(100),UU(100),ENG(100),ELL(100),
     2              AC(100)
      DIMENSION  YY(6,3)
C
C COMBINE NSOL INDEPENDENT SOLUTIONS
C
      DO  10  LL=1,6
      YN(LL,N)=0.0
   10 CONTINUE
C
C
      DO  200  LL=1,L
      WS=0.0
      DO  100  NS=1,NSOL
      WS=WS+YY(LL,NS)
  100 CONTINUE
      YN(LL,N)=WS
  200 CONTINUE
C
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE  ENERGY
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/MODEL /
     1              RAD,TMASS,GS,RHO0,GA,
     2              HI(500),RHOI(500),VPI(500),VSI(500),NAME(20),
     3              XII(500),PHII(500),ETAI(500),qbetai(500),
     4              D(1000),H(1000),RHO(1000),ELAMB(1000),EMU(1000),
     5              GR(1000),XI(1000),PHI(1000),ETA(1000),qbeta(1000),
     6              DMAX,ICONT,NLAY,ISO,LAME,NDIV,ILAY,idpmi,idpmax
      COMMON/VALUE /MODE,IERROR,I,ISTEP,J,K,L,N,NSOL,ISUM,KMAX,
     1              NMAX,IBOTM,ITOP,LY,LD,ACR,ELLIP,
     2              RMAX,RC,HC,H1,H2,ALP,WN,WN2,C,C2,FRQ,FRQ2,T,
     3              U,ENGY,DELTA,
     4              A(6,6),F(20),Y(6),repon,njob,nper,jyy,jdd,npar,
     5			ncouc,hdeb,hmin,hmax
      COMMON/SOL   /YN(6,1000),YB(6,3,20),SUM(20),Q(3,21),
     1              WNB(100),TT(100),CC(100),UU(100),ENG(100),ELL(100),
     2              AC(100)
C
C PURPOSE -   TO CALCULATE ENERGY INTEGRALS USING THREE ORDINATES,
C     (I,N), (I-2,N-1), AND (I-4,N-2).
C             ON EXIT INPUT VALUES ARE RETURNED TO RC AND I.
C SUBROUTINE REQUIRED - INTGND
C
      SAVE=RC
      II=I
      H1=H(I)+H(I-1)
      H2=H(I-2)+H(I-3)
      W1=(H1+H2)/6.0D+0
      ALP=H2/H1
C
      W2    =W1*(2.0D+0-ALP)
      RC=RAD-D(I)
      DO  10  LL=1,L
      Y(LL)=YN(LL,N)
   10 CONTINUE
C
      CALL  INTGND
C
      DO  100  IS=1,ISUM
      SUM(IS)=SUM(IS)+W2*F(IS)
  100 CONTINUE
C
      W2    =W1*(2.0D+0+ALP+1.0D+0/ALP)
      I=I-2
      RC=RAD-D(I)
      DO  20  LL=1,L
      Y(LL)=YN(LL,N-1)
   20 CONTINUE
C
      CALL  INTGND
C
      DO  200  IS=1,ISUM
      SUM(IS)=SUM(IS)+W2*F(IS)
  200 CONTINUE
C
      W2    =W1*(2.0D+0-1.0D+0/ALP)
      I=I-2
      RC=RAD-D(I)
      DO  30  LL=1,L
      Y(LL)=YN(LL,N-2)
   30 CONTINUE
C
      CALL  INTGND
C
      DO  300  IS=1,ISUM
      SUM(IS)=SUM(IS)+W2*F(IS)
  300 CONTINUE
C
C EXIT
C
 1000 CONTINUE
      RC=SAVE
      I=II
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE  DERIV
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/MODEL /
     1              RAD,TMASS,GS,RHO0,GA,
     2              HI(500),RHOI(500),VPI(500),VSI(500),NAME(20),
     3              XII(500),PHII(500),ETAI(500),qbetai(500),
     4              D(1000),H(1000),RHO(1000),ELAMB(1000),EMU(1000),
     5              GR(1000),XI(1000),PHI(1000),ETA(1000),qbeta(1000),
     6              DMAX,ICONT,NLAY,ISO,LAME,NDIV,ILAY,idpmi,idpmax
      COMMON/VALUE /MODE,IERROR,I,ISTEP,J,K,L,N,NSOL,ISUM,KMAX,
     1              NMAX,IBOTM,ITOP,LY,LD,ACR,ELLIP,
     2              RMAX,RC,HC,H1,H2,ALP,WN,WN2,C,C2,FRQ,FRQ2,T,
     3              U,ENGY,DELTA,
     4              A(6,6),F(20),Y(6),repon,njob,nper,jyy,jdd,npar,
     5			ncouc,hdeb,hmin,hmax
C
C PURPOSE - TO CALCULATE DERIVATIVES OF EIGENVALUE
C SUBROUTINE REQUIRED - DGINTG, DCDP
C
      DG3=0.0
      N=0
      I=ITOP-1
C
  100 CONTINUE
      N=N+1
      I=I+1
      RC=RAD-D(I)
      IF(MODE.EQ.6)  F3=DGINTG(RC,I,N)
      CALL  DCDP(DG3)
      IF(N.GE.NMAX)  GO TO  1000
  200 CONTINUE
      IF(MODE.NE.6)  GO TO  250
C
      H1=H(I+4)+H(I+3)
      H2=H(I+2)+H(I+1)
      F2=DGINTG(RC-H2,I+2,N+1)
      F1=DGINTG(RC-H1-H2,I+4,N+2)
      ALP=H2/H1
      W1    =-(ALP*ALP)/(1.0D+0+ALP)
      W2    =3.0D+0+ALP
      W3    =2.0D+0+1.0D+0/(1.0D+0+ALP)
      DG2   =DG3+(H2/6.0D+0)*(W1*F1+W2*F2+W3*F3)
      W1    =2.0D+0-ALP
      W2    =2.0D+0+ALP+1.0D+0/ALP
      W3    =2.0D+0-1.0D+0/ALP
      DG1   =DG3+((H1+H2)/6.0D+0)*(W1*F1+W2*F2+W3*F3)
C
  250 CONTINUE
      I=I+2
      N=N+1
      RC=RAD-D(I)
      CALL  DCDP(DG2)
      I=I+2
      N=N+1
      RC=RAD-D(I)
      CALL  DCDP(DG1)
C
  300 CONTINUE
      F3=F1
      DG3=DG1
      IF(N.GE.NMAX)  GO TO  1000
      IF(H(I+1))  200,100,200
C
C EXIT
C
 1000 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      FUNCTION  DGINTG(RR,II,NN)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/MODEL /
     1              RAD,TMASS,GS,RHO0,GA,
     2              HI(500),RHOI(500),VPI(500),VSI(500),NAME(20),
     3              XII(500),PHII(500),ETAI(500),qbetai(500),
     4              D(1000),H(1000),RHO(1000),ELAMB(1000),EMU(1000),
     5              GR(1000),XI(1000),PHI(1000),ETA(1000),qbeta(1000),
     6              DMAX,ICONT,NLAY,ISO,LAME,NDIV,ILAY,idpmi,idpmax
      COMMON/VALUE /MODE,IERROR,I,ISTEP,J,K,L,N,NSOL,ISUM,KMAX,
     1              NMAX,IBOTM,ITOP,LY,LD,ACR,ELLIP,
     2              RMAX,RC,HC,H1,H2,ALP,WN,WN2,C,C2,FRQ,FRQ2,T,
     3              U,ENGY,DELTA,
     4              A(6,6),F(20),Y(6),repon,njob,nper,jyy,jdd,npar,
     5			ncouc,hdeb,hmin,hmax
      COMMON/SOL   /YN(6,1000),YB(6,3,20),SUM(20),Q(3,21),
     1              WNB(100),TT(100),CC(100),UU(100),ENG(100),ELL(100),
     2              AC(100)
C
C PURPOSE - TO CALCULATE THE INTEGRAND OF DG/DRHO
C
      RH=RHO(II)
      EL=ELAMB(II)
      PH    =1.0D+0+PHI(II)
      ET    =1.0D+0+ETA(II)
      EP=ET/PH
      GC=GR(II)
      R2=RR*RR
      Y1=YN(1,NN)
      Y3=YN(3,NN)
      IF(EMU(II).NE.0.0)  GO TO  100
C LIQUID LAYER
      Y3    =((RH*GC/RR-2.0D+0*EL*(1.0D+0-ET*EP)/R2)*Y1
     1      -EP*YN(2,NN)/RR-RH*Y3/RR)
     2      /(FRQ2*RH-WN2*EL*(1.0D+0-ET*EP)/R2)
C
  100 CONTINUE
      DGINTG=(GA*RH*(2.0D+0*Y1-WN2*Y3)*Y1)/(RR*R2)
C
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE  DIVIDE
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/MODEL /
     1              RAD,TMASS,GS,RHO0,GA,
     2              HI(500),RHOI(500),VPI(500),VSI(500),NAME(20),
     3              XII(500),PHII(500),ETAI(500),qbetai(500),
     4              D(1000),H(1000),RHO(1000),ELAMB(1000),EMU(1000),
     5              GR(1000),XI(1000),PHI(1000),ETA(1000),qbeta(1000),
     6              DMAX,ICONT,NLAY,ISO,LAME,NDIV,ILAY,idpmi,idpmax
C
C PURPOSE - TO DIVIDE ONE STEP INTO NDIV STEPS
C           LAYERS BELOW DMAX REMAIN INTACT
C
      N=NDIV-1
      DD=0.0
      I=0
      J=0
      IF(NDIV-1)  500,500,10
   10 CONTINUE
      DD=HI(1)
      HH=DD
      DX    =1.0D+0/DFLOAT(NDIV)
C
  100 CONTINUE
      I=I+1
      J=J+1
      D(J)=DD
      H(J)=HH
      RHO(J)=RHOI(I)
      ELAMB(J)=VPI(I)
      EMU(J)=VSI(I)
      XI(J)=XII(I)
      PHI(J)=PHII(I)
      ETA(J)=ETAI(I)
      qbeta(j)=qbetai(i)
      IF(I-ILAY)  110,1000,1000
  110 CONTINUE
      IF(DD-DMAX)  120,500,500
C
C BACKWARD INTERPOLATION
C
  120 CONTINUE
      D1=HI(I+1)
      IF(I+1-ILAY)  130,400,400
  130 CONTINUE
      D2=HI(I+2)
      IF(D2)  140,400,140
  140 CONTINUE
C
      Q1=D1/(D1+D2)
      Q2=-D1/D2
      Q3=-Q1*Q2
      Q4    =-1.0D+0-Q1
      Q5    =1.0D+0-Q2
      Q6=-Q3
C
      AR=Q1*RHOI(I)+Q2*RHOI(I+1)+Q3*RHOI(I+2)
      BR=Q4*RHOI(I)+Q5*RHOI(I+1)+Q6*RHOI(I+2)
      CR=RHOI(I)
      AL=Q1*VPI(I)+Q2*VPI(I+1)+Q3*VPI(I+2)
      BL=Q4*VPI(I)+Q5*VPI(I+1)+Q6*VPI(I+2)
      CL=VPI(I)
      AM=Q1*VSI(I)+Q2*VSI(I+1)+Q3*VSI(I+2)
      BM=Q4*VSI(I)+Q5*VSI(I+1)+Q6*VSI(I+2)
      CM=VSI(I)
      AX=Q1*XII(I)+Q2*XII(I+1)+Q3*XII(I+2)
      BX=Q4*XII(I)+Q5*XII(I+1)+Q6*XII(I+2)
      CX=XII(I)
      AP=Q1*PHII(I)+Q2*PHII(I+1)+Q3*PHII(I+2)
      BP=Q4*PHII(I)+Q5*PHII(I+1)+Q6*PHII(I+2)
      CP=PHII(I)
      AE=Q1*ETAI(I)+Q2*ETAI(I+1)+Q3*ETAI(I+2)
      BE=Q4*ETAI(I)+Q5*ETAI(I+1)+Q6*ETAI(I+2)
      CE=ETAI(I)
      aq=Q1*qbetai(I)+Q2*qbetai(I+1)+Q3*qbetai(I+2)
      bq=Q4*qbetai(I)+Q5*qbetai(I+1)+Q6*qbetai(I+2)
      cq=qbetai(I)
      D3=D1
C
  200 CONTINUE
      X=DX
      HH=D3*DX
C
      DO  300  NN=1,N
      J=J+1
      DD=DD+HH
      D(J)=DD
      H(J)=HH
      RHO(J)  =X*(AR*X+BR)+CR
      ELAMB(J)=X*(AL*X+BL)+CL
      EMU(J)  =X*(AM*X+BM)+CM
      XI(J)   =X*(AX*X+BX)+CX
      PHI(J)  =X*(AP*X+BP)+CP
      ETA(J)  =X*(AE*X+BE)+CE
      qbeta(j)=x*(aq*x+bq)+cq
      X=X+DX
  300 CONTINUE
      I=I+1
      J=J+1
      DD=DD+HH
      D(J)=DD
      H(J)=HH
      RHO(J)=RHOI(I)
      ELAMB(J)=VPI(I)
      EMU(J)=VSI(I)
      XI(J)=XII(I)
      PHI(J)=PHII(I)
      ETA(J)=ETAI(I)
      qbeta(j)=qbetai(i)
C
      IF(I-ILAY)  310,1000,1000
  310 CONTINUE
      IF(DD-DMAX)  320,500,500
  320 CONTINUE
      HH=HI(I+1)
      IF(HH)  330,100,330
  330 CONTINUE
C
C FORWARD INTERPOLATION
C
      D1=HI(I)
      IF(D1)  340,400,340
  340 CONTINUE
      D2=HI(I+1)
C
      Q2=-D2/D1
      Q3=D2/(D1+D2)
      Q1=-Q2*Q3
      Q4=-Q1
      Q5    =-1.0D+0-Q2
      Q6    =1.0D+0-Q3
C
      AR=Q1*RHOI(I-1)+Q2*RHOI(I)+Q3*RHOI(I+1)
      BR=Q4*RHOI(I-1)+Q5*RHOI(I)+Q6*RHOI(I+1)
      CR=RHOI(I)
      AL=Q1*VPI(I-1)+Q2*VPI(I)+Q3*VPI(I+1)
      BL=Q4*VPI(I-1)+Q5*VPI(I)+Q6*VPI(I+1)
      CL=VPI(I)
      AM=Q1*VSI(I-1)+Q2*VSI(I)+Q3*VSI(I+1)
      BM=Q4*VSI(I-1)+Q5*VSI(I)+Q6*VSI(I+1)
      CM=VSI(I)
      AX=Q1*XII(I-1)+Q2*XII(I)+Q3*XII(I+1)
      BX=Q4*XII(I-1)+Q5*XII(I)+Q6*XII(I+1)
      CX=XII(I)
      AP=Q1*PHII(I-1)+Q2*PHII(I)+Q3*PHII(I+1)
      BP=Q4*PHII(I-1)+Q5*PHII(I)+Q6*PHII(I+1)
      CP=PHII(I)
      AE=Q1*ETAI(I-1)+Q2*ETAI(I)+Q3*ETAI(I+1)
      BE=Q4*ETAI(I-1)+Q5*ETAI(I)+Q6*ETAI(I+1)
      CE=ETAI(I)
      aq=Q1*qbetai(I-1)+Q2*qbetai(I)+Q3*qbetai(I+1)
      bq=Q4*qbetai(I-1)+Q5*qbetai(I)+Q6*qbetai(I+1)
      cq=qbetai(I)
      D3=D2
      GO TO  200
C
C LINEAR INTERPOLATION
C
  400 CONTINUE
      AR=0.0
      BR=RHOI(I+1)-RHOI(I)
      CR=RHOI(I)
      AL=0.0
      BL=VPI(I+1)-VPI(I)
      CL=VPI(I)
      AM=0.0
      BM=VSI(I+1)-VSI(I)
      CM=VSI(I)
      AX=0.0
      BX=XII(I+1)-XII(I)
      CX=XII(I)
      AP=0.0
      BP=PHII(I+1)-PHII(I)
      CP=PHII(I)
      AE=0.0
      BE=ETAI(I+1)-ETAI(I)
      CE=ETAI(I)
      aq=0.0
      bq=qbetai(I+1)-qbetai(I)
      cq=qbetai(I)
      D3=HI(I+1)
      GO TO  200
C
  500 CONTINUE
      I=I+1
      DO  600  II=I,ILAY
      J=J+1
      HH=HI(II)
      DD=DD+HH
      D(J)=DD
      H(J)=HH
      RHO(J)=RHOI(II)
      ELAMB(J)=VPI(II)
      EMU(J)=VSI(II)
      XI(J)=XII(II)
      PHI(J)=PHII(II)
      ETA(J)=ETAI(II)
      qbeta(j)=qbetai(ii)
  600 CONTINUE
C
C EXIT
C
 1000 CONTINUE
      NLAY=J
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE  GRAVTY
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/MODEL /
     1              RAD,TMASS,GS,RHO0,GA,
     2              HI(500),RHOI(500),VPI(500),VSI(500),NAME(20),
     3              XII(500),PHII(500),ETAI(500),qbetai(500),
     4              D(1000),H(1000),RHO(1000),ELAMB(1000),EMU(1000),
     5              GR(1000),XI(1000),PHI(1000),ETA(1000),qbeta(1000),
     6              DMAX,ICONT,NLAY,ISO,LAME,NDIV,ILAY,idpmi,idpmax
      DATA PEI/12.56637061435917D+00/
      DATA  G/6.67E-8/
C
C PURPOSE - TO COMPUTE THE GRAVITY DISTRIBUTION
C INPUT - TMASS=TOTAL MASS OF THE EARTH IN C.G.S.
C         GS   =SURFACE VALUE OF GRAVITY IN C.G.S.
C         RHO0 =DENSITY OF THE INNERMOST HOMOGENEOUS SPHERE
C IF(TMASS.EQ.0)  TMASS MAY BE CALCULATED FROM GS
C IF(GS   .EQ.0)  INNER SPHERE R.LT.R(NLAY) IS ASSUMED TO BE HOMOGENEOUS
C                 WITH DENSITY OF RHO0
C IF(RHO0 .EQ.0)  RHO(NLAY) IS ASSIGNED TO RHO0
C OUTPUT - GR  =GRAVITY IN (10*5) C.G.S.
C
      GA=PEI*G
      I=NLAY
      EM=0.0
      GR(I)=0.0
      F3=RHO(I)*((RAD-D(I))**2)
  100 CONTINUE
      IF(H(I))  200,110,200
  110 CONTINUE
      I=I-1
      GR(I)=EM
      F3=RHO(I)*((RAD-D(I))**2)
      IF(I-1)  500,500,200
  200 CONTINUE
      F1=F3
      H1=H(I)
      H2=H(I-1)
      A=H2/H1
      F2=RHO(I-1)*((RAD-D(I-1))**2)
      F3=RHO(I-2)*((RAD-D(I-2))**2)
      C1    =2.0D+0+A/(1.0D+0+A)
      C2    =3.0D+0+1.0D+0/A
      C3    =-1.0D+0/(A*(1.0D+0+A))
      GR(I-1)=GR(I)+(H1/6.0D+0)*(C1*F1+C2*F2+C3*F3)
      C1    =2.0D+0-A
      C2    =2.0D+0+A+1.0D+0/A
      C3    =2.0D+0-1.0D+0/A
      EM    =GR(I)+((H1+H2)/6.0D+0)*(C1*F1+C2*F2+C3*F3)
      GR(I-2)=EM
      I=I-2
      IF(I-1)  500,500,100
C
  500 CONTINUE
      IF(TMASS)  510,520,510
  510 CONTINUE
      EM    =(TMASS*1.0D-15)/PEI-GR(1)
      RHO0  =(3.0D+0*EM)/((RAD-D(NLAY))**3)
      GO TO  570
  520 CONTINUE
      IF(GS)  530,540,530
  530 CONTINUE
      TMASS =((PEI*GS*RAD*RAD)/GA)*1.0D+10
      GO TO  510
  540 CONTINUE
      IF(RHO0)  550,560,550
  550 CONTINUE
      EM=RHO0*((RAD-D(NLAY))**3)/3.0D+0
      TMASS =PEI*(EM+GR(1))*1.0D+15
      GO TO  570
  560 CONTINUE
      RHO0=RHO(NLAY)
      GO TO  550
  570 CONTINUE
      DO  600  I=1,NLAY
      GR(I)=(GA*(GR(I)+EM))/((RAD-D(I))**2)
  600 CONTINUE
      GS    =GR(1)*1.0D+5
C
C EXIT
C
 1000 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE  DLIST(JY,ISET,ep,aprm,incr)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/MODEL /
     1              RAD,TMASS,GS,RHO0,GA,
     2              HI(500),RHOI(500),VPI(500),VSI(500),NAME(20),
     3              XII(500),PHII(500),ETAI(500),qbetai(500),
     4              D(1000),H(1000),RHO(1000),ELAMB(1000),EMU(1000),
     5              GR(1000),XI(1000),PHI(1000),ETA(1000),qbeta(1000),
     6              DMAX,ICONT,NLAY,ISO,LAME,NDIV,ILAY,idpmi,idpmax
      COMMON/VALUE /MODE,IERROR,I,ISTEP,J,K,L,N,NSOL,ISUM,KMAX,
     1              NMAX,IBOTM,ITOP,LY,LD,ACR,ELLIP,
     2              RMAX,RC,HC,H1,H2,ALP,WN,WN2,C,C2,FRQ,FRQ2,T,
     3              U,ENGY,DELTA,
     4              A(6,6),F(20),Y(6),repon,njob,nper,jyy,jdd,npar,
     5			ncouc,hdeb,hmin,hmax
      COMMON/SOL/  YN(6,1000),YB(6,3,20),SUM(20),Q(3,21),
     1             WNB(100),TT(100),CC(100),UU(100),ENG(100),ELL(100),
     2             AC(100)
      DIMENSION  YY(7)
      DIMENSION ADP(100,7), EP(100), PARAM(7),APRM(100,7)
	character*80 titre
      save ipass
      data ipass/0/
      DATA IEP/0/
      DATA PI/3.1415926536D+00/

      if(ipass.ne.1) then
c     On ne passe ici qu'une fois par run ==> se limiter a un seul modele
c     ndisc = nb de discontinuites dans le modele
c  On a: nlay=ilay+(ndiv-1)*(ilay-1-ndisc) ==>
	 ndisc=ilay-1-(nlay-ilay)/(ndiv-1)
c     indice=(indice initial)*NDIV - (NDIV-1)*(1+nb de profondeurs doublees)
	idpmi=idpmi*ndiv-(ndiv-1)*(1+ndisc)
	idpmax=idpmax*ndiv-(ndiv-1)*(1+ndisc)
	
      npar=7
	ncouc=(idpmax-idpmi)/ndiv +2
      if(ncouc.gt.100) stop 'depassement dimension EP()'
      WRITE(titre,201)  NAME
      write(15,'(a)') titre(:lnblnk(titre))
      write(15,*) njob
      write(15,*) npar, ncouc
      write(15,*) nper
      ipass=1
      endif

c  On réinitialise les couches min-max et les der.part
c     a chaque passage dans cette routine
      IDPMIN=IDPMI
      INCR=1
	do 275 i=1,ncouc
      ep(i)=0.
	do 275 ll=1,7
	adp(i,ll)=0.
  275	continue
C
C PRINT AND/OR PUNCH DERIVATIVES
C
C TITLE
C
      SWN=WN
      ST=T
      SC=C
      SU=U
      SENG=ENGY
      SR=RAD
      SG=GS
      SM=TMASS
      SRH=RHO0
      MODE1=MODE+10
C     IF(ISET.LE.7)  GO TO  1
      GO TO 1
 6457 WRITE(ISET,5)  MODE1,SWN,ST,SC,SU,SENG,NMAX,IBOTM,ITOP
    5 FORMAT (I5,5E13.5,3I10)
    1 CONTINUE
      IF(JY-2)  100,200,100
C
C PRINT
C
  100 CONTINUE
C     WRITE(14,101)
  101 FORMAT(1X,4HK(N),5X,1HT,8X,1HC,6X,1HU,7X,3HK.E,4X,
     1 5HY3/Y1,3X,5HERROR)
C     WRITE(14,102)  WN,T,C,U,ENGY,ELLIP,ACR
  102 FORMAT(1X,F4.0,F11.5,2F7.4,E10.4,F6.4,E10.4,2H $)
C     IF(JY-2)  300,200,200
C
C PUNCH
C
  200 CONTINUE
C     WRITE(10,201)  NAME
  201 FORMAT(20A4)
C     WRITE(10,202)  NLAY,RAD,SM,GS,RHO0,ISO,LAME
  202 FORMAT(I5,F10.5,E15.7,F10.5,F10.6,2I2)
C     WRITE(10,203)  MODE1,SWN,ST,SC,SU,SENG
  203 FORMAT(I5,5E15.7)
C     WRITE(10,204)  NMAX,IBOTM,ITOP
  204 FORMAT(3I5)
C
  300 CONTINUE
      IF(JY-2)  310,320,310
  310 CONTINUE
      IF(MODE.GE.5)  GO TO  315
      WRITE(14,311)
  311 FORMAT(1H ,3X,5HDEPTH,8X,7HDC/DRHO,11X,5HDC/DP,12X,5HDC/DS,11X,
     1       6HDC/DXI,11X,7HDC/DPHI,10X,7HDC/DETA)
      GO TO  320
  315 CONTINUE
      WRITE(14,316)
  316 FORMAT(1H ,3X,5HDEPTH,8X,7HDT/DRHO,11X,5HDT/DP,12X,5HDT/DS,11X,
     1       6HDT/DXI,11X,7HDT/DPHI,10X,7HDT/DETA)
      IF(JY.GE.2) WRITE(10,'(41H ATTENTION: DERIVEES DE LA PERIODE PROPR
     *E)')
C
  320 CONTINUE
      N=0
      I=ITOP-1
      alpha=0.
C
  500 CONTINUE
      I=I+1
      n=n+1
  600 CONTINUE
      DO  601  LL=1,6
  601 YY(LL)=YN(LL,N)
      yy(7)=0.
      DD=D(I)
C     IF(ISET.GT.7)  WRITE(ISET,10) DD,(YY(LL),LL=1,6)
   10 FORMAT (7E13.5)
      IF(JY-2)  700,800,700
C
C PRINT
C
  700 CONTINUE
      WRITE(14,701)  DD,(YY(LL),LL=1,6)
  701 FORMAT(1H ,F10.5,6(3X,E14.7))
      IF(JY-2)  900,800,800
C
C PUNCH
C
  800 CONTINUE

C PUNCH DES DERIVEES PART. AU FORMAT UTILISABLE PAR dgk et estlin
c	 yy = der. part. relatives pour une couche de 1km
c	adp = der. part.  absolues pour la couche entiere
C   (X/C)*(DC/DX) A T=CSTE, POUR LA COUCHE ENTIERE
C   FONCTIONNE AVEC SAITO MANTEAU, QUI SORT BIEN DC/DX ET NON DT/DX
C  IDPMIN = INDICE DE LA COUCHE-MODELE (APRES DIVIDE) OU DOIT
C           COMMENCER LE CALCUL DES DERIVEES PART.
C  IDPMAX = INDICE DE LA COUCHE OU FINIT LE CALCUL
      PARAM(1)=RHO(I)
      PARAM(2)=SQRT((ELAMB(I)+2.*EMU(I))/RHO(I))
      PARAM(3)=SQRT(EMU(I)/RHO(I))
      PARAM(4)=1.+XI(I)
      PARAM(5)=1.+PHI(I)
      PARAM(6)=1.+ETA(I)

c     calcul du facteur de qualite a la periode T a
c     partir du facteur de qualite des ondes S
c     on calcule en fait le coefficient de l'exponentielle
c     alpha = pi/(U*Q(temporel)*T)= pi/(C*Q(spatial)*T)
c   Attention: a ce niveau, YY(2) et YY(3) contiennent les
c              derivees partielles sans dimension de VP et VS
c              pour une couche de 1 km.
c              qbeta contient l'inverse du facteur de qualite
c      On integre sur tous les points ou la d.p. yy existe
c         mais on n'ecrit les derivees partielles qu'aux
c         profondeurs du modele entre (ce qui est fait en 
c         modifiant idpmin a chaque passage).

      if(i.eq.1) then
      epaiss=d(i+1)
      else
          if(i.eq.nlay) then
          epaiss=d(i)/ndiv
          else
*JJ          epaiss=(D(I+2)-D(I-2))/2.
C l'indice I augmente de 2 en 2, cf etiq 900 + qques lignes
          epaiss=D(I+1)-D(I-1)
          endif
      endif

	if(mode.eq.3) then
c	onde de Love, on remplace Vsv par Vsh
      alpha=alpha+pi*epaiss*(4./3.*param(3)*param(3)*param(4)/
     *            param(2)/param(2)*yy(2)+2*yy(4))
     *            *qbeta(i)/(c*t)
	else
      alpha=alpha+pi*epaiss*(4./3.*param(3)*param(3)/
     *                            param(2)/param(2)*yy(2)+yy(3))
     *            *qbeta(i)/(c*t)
	endif

      IF(I.LT.IDPMIN.OR.I.GT.IDPMAX) GO TO 1502

      INCR=INCR+1
*JJ       write(67,*) 'ncouc',ncouc,' idpmin',idpmin,' i',i,
*JJ      & ' idpmax',idpmax,' incr',incr,' n=',n,' nmax=',nmax

C      INCR=1 DERIV. PART. BIDON CORRESPONDANT A LA PARTIE
C             SUPERFICIELLE DU MODELE, LA OU L'INVERSION
C             N'AGIT PAS. ON A BESOIN DE L'EPAISSEUR DE CETTE
C             COUCHE POUR CALER LA PROFONDEUR DANS INV
      IF(INCR.EQ.2) then
	EP(1)=D(I-NDIV)+(D(I)-D(I-NDIV))/2.
	hmin=d(i)
	hdeb=d(i-ndiv)
	endif

      EP(INCR)=(D(I+NDIV)-D(I-NDIV))/2.
	hmax=d(i)

      idpmin=idpmin+ndiv

      aprm(incr,7)=qbeta(i)

c     adp(incr,7) contient une valeur qui permet de calculer
c     dgk/dq par un algorithme semblable au dgk/dm

	if(mode.eq.3) then
c	onde de Love, on remplace Vsv par Vsh
      adp(incr,7)=0.5*c*ep(incr)*(4./3.*param(3)*param(3)*param(4)/
     *            param(2)/param(2)*yy(2)+2.*yy(4))
	else
      adp(incr,7)=0.5*c*ep(incr)*(4./3.*param(3)*param(3)/
     *                            param(2)/param(2)*yy(2)+yy(3))
	endif


      DO 1501 LL=1,6
      LLL=LL
      IF(LL.EQ.1) LLL=3
      IF(LL.EQ.2) LLL=1
      IF(LL.EQ.3) LLL=2
C     ORDRE DES DER. PART. : VS,RO,VP,XI,PHI,ETA,qbeta
C   ADP CONTIENT MAINTENANT LA DERIVEE PARTIELLE "ABSOLUE"
C       POUR UNE COUCHE ENTIERE
      ADP(INCR,LL)=YY(LLL)*EP(INCR)*C/PARAM(LLL)
      APRM(INCR,LL)=PARAM(LLL)
 1501 CONTINUE
 1502 CONTINUE
C     WRITE(15,801)  DD,(YY(LL),LL=1,6)
  801 FORMAT(5E15.7)
C

  900 CONTINUE
c     Cas ou le calcul s'arrete car on atteint la prof max d'integration
c      On sort par n.ge.nmax quand la profondeur de debut d'integration
c      est inferieure a la profondeur max demandee pour les der.part (idpmax)
c  NB: en sortant par n.ge.nmax, on risque d'avoir incr < ncouc
c      et donc des der.part non calculees entre incr+1 et ncouc.
c      il est donc primordial que le tableau adp() soit remis a zero
c      a chaque passage dans DLIST, car on ecrit adp de 1 a ncouc.
c     Avant 2002, cette initialisation n'etait faite que la 1ere fois (bug).
      IF(N.GE.NMAX)  GO TO  1000
c     Cas ou l'on arrive a une discontinuite du modele
      IF(H(I+1).EQ.0.0)  GO TO  500
c     Cas standard
      I=I+2
      n=n+1
      GO TO  600
C
C EXIT
C
 1000 CONTINUE
      IF(JY-2) 2000,1600,1600
 1600 continue
c	Cq= Vitesse de phase incluant l'effet de l'attenuation
c	(Kanamori & Anderson, 1977, Rev.Geophys.Sp.Phys.)
c	pour une periode de reference de 1 seconde
c	(en ondes de surface, on a generalement T > 1s , donc
c	la correction est negative et Cq < C ).
	Tr=1.
	Cq= C * ( 1.+ C*T*alpha*log(Tr/T)/(pi*pi) )
      WRITE(15,'(4X,3F10.4,e14.6,f10.4,2f6.0)')T,C,U,alpha,Cq,hmin,hmax
      DO 684 LL=1,7
      WRITE(15,'(5e16.8)') (ADP(I,LL),I=1,ncouc)
  684 CONTINUE
 2000 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE  TTL1(MODE)
      IMPLICIT REAL*8(A-H,O-Z)
C
      GO TO  (100,200,300,400,500,600),MODE
C
C      LOVE WAVE
C
  100 CONTINUE
      WRITE(14,101)
  101 FORMAT(1H ,22HLOVE WAVE (FLAT MODEL))
      GO TO  1000
C
C      RAYLEIGH WAVE
C
  200 CONTINUE
      WRITE(14,201)
  201 FORMAT(1H ,26HRAYLEIGH WAVE (FLAT MODEL))
      GO TO  1000
C
C      MANTLE LOVE WAVE
C
  300 CONTINUE
      WRITE(14,301)
  301 FORMAT(1H ,34HMANTLE LOVE WAVE (SPHERICAL MODEL))
      GO TO  1000
C
C      MANTLE RAYLEIGH WAVE
C
  400 CONTINUE
      WRITE(14,401)
  401 FORMAT(1H ,38HMANTLE RAYLEIGH WAVE (SPHERICAL MODEL))
      GO TO  1000
C
C      TORSIONAL OSCILLATION
C
  500 CONTINUE
      WRITE(14,501)
  501 FORMAT(1H ,21HTORSIONAL OSCILLATION)
      GO TO  1000
C
C      SPHEROIDAL OSCILLATION
C
  600 CONTINUE
      WRITE(14,601)
  601 FORMAT(1H ,22HSPHEROIDAL OSCILLATION)
      GO TO  1000
C
 1000 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE  TABLE(LAY,DEPTH,ISTEP,idisc)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/MODEL /
     1              RAD,TMASS,GS,RHO0,GA,
     2              HI(500),RHOI(500),VPI(500),VSI(500),NAME(20),
     3              XII(500),PHII(500),ETAI(500),qbetai(500),
     4              D(1000),H(1000),RHO(1000),ELAMB(1000),EMU(1000),
     5              GR(1000),XI(1000),PHI(1000),ETA(1000),qbeta(1000),
     6              DMAX,ICONT,NLAY,ISO,LAME,NDIV,ILAY,idpmi,idpmax
C
C PURPOSE -   TO FIND  LAY  SUCH THAT  ABS(DEPTH-D(LAY))=MIN.
C     IF THERE ARE TWO SOLUTIONS, LAY TAKES THE LARGER ONE.
C
ccc	modification jjl (20 April 1990): on prend la couche
ccc	correspondant a la profondeur demandee ou immediatement
ccc	plus grande

c	idisc = nbre de discontinuites entre la surface et DEPTH

      DIF=1.0E+20
      I=1
      J=0
	idisc=0
	ijj=0
  100 CONTINUE
      IF(H(I+1))  200,110,200
  110 CONTINUE
	if((ijj-2*(ijj/2)).ne.0) then
       write(*,*)
       write(*,*)
     & 'STOP: il faut un nb pair de pts-modele entre 2 discontinuites'
       write(13,*)
     & 'STOP: il faut un nb pair de pts-modele entre 2 discontinuites'
       write(14,*)
     & 'STOP: il faut un nb pair de pts-modele entre 2 discontinuites'
       write(15,*)
     & 'STOP: il faut un nb pair de pts-modele entre 2 discontinuites'
          stop
	endif
	idisc=idisc+1
	ijj=0
      I=I+1
  200 CONTINUE
	ijj=ijj+1
      DIFF=D(I)-DEPTH
ccc      IF(DIFF)  250,250,300
	if(diff) 250,310,310
  250 CONTINUE
      DIF=DIFF
      J=I
      I=I+ISTEP
      IF(I-NLAY)  100,200,1000
  300 CONTINUE
      IF(ABS(DIFF)-ABS(DIF))  310,310,1000
  310 CONTINUE
      J=I
C
C EXIT
C
 1000 CONTINUE
      LAY=J
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE  YLIST(JY,ISET,idisc)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/MODEL /
     1              RAD,TMASS,GS,RHO0,GA,
     2              HI(500),RHOI(500),VPI(500),VSI(500),NAME(20),
     3              XII(500),PHII(500),ETAI(500),qbetai(500),
     4              D(1000),H(1000),RHO(1000),ELAMB(1000),EMU(1000),
     5              GR(1000),XI(1000),PHI(1000),ETA(1000),qbeta(1000),
     6              DMAX,ICONT,NLAY,ISO,LAME,NDIV,ILAY,idpmi,idpmax
      COMMON/VALUE /MODE,IERROR,I,ISTEP,J,K,L,N,NSOL,ISUM,KMAX,
     1              NMAX,IBOTM,ITOP,LY,LD,ACR,ELLIP,
     2              RMAX,RC,HC,H1,H2,ALP,WN,WN2,C,C2,FRQ,FRQ2,T,
     3              U,ENGY,DELTA,
     4              A(6,6),F(20),Y(6),repon,njob,nper,jyy,jdd,npar,
     5			ncouc,hdeb,hmin,hmax
      COMMON/SOL   /YN(6,1000),YB(6,3,20),SUM(20),Q(3,21),
     1              WNB(100),TT(100),CC(100),UU(100),ENG(100),ELL(100),
     2              AC(100)
      DIMENSION  YY(6)
      data ipass/0/
      save ipass
c
c  Au premier passage, impression du modele sur l'unite 13
c
	 ndisc2=idisc
	 if((mode.eq.1.or.mode.eq.3.or.mode.eq.5).and.emu(2).eq.0.)
     &                      ndisc2=idisc-1
      if(ipass.eq.0)then
           call listm(ndisc2)
           ipass=1
      endif
C
C PRINT AND/OR PUNCH EIGENFUNCTIONS
C
C TITLE
C
      W=1.0
      SWN=WN
      ST=T
      SC=C
      SU=U
      srep=repon
      SENG=ENGY
      SR=RAD
      SG=GS
      SM=TMASS
      SRH=RHO0
C
C PRINT
C
  100 CONTINUE
      WRITE(14,101)
  101 FORMAT(1X,3X,1HT,8X,1HC,6X,1HU,2x,4hK(N),7X,3HK.E,4X,
     1 5HY3/Y1,3X,5HERROR)
      WRITE(14,102)  T,C,U,wn,ENGY,ELLIP,ACR
  102 FORMAT(1X,F8.3,1x,2F8.4,1x,f4.0,E10.4,F6.4,E10.4,2H $)
      IF(JY-2)  300,200,200
C
C PUNCH
C
  200 CONTINUE
      WRITE(13,201)  NAME
  201 FORMAT(20A4)
      WRITE(13,202)  NLAY,RAD,SM,GS,RHO0,ISO,LAME
  202 FORMAT(I5,F10.5,E15.7,F10.5,F10.6,2I2)
      WRITE(13,203)  MODE,SWN,ST,SC,SU,srep,SENG
  203 FORMAT(I5,6E15.7)
c-------------------- obsolete ci-dessous -------------------------
c     ndisc = nb de discontinuites dans le modele
***      ndisc=2*nmax-ibotm-1
c  On a: nlay=ilay+(ndiv-1)*(ilay-1-ndisc) ==>
	 ndisc=ilay-1-(nlay-ilay)/(ndiv-1)
c     ndep = nb de deplacements ecrits dans le fichier-deplacements
c	on ne garde qu'1 couche sur (ndiv/2) du modele apres ndiv
***      ndep=((ibotm-itop+1)+ndisc*(ndiv-1)-1)*2/ndiv - ndisc + 1
c-------------------- obsolete ci-dessus -------------------------
	 ndep=(ibotm-itop-ndisc2)*2/ndiv+ndisc2+1
      write(13,204)  ndep, ibotm, itop, ilay, nlay, ndisc, ndisc2,idisc
c     WRITE(13,204)  NMAX,IBOTM,ITOP
  204 FORMAT(16I5)
C
  300 CONTINUE
      IF(JY.EQ.0)  GO TO  320
      IF(JY-2)  310,320,310
  310 CONTINUE
      WRITE(14,311)
  311 FORMAT(1H ,3X,5HDEPTH,10X,2HY1,15X,2HY2,15X,2HY3,15X,2HY4,15X,
     1       2HY5,15X,2HY6)
  320 CONTINUE
C**** FOLLOWING CARDS HAVE BEEN INSERTED OR CHANGED BY BOB NORTH********
      DPREV=-100.
C**** END OF INSERTION/REPLACEMENT *************************************
      N=0
      I=ITOP-1
C
C
  500 CONTINUE
      I=I+1
      n=n+1
      GO TO  (510,520,510,520,510,560),MODE
C
C SET L VALUE
C
  510 CONTINUE
      L=2
      GO TO  600
  520 CONTINUE
      L=4
      IF(EMU(I).EQ.0.0)  L=2
      GO TO  600
  560 CONTINUE
      L=6
      IF(EMU(I).EQ.0.0)  L=4
      GO TO  600
C
  600 CONTINUE
      IF(MODE.GE.3)  W=1.0/(RAD-D(I))
      DO  601  LL=1,L
      YY(LL)=W*YN(LL,N)
  601 CONTINUE
      DD=D(I)
      IF(JY.EQ.0)  GO TO  900
      IF(JY-2)  700,800,700
C
C PRINT
C
  700 CONTINUE
      WRITE(14,701)  DD,(YY(LL),LL=1,L)
  701 FORMAT(1H ,F10.5,6(3X,E14.7))
      IF(JY-2)  900,800,800
C
C PUNCH
C
  800 CONTINUE
      WRITE(13,801)  L,DD,(YY(LL),LL=1,L)
  801 FORMAT(I5,(5E15.7))
C
  900 CONTINUE
      IF(N.GE.NMAX)  GO TO  1000
      IF(H(I+1).EQ.0.0)  GO TO  500
      I=I+ndiv/2
      n=n+ndiv/4
      GO TO  600
C
C EXIT
C
 1000 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE  INTERP(AX,AP,AC1,ADC,AC2,AEPSIL,JUMP)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  X(100),Y(100)
C     DIMENSION TTO(50)
C PURPOSE - TO FIND A ROOT OF Y(X,P)=0 BETWEEN X=C1(DC)C2
C           P IS A PARAMETER
C OUTPUT - XX=ROOT
C          JUMP=0  IF ROOT NOT FOUND
C              =1  IF ROOT FOUND
      P=AP
      C1=AC1
      C2=AC2
      DC=ADC
      EPSIL=AEPSIL
      M=0
      JUMP=0
C     WRITE(14,101)
  101 FORMAT(' INTERP VERSION ORIGINALE SAITO.GAULON')
      WRITE(14,103) YC,XC,P,IER
  103 FORMAT(1X,3E15.6,I10)
      XC=C1
  100 CONTINUE
      CALL  FUNCT(YC,XC,AP,IER)
C
      IF(IER)  1600,1050,1600
 1050 CONTINUE
C     WRITE(14,1051)  XC,YC
 1051 FORMAT(1H ,F12.5,5X,E12.5)
      IF(EPSIL.LT.1.E-10) GO TO 1400   
 1060 CONTINUE
      IF(M-1)  1100,1200,1310
C
C FIRST TRIAL
C
 1100 CONTINUE
      M=1
      X1=XC
      Y1=YC
      X(1)=XC
      Y(1)=YC
      IF(YC)  1400,1500,1400
C A ZERO CROSSING HAS NOT BEEN REACHED
 1200 CONTINUE
      IF(Y1)  1210,1500,1220
 1210 CONTINUE
      IF(YC)  1100,1500,1300
 1220 CONTINUE
      IF(YC)  1300,1500,1100
C
C A ZERO CROSSING IS FOUND
C
 1300 CONTINUE
      X2=XC
      Y2=YC
      GO TO  1350
C TO CHECK IF YC IS BETWEEN Y1 AND Y2
 1310 CONTINUE
      IF(Y1-YC)  1320,1320,1330
 1320 CONTINUE
      IF(Y2-YC)  1340,1340,1350
 1330 CONTINUE
      IF(Y2-YC)  1350,1340,1340
C YC IS OUTSIDE OF (Y1,Y2)
 1340 CONTINUE
      IF(YC)  1342,1342,1344
 1342 CONTINUE
      IF(Y1)  1343,1343,1345
 1343 CONTINUE
      X1=XC
      Y1=YC
      GO TO  1346
 1344 CONTINUE
      IF(Y1)  1345,1343,1343
 1345 CONTINUE
      X2=XC
      Y2=YC
 1346 CONTINUE
C     X(1)=X1
C     Y(1)=Y1
C     X(2)=X2
C     Y(2)=Y2
C     M=1
C     GO TO  1355
C INTERPOLATION BY A LAGRANGE FORMULA
C
 1350 CONTINUE
      X(M+1)=XC
      Y(M+1)=YC
 1355 CONTINUE
C     DO  1360  KK=1,M
C     X(M-KK+1)=(-Y(M-KK+1)*X(M-KK+2)+Y(M+1)*X(M-KK+1))
C    1         /(Y(M+1)-Y(M-KK+1))
 1360 CONTINUE
      IF(Y(M+1).EQ.0.OR.Y(M).EQ.0.) GOTO 1500
      IF((Y(M+1)/DABS(Y(M+1)))*(Y(M)/DABS(Y(M))))  20,1500,30
   30 X(M) = X(M-1)
      Y(M) = Y(M-1)
   20 X(M+2) = X(M) -(X(M+1) - X(M)) / (Y(M+1)-Y(M))*Y(M)
      ERROR=(XC-X(M+2))/XC
C     ERROR=(XC-X(1))/XC
      IF(DABS(ERROR)-EPSIL) 1500,1500,1370
 1370 CONTINUE
C     XC=X(1)
      XC =  X(M+2)
      M=M+1
      IF(M.GE.98) GO TO 2000
      GO TO  1410
C INCREASE XC BY DC
C
 1400 CONTINUE
      XC=XC+DC
C CHECK IF XC IS IN THE RANGE (C1,C2)
 1410 CONTINUE
      IF((XC-C1)*(XC-C2))  100,100,1600
C
C A ROOT IS FOUND
C
 1500 CONTINUE
      JUMP=1
      AX=XC
      GO TO  1000
C
C A ROOT IS NOT FOUND
C
 1600 CONTINUE
      WRITE(14,1601)
 1601 FORMAT(1H ,7HNO ROOT)
      AX=C2
      GO TO  1000
 2000 WRITE(14,52) M,XC
      AX=XC
   52 FORMAT('$ DEPASSEMENT TAILLE Y DANS INTERP, M=',I4,' XC=',F10.4)
 1000 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE  MLIST
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/MODEL /
     1              RAD,TMASS,GS,RHO0,GA,
     2              HI(500),RHOI(500),VPI(500),VSI(500),NAME(20),
     3              XII(500),PHII(500),ETAI(500),qbetai(500),
     4              D(1000),H(1000),RHO(1000),ELAMB(1000),EMU(1000),
     5              GR(1000),XI(1000),PHI(1000),ETA(1000),qbeta(1000),
     6              DMAX,ICONT,NLAY,ISO,LAME,NDIV,ILAY,idpmi,idpmax
C
C PURPOSE - TO PRINT OUT MODEL
C
      WRITE(14,1001)
 1001 FORMAT(1H ,3X,5HDEPTH,8X,3HRHO,8X,2HVP,9X,2HVS,9X,1HG,10X,2HXI,
     1       8X,3HPHI,7X,3HETA)
      DO  1500  I=1,NLAY
      DC=D(I)
      S=EMU(I)/RHO(I)
      P=ELAMB(I)/RHO(I)+2.0*S
      P=SQRT(P)
      S=SQRT(S)
      WRITE(14,1501) DC,RHO(I),P,S,GR(I),XI(I),PHI(I),ETA(I),qbeta(i)
1501     FORMAT(1H.,F10.5,3F11.5,5PF11.3,0PF11.5,3F10.5)
1500      CONTINUE
C
      RETURN
       END
C            ()manteau.f:1.1    2/28/86
 
C ----------------------------------------------------------------------
      SUBROUTINE INIVAL
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/MODEL /
     1              RAD,TMASS,GS,RHO0,GA,
     2              HI(500),RHOI(500),VPI(500),VSI(500),NAME(20),
     3              XII(500),PHII(500),ETAI(500),qbetai(500),
     4              D(1000),H(1000),RHO(1000),ELAMB(1000),EMU(1000),
     5              GR(1000),XI(1000),PHI(1000),ETA(1000),qbeta(1000),
     6              DMAX,ICONT,NLAY,ISO,LAME,NDIV,ILAY,idpmi,idpmax
      COMMON/VALUE /MODE,IERROR,I,ISTEP,J,K,L,N,NSOL,ISUM,KMAX,
     1              NMAX,IBOTM,ITOP,LY,LD,ACR,ELLIP,
     2              RMAX,RC,HC,H1,H2,ALP,WN,WN2,C,C2,FRQ,FRQ2,T,
     3              U,ENGY,DELTA,
     4              A(6,6),F(20),Y(6),repon,njob,nper,jyy,jdd,npar,
     5			ncouc,hdeb,hmin,hmax
      COMMON/SOL/  YN(6,1000),YB(6,3,20),SUM(20),Q(3,21),
     1             WNB(100),TT(100),CC(100),UU(100),ENG(100),ELL(100),
     2             AC(100)
C
C INITIAL VALUES
C MANTLE WAVE
C Y3+WNN*Y4=0
C SUBROUTINE REQUIRED - ZN
C
      RH=RHO(I)
      EL=ELAMB(I)
      EM=EMU(I)
      RVP   =EL+2.0D+0*EM
      P=RVP/RH
      S=EM/RH
      RB=RAD-D(IBOTM)
      RT=RAD-D(ITOP)
      IF(J-1)  30,20,10
   10 CONTINUE
      DO  11  IS=1,ISUM
      SUM(IS)=0.0
   11 CONTINUE
      GO TO  30
   20 CONTINUE
      DO  21  KK=1,KMAX
      DO  21  LL=1,6
      YB(LL,3,KK)=YB(LL,1,KK)
   21 CONTINUE
C
   30 CONTINUE
C
      GO TO  (900,900,300,400,300,900),MODE
C
C MANTLE LOVE WAVE (SPHERICAL MODEL)
C
  300 CONTINUE
      NSOL=1
      L=2
      ISUM=3
      IF(J)  310,350,350
C
  310 CONTINUE
      IF(EM.EQ.0.0)  GO TO  120
C
C SOLID LAYER
C
      YB(1,1,1)=RB
      XB=(FRQ*RB)/DSQRT(S)
      YB(2,1,1)=EM*((WN-1.0D+0)-ZN(XB,WN,20))
      GO TO  1000
C
C LIQUID LAYER
C FIND A SOLID LAYER
C
  120 CONTINUE
      IF(I.LT.2)  GO TO  995
      DO  122  II=2,I
      III=I-II+2
      IF(EMU(III-1).NE.0.0)  GO TO  123
  122 CONTINUE
      GO TO  995
  123 CONTINUE
      IF(H(III).NE.0.0)  GO TO  995
      I=III
C
      YB(1,1,1)=1.0D+0
      YB(2,1,1)=0.0
      GO TO  1000
C
C
  350 CONTINUE
      Q(1,KMAX-1)=RT/YB(1,1,KMAX)
      IF(KMAX.LE.2)  GO TO  152
      K1=KMAX-2
      DO  151  KK=1,K1
      KKK=KMAX-KK-1
      Q(1,KKK)=Q(1,KKK+1)
  151 CONTINUE
  152 CONTINUE
      ELLIP=0.0
      GO TO  290
C
C MANTLE RAYLEIGH WAVE (SPHERICAL MODEL)
C
  400 CONTINUE
      IF(J-1)  410,430,460
C
  410 CONTINUE
      IF(EM.EQ.0.0)  GO TO  420
C
C SOLID LAYER
C
      NSOL=1
      ISUM=3
      L=5
      XA=(FRQ*RB)/DSQRT(P)
      XBB=(FRQ2*RB*RB)/S
      XB=DSQRT(XBB)
      ZA=ZN(XA,WN,20)
      ZB=ZN(XB,WN,20)
      WN1=WN+1.0D+0
      ER=EM/RB
C
      YB(1,1,1)=0.5D+0*(WN1*ZA+WN*ZB-ZA*ZB)/WN1
      YB(2,1,1)=ER*ER*(XBB*(-0.5D+0*XBB+(WN-1.0D+0)*(2.0D+0*WN+1.0D+0))
     1         +2.0D+0*(-WN1*(WN2-2.0D+0)+XBB)*ZA
     2         +(-2.0D+0*WN*(WN2-2.0D+0)+XBB)*ZB
     3         +2.0D+0*(WN2-2.0D+0)*ZA*ZB)/WN1
      YB(3,1,1)=ER*WN*(-0.5D+0*XBB+WN1*ZA+WN*ZB-ZA*ZB)
      YB(4,1,1)=ER*(0.5D+0*WN*XBB-(0.5D+0*XBB+WN1)*ZA-WN*ZB+ZA*ZB)/WN1
      YB(5,1,1)=2.0D+0*ER*(WN1*(-0.25D+0*XBB+ZA)
     1         +(WN+0.25D+0*XBB)*ZB-ZA*ZB)/WN1
      GO TO  1000
C
C LIQUID LAYER
C
  420 CONTINUE
      NSOL=1
      L=2
      ISUM=3
      XA=(FRQ*RB)/DSQRT(P)
      YB(1,1,1)=WN-ZN(XA,WN,20)
      YB(2,1,1)=-FRQ2*RH*RB
      GO TO  1000
C
  430 CONTINUE
      IF(EM.EQ.0.0)  GO TO  420
C
C SOLID LAYER
C
      NSOL=2
      L=4
      ISUM=3
      XA=(FRQ*RB)/DSQRT(P)
      XBB=(FRQ2*RB*RB)/S
      XB=DSQRT(XBB)
      ZA=ZN(XA,WN,20)
      ZB=ZN(XB,WN,20)
      WN1=WN+1.0D+0
      YB(1,1,1)=WN-0.5D+0*ZA
      YB(2,1,1)=(EM/RB)*(2.0D+0*WN*(WN-1.0D+0)-0.5D+0*XBB+2.0D+0*ZA
     1         -WN*ZB)
      YB(3,1,1)=1.0D+0-(0.5D+0*ZB)/WN1
      YB(4,1,1)=(EM/RB)*(2.0D+0*(WN-1.0D+0)-(0.5D+0*XBB)/WN1
     1         -ZA+ZB/WN1)
      YB(1,2,1)=-0.5D+0*ZA
      YB(2,2,1)=(EM/RB)*(-0.5D+0*XBB+2.0D+0*ZA+WN*ZB)
      YB(3,2,1)=(0.5D+0*ZB)/WN1
      YB(4,2,1)=(EM/RB)*((0.5D+0*XBB)/WN1-ZA-ZB/WN1)
      GO TO  1000
C
  460 CONTINUE
      IF(EMU(ITOP))  461,465,461
C
C SOLID SURFACE LAYER
C
  461 CONTINUE
      NSOL=2
      WS=YB(3,3,KMAX)**2+WN2*(YB(4,3,KMAX)**2)
      Q(1,KMAX-1)= RT*(YB(2,2,KMAX)*YB(3,3,KMAX)
     1           +WN2*YB(4,2,KMAX)*YB(4,3,KMAX))/WS
      Q(2,KMAX-1)=-RT*(YB(2,1,KMAX)*YB(3,3,KMAX)
     1           +WN2*YB(4,1,KMAX)*YB(4,3,KMAX))/WS
      ELLIP=YB(3,3,KMAX)*(YB(5,3,KMAX)-YB(4,3,KMAX))/WS
      GO TO  270
C
C LIQUID SURFACE LAYER
C
  465 CONTINUE
      NSOL=1
      Q(1,KMAX-1)=RT/YB(1,3,KMAX)
      ELLIP=-YB(3,3,KMAX-1)/(WN2*YB(4,3,KMAX-1))
C
  270 CONTINUE
      IF(KMAX.LE.2)  GO TO  279
      K1=KMAX-2
      DO  271  KK=1,K1
      KKK=KMAX-KK-1
      IF(NSOL-2)  272,273,273
  272 CONTINUE
      NSOL=2
      Q(1,KKK)= YB(4,2,KKK+1)*Q(1,KKK+1)
      Q(2,KKK)=-YB(4,1,KKK+1)*Q(1,KKK+1)
      GO TO  271
  273 CONTINUE
      NSOL=1
      Q(1,KKK)=Q(1,KKK+1)
  271 CONTINUE
C
  279 CONTINUE
      L=2*NSOL
C
C
  290 CONTINUE
      DO  275  NS=1,NSOL
      DO  275  LL=1,L
      YB(LL,NS,1)=Q(NS,1)*YB(LL,NS,1)
  275 CONTINUE
      GO TO  1000
C
C ERROR
C
  900 CONTINUE
      WRITE(14,901)  MODE
  901 FORMAT(1H0,5HMODE=,I2,3X,19HINIVAL(MANTLE WAVE),/,/)
      IERROR=3
      GO TO  1000
  995 CONTINUE
      WRITE(14,996)  NAME
  996 FORMAT(1H0,5HMODEL,2X,20A4//)
      IF(IERROR.LE.1)  IERROR=2
      GO TO  1000
C
C EXIT
C
 1000 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      FUNCTION  ZN(X,FN,M)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL  FIRST
C
C PURPOSE - TO COMPUTE  X*J   (X)/J (X)
C                          N+1     N
C
      FIRST=.TRUE.
      M1=M
    1 CONTINUE
      FM=M1
      XX=X*X
      C     =2.0D+0*(FN+FM)+1.0D+0
      CC    =C+2.0D+0*FM
      Z=0.0
      ZN=0.0
C
      DO  100  MM=1,M1
      Z=XX/(C-Z)
      ZN    =XX/(CC-2.0D+0-(XX/(CC-ZN)))
      C=C-2.0D+0
      CC=CC-4.0D+0
  100 CONTINUE
C
      E=(ZN-Z)/ZN
      IF(ABS(E).LE.1.0E-05)  GO TO  1000
      IF(.NOT.FIRST)  GO TO  900
      FIRST=.FALSE.
      M1=2*M
      GO TO  1
  900 CONTINUE
      WRITE(14,999)  ZN,Z,X,FN,M1
  999 FORMAT(1H ,30X,13HWARNING IN ZN,3X,3HZN=,D13.5,
     1       3X,2HZ=,D13.5,3X,2HX=,D13.5,3X,2HN=,D13.5,3X,2HM=,I3)
C
 1000 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE  CONT
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/MODEL /
     1              RAD,TMASS,GS,RHO0,GA,
     2              HI(500),RHOI(500),VPI(500),VSI(500),NAME(20),
     3              XII(500),PHII(500),ETAI(500),qbetai(500),
     4              D(1000),H(1000),RHO(1000),ELAMB(1000),EMU(1000),
     5              GR(1000),XI(1000),PHI(1000),ETA(1000),qbeta(1000),
     6              DMAX,ICONT,NLAY,ISO,LAME,NDIV,ILAY,idpmi,idpmax
      COMMON/VALUE /MODE,IERROR,I,ISTEP,J,K,L,N,NSOL,ISUM,KMAX,
     1              NMAX,IBOTM,ITOP,LY,LD,ACR,ELLIP,
     2              RMAX,RC,HC,H1,H2,ALP,WN,WN2,C,C2,FRQ,FRQ2,T,
     3              U,ENGY,DELTA,
     4              A(6,6),F(20),Y(6),repon,njob,nper,jyy,jdd,npar,
     5			ncouc,hdeb,hmin,hmax
      COMMON/SOL/  YN(6,1000),YB(6,3,20),SUM(20),Q(3,21),
     1             WNB(100),TT(100),CC(100),UU(100),ENG(100),ELL(100),
     2             AC(100)
      COMMON/Y123 /Y1(6,3),Y2(6,3),Y3(6,3)
C
C CONTINUITY AT SOLID-LIQUID BOUNDARIES
C SURFACE WAVE
C Y3+WNN*Y4=0
C
      DO  10  NS=1,NSOL
      DO  10  LL=1,L
      YB(LL,NS,K)=Y1(LL,NS)
   10 CONTINUE
C
      GO TO  (100,200,100,200,100,900),MODE
C
C LOVE WAVE TYPE
C
  100 CONTINUE
      IF(EMU(I))  1000,900,1000
C
C RAYLEIGH WAVE TYPE
C
  200 CONTINUE
      IF(J-1)  210,220,230
  210 CONTINUE
      IF(EMU(I))  211,212,211
C
C LIQUID TO SOLID
C
  211 CONTINUE
      NSOL=1
      L=5
      Y1(5,1)=-Y1(2,1)
      Y1(2,1)=0.0
      Y1(3,1)=0.0
      Y1(4,1)=0.0
      GO TO  1000
C
C SOLID TO LIQUID
C
  212 CONTINUE
      NSOL=1
      L=2
      Y1(1,1)=Y1(4,1)
      GO TO  1000
C
C
  220 CONTINUE
      IF(EMU(I))  221,222,221
C
C LIQUID TO SOLID
C
  221 CONTINUE
      NSOL=2
      L=4
      Y1(3,1)=0.0
      Y1(4,1)=0.0
      Y1(1,2)=0.0
      Y1(2,2)=0.0
      Y1(3,2)=1.0D+0
      Y1(4,2)=0.0
      GO TO  1000
C
C SOLID TO LIQUID
C
  222 CONTINUE
      NSOL=1
      L=2
      Y1(1,1)=YB(4,3,K)
      Y1(2,1)=YB(2,3,K)
      GO TO  1000
C
  230 CONTINUE
      IF(EMU(I))  232,231,232
C
C SOLID TO LIQUID
C
  231 CONTINUE
      NSOL=1
      L=2
      Y1(1,1)=Y1(1,1)+Y1(1,2)
      Y1(2,1)=Y1(2,1)+Y1(2,2)
      GO TO  1000
C
C LIQUID TO SOLID
C
  232 CONTINUE
      NSOL=2
      L=4
      Y1(3,1)=0.0
      Y1(4,1)=0.0
      Y1(1,2)=0.0
      Y1(2,2)=0.0
      Y1(3,2)=Q(2,K)
      Y1(4,2)=0.0
      GO TO  1000
C
C ERROR
C
  900 CONTINUE
      WRITE(14,901)  MODE
  901 FORMAT(1H0,5HMODE=,I2,3X,18HCONT(SURFACE WAVE),/,/)
      IERROR=3
      GO TO  1000
C
C EXIT
C
 1000 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE  DIFCOE
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/MODEL /
     1              RAD,TMASS,GS,RHO0,GA,
     2              HI(500),RHOI(500),VPI(500),VSI(500),NAME(20),
     3              XII(500),PHII(500),ETAI(500),qbetai(500),
     4              D(1000),H(1000),RHO(1000),ELAMB(1000),EMU(1000),
     5              GR(1000),XI(1000),PHI(1000),ETA(1000),qbeta(1000),
     6              DMAX,ICONT,NLAY,ISO,LAME,NDIV,ILAY,idpmi,idpmax
      COMMON/VALUE /MODE,IERROR,I,ISTEP,J,K,L,N,NSOL,ISUM,KMAX,
     1              NMAX,IBOTM,ITOP,LY,LD,ACR,ELLIP,
     2              RMAX,RC,HC,H1,H2,ALP,WN,WN2,C,C2,FRQ,FRQ2,T,
     3              U,ENGY,DELTA,
     4              A(6,6),F(20),Y(6),repon,njob,nper,jyy,jdd,npar,
     5			ncouc,hdeb,hmin,hmax
C
C COEFFICIENT CATRIX FOR DIFFERENTIAL EQUATIONS
C MANTLE WAVES
C ANISOTROPIC MODEL
C Y3+WNN*Y4=0
C
      DO  1  LL=1,6
      DO  1  LLL=1,6
      A(LL,LLL)=0.0
    1 CONTINUE
C
      RH=RHO(I)
      EL=ELAMB(I)
      EM=EMU(I)
      XC    =1.0D+0+XI(I)
      RVP   =EL+2.0D+0*EM
      PH    =1.0D+0+PHI(I)
      ET    =1.0D+0+ETA(I)
      EP=ET/PH
      R1    =1.0D+0/RC
      R2=R1*R1
C
      GO TO  (1000,1000,300,400,300,1000),MODE
C
C MANTLE LOVE WAVE
C TOROIDAL OSCILLATIONS
C
  300 CONTINUE
      IF(EM)  310,1000,310
  310 CONTINUE
      A(1,1)=2.0D+0*R1
      A(1,2)=1.0D+0/EM
      A(2,1)=-FRQ2*RH+R2*(WN2-2.0D+0)*XC*EM
      A(2,2)=-A(1,1)
      GO TO  1000
C
C MANTLE RAYLEIGH WAVES
C
  400 CONTINUE
      IF(EM)  410,450,410
  410 CONTINUE
      IF(J)  420,430,430
C
C DIFFERENTIAL EQUATION FOR CHARACTERISTIC EQUATION
C
  420 CONTINUE
      A(1,1)=R1*(3.0D+0-(2.0D+0*EP*EL)/RVP)
      A(1,4)=1.0D+0/EM
      A(1,5)=-1.0D+0/(PH*RVP)
      A(2,2)=-A(1,1)
      A(2,3)=4.0D+0*R2*(RVP-XC*EM-(ET*EP*EL*EL)/RVP)
      A(2,4)=-FRQ2*RH+A(2,3)
      A(2,5)=FRQ2*RH-R2*(WN2*(RVP-ET*EP*EL*EL/RVP)-2.0D+0*XC*EM)
      A(4,3)=-2.0D+0*R1*EP*EL/RVP
      A(4,4)=-R1*(1.0D+0+2.0D+0*EP*EL/RVP)
      A(5,3)=-2.0D+0*R1
      A(3,1)=-0.5D+0*WN2*A(2,3)
      A(3,4)=-0.5D+0*WN2*A(5,3)
      A(3,5)=-0.5D+0*WN2*A(4,3)
      A(4,1)=-A(2,5)
      A(4,2)=-A(1,5)
      A(5,1)=-A(2,4)
      A(5,2)=-A(1,4)
      A(5,5)=-A(4,4)
      GO TO  1000
C
C DIFFERENTIAL EQUATION FOR EIGENFUNCTIONS
C
  430 CONTINUE
      A(1,1)=R1*(1.0D+0-(2.0D+0*EP*EL)/RVP)
      A(1,2)=1.0D+0/(PH*RVP)
      A(3,1)=-R1
      A(3,3)=2.0D+0*R1
      A(3,4)=1.0D+0/EM
      A(4,1)=-2.0D+0*R2*(RVP-XC*EM-(ET*EP*EL*EL)/RVP)
      A(4,2)=-R1*EP*EL/RVP
      A(4,3)=-FRQ2*RH+R2*(WN2*(RVP-ET*EP*EL*EL/RVP)-2.0D+0*XC*EM)
      A(2,1)=-FRQ2*RH-2.0D+0*A(4,1)
      A(1,3)=-WN2*A(4,2)
      A(2,2)=-A(1,1)
      A(2,3)=WN2*A(4,1)
      A(2,4)=-WN2*A(3,1)
      A(4,4)=-A(3,3)
      GO TO  1000
C
C LIQUID LAYER
C
  450 CONTINUE
      A(1,1)=R1*(1.0D+0-2.0D+0*EP)
      A(1,2)=1.0D+0/(PH*EL)
      A(4,1)=-2.0D+0*R2*EL*(1.0D+0-ET*EP)
      A(4,2)=-R1*EP
      A(4,3)=-FRQ2*RH+R2*WN2*EL*(1.0D+0-ET*EP)
      A(1,3)=-WN2*A(4,2)
      A(2,1)=-FRQ2*RH-2.0D+0*A(4,1)
      A(2,3)=WN2*A(4,1)
C
      A(1,1)=A(1,1)-(A(1,3)*A(4,1))/A(4,3)
      A(1,2)=A(1,2)-(A(1,3)*A(4,2))/A(4,3)
      A(2,1)=A(2,1)-(A(2,3)*A(4,1))/A(4,3)
      A(2,2)=-A(1,1)
      GO TO  1000
C
C EXIT
C
 1000 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE  INTGND
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/MODEL /
     1              RAD,TMASS,GS,RHO0,GA,
     2              HI(500),RHOI(500),VPI(500),VSI(500),NAME(20),
     3              XII(500),PHII(500),ETAI(500),qbetai(500),
     4              D(1000),H(1000),RHO(1000),ELAMB(1000),EMU(1000),
     5              GR(1000),XI(1000),PHI(1000),ETA(1000),qbeta(1000),
     6              DMAX,ICONT,NLAY,ISO,LAME,NDIV,ILAY,idpmi,idpmax
      COMMON/VALUE /MODE,IERROR,I,ISTEP,J,K,L,N,NSOL,ISUM,KMAX,
     1              NMAX,IBOTM,ITOP,LY,LD,ACR,ELLIP,
     2              RMAX,RC,HC,H1,H2,ALP,WN,WN2,C,C2,FRQ,FRQ2,T,
     3              U,ENGY,DELTA,
     4              A(6,6),F(20),Y(6),repon,njob,nper,jyy,jdd,npar,
     5			ncouc,hdeb,hmin,hmax
C
C TO COMPUTE INTEGRANDS OF ENERGY INTEGRALS
C MANTLE WAVE (SPHERICAL MODEL)
C
      DO  1  IS=1,ISUM
      F(IS)=0.0
    1 CONTINUE
C
      RH=RHO(I)
      EL=ELAMB(I)
      EM=EMU(I)
      XC    =1.0D+0+XI(I)
      PH    =1.0D+0+PHI(I)
      ET    =1.0D+0+ETA(I)
      RVP   =EL+EM*2.0D+0
      EP=ET/PH
      R2=RC*RC
C
      GO TO  (1000,1000,300,400,300,1000),MODE
C
C MANTLE LOVE WAVE (SPHERICAL MODEL)
C
  300 CONTINUE
      IF(EM.EQ.0.0)  GO TO  1000
C         IF(Y(1).LT.1E-15) Y(1)=0.
C         IF(Y(2).LT.1E-15) Y(2)=0.
      F(1)=RH*Y(1)*Y(1)
      F(3)=XC*EM*((Y(1)/RC)**2)
      F(2)=(WN2-2.0D+0)*F(3)+(Y(2)*Y(2))/EM
      F(3)=RAD*RAD*F(3)
      GO TO  1000
C
C MANTLE RAYLEIGH WAVE (SPHERICAL MODEL)
C
  400 CONTINUE
      IF(EM.EQ.0.0)  GO TO  450
C
C SOLID LAYER
C
      W=RVP-EP*ET*EL*EL/RVP
      WW    =2.0D+0*Y(1)-WN2*Y(3)
      Y33=Y(3)*Y(3)
      Y44=Y(4)*Y(4)
      F(1)=RH*(Y(1)*Y(1)+WN2*Y33)
      F(2)=(W-XC*EM)*WW*WW/R2
     1    +Y(2)*Y(2)/(PH*RVP)
     2    +WN2*(WN2-2.0D+0)*XC*EM*Y33/R2
     3    +WN2*Y44/EM
      F(3)=-FRQ2*RH*Y33
     1    -2.0D+0*W*Y(3)*WW/R2
     2    -2.0D+0*EP*EL*Y(2)*Y(3)/(RVP*RC)
     3    +2.0D+0*XC*EM*Y(3)*(2.0D+0*Y(1)-Y(3))/R2
     4    +Y44/EM
      F(3)=RAD*RAD*F(3)
      GO TO  1000
C
C LIQUID LAYER
C
  450 CONTINUE
      W     =EL*(1.0D+0-EP*ET)
      A(4,3)=FRQ2*RH-WN2*W/R2
      A(4,1)=-2.0D+0*W/R2
      A(4,2)=-EP/RC
      Y3=(A(4,1)*Y(1)+A(4,2)*Y(2))/A(4,3)
      WW    =2.0D+0*Y(1)-WN2*Y3
      Y33=Y3*Y3
      F(1)=RH*(Y(1)*Y(1)+WN2*Y33)
      F(2)=W*WW*WW/R2
     1    +Y(2)*Y(2)/(PH*EL)
      F(3)=-FRQ2*RH*Y33
     1    -2.0D+0*W*Y3*WW/R2
     2    -2.0D+0*EP*Y(2)*Y3/RC
      F(3)=RAD*RAD*F(3)
      GO TO  1000
C
C EXIT
C
 1000 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE  DENERG
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/MODEL /
     1              RAD,TMASS,GS,RHO0,GA,
     2              HI(500),RHOI(500),VPI(500),VSI(500),NAME(20),
     3              XII(500),PHII(500),ETAI(500),qbetai(500),
     4              D(1000),H(1000),RHO(1000),ELAMB(1000),EMU(1000),
     5              GR(1000),XI(1000),PHI(1000),ETA(1000),qbeta(1000),
     6              DMAX,ICONT,NLAY,ISO,LAME,NDIV,ILAY,idpmi,idpmax
      COMMON/VALUE /MODE,IERROR,I,ISTEP,J,K,L,N,NSOL,ISUM,KMAX,
     1              NMAX,IBOTM,ITOP,LY,LD,ACR,ELLIP,
     2              RMAX,RC,HC,H1,H2,ALP,WN,WN2,C,C2,FRQ,FRQ2,T,
     3              U,ENGY,DELTA,
     4              A(6,6),F(20),Y(6),repon,njob,nper,jyy,jdd,npar,
     5			ncouc,hdeb,hmin,hmax
      COMMON/SOL/  YN(6,1000),YB(6,3,20),SUM(20),Q(3,21),
     1             WNB(100),TT(100),CC(100),UU(100),ENG(100),ELL(100),
     2             AC(100)
C
C ENERGY INTEGRALS IN THE INNERMOST HOMOGENEOUS SPHERE
C SINCE THE CONTRIBUTION FROM THE SPHERE R.LT.R(IBOTM) IS ASSUMED TO BE
C  NEGLIGIBLE, THIS PROGRAM COMPUTES ONLY THE GROUP VELOCITY AND THE
C  ACCURACY OF THE VARIATIONAL EQUATION
C MANTLE WAVES
C
      pi=3.141592653
      U=SUM(3)/(C*SUM(1))
      ENGY=FRQ2*SUM(1)
      ACR=(SUM(2)/ENGY)-1.0
      WNB(LY)=WN
      TT(LY)=T
      CC(LY)=C
      UU(LY)=U
      ENG(LY)=ENGY
      AC(LY)=ACR
      ELL(LY)=ELLIP
      repon=1./sum(3)
C
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE  DCDP(DG)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/MODEL /
     1              RAD,TMASS,GS,RHO0,GA,
     2              HI(500),RHOI(500),VPI(500),VSI(500),NAME(20),
     3              XII(500),PHII(500),ETAI(500),qbetai(500),
     4              D(1000),H(1000),RHO(1000),ELAMB(1000),EMU(1000),
     5              GR(1000),XI(1000),PHI(1000),ETA(1000),qbeta(1000),
     6              DMAX,ICONT,NLAY,ISO,LAME,NDIV,ILAY,idpmi,idpmax
      COMMON/VALUE /MODE,IERROR,I,ISTEP,J,K,L,N,NSOL,ISUM,KMAX,
     1              NMAX,IBOTM,ITOP,LY,LD,ACR,ELLIP,
     2              RMAX,RC,HC,H1,H2,ALP,WN,WN2,C,C2,FRQ,FRQ2,T,
     3              U,ENGY,DELTA,
     4              A(6,6),F(20),Y(6),repon,njob,nper,jyy,jdd,npar,
     5			ncouc,hdeb,hmin,hmax
      COMMON/SOL/  YN(6,1000),YB(6,3,20),SUM(20),Q(3,21),
     1             WNB(100),TT(100),CC(100),UU(100),ENG(100),ELL(100),
     2             AC(100)
C
C PARTIAL DERIVATIVES W.R.T. RHO, VP, VS, XI, PHI, AND ETA
C MANTLE WAVE (SPHERICAL MODEL)
C
      DRO=0.0
      DVP=0.0
      DVS=0.0
      DXC=0.0
      DPH=0.0
      DET=0.0
      RH=RHO(I)
      EL=ELAMB(I)
      EM=EMU(I)
      RVP   =EL+EM*2.0D+0
      XC    =1.0D+0+XI(I)
      PH    =1.0D+0+PHI(I)
      ET    =1.0D+0+ETA(I)
      EP=ET/PH
      Y1=YN(1,N)
      Y2=YN(2,N)
      Y3=YN(3,N)
      Y4=YN(4,N)
      Y11=Y1*Y1
      Y22=Y2*Y2
      R2=RC*RC
      CU    =0.5D+0*C/(ENGY*U)
C
      GO TO  (1000,1000,300,400,300,1000),MODE
C
C MANTLE LOVE WAVE (SPHERICAL MODEL)
C
  300 CONTINUE
      IF(EM.EQ.0.0)  GO TO  900
      DRO=-FRQ2*RH*Y11+Y22/EM+(WN2-2.0D+0)*XC*EM*Y11/R2
      DVS=2.0D+0*(Y22/EM+(WN2-2.0D+0)*XC*EM*Y11/R2)
      DXC=(WN2-2.0D+0)*XC*EM*Y11/R2
      GO TO  900
C
C MANTLE RAYLEIGH WAVE (SPHERICAL MODEL)
C
  400 CONTINUE
      IF(EM.EQ.0.0)  GO TO  450
C
C SOLID LAYER
C
      WW    =2.0D+0*Y1-WN2*Y3
      Y33=Y3*Y3
      Y44=Y4*Y4
      DRO=-FRQ2*RH*(Y11+WN2*Y33)
     1   +(RVP-XC*EM-(EP*ET*EL*EL)/RVP)*WW*WW/R2
     2   +Y22/(PH*RVP)
     3   +WN2*(WN2-2.0D+0)*XC*EM*Y33/R2
     4   +WN2*Y44/EM
      DVP=2.0D+0*(((Y2+2.0D+0*ET*EM*WW/RC)**2)/(PH*RVP)
     1   +RVP*(1.0D+0-EP*ET)*WW*WW/R2)
      DVS=2.0D+0*(((-XC*EM+4.0D+0*EP*ET*EL*EM/RVP)*WW/R2
     1   -4.0D+0*EP*EM*Y2/(RVP*RC))*WW
     2   +WN2*(WN2-2.0D+0)*XC*EM*Y33/R2
     3   +WN2*Y44/EM)
      DXC=XC*EM*(-WW*WW+WN2*(WN2-2.0D+0)*Y33)/R2
      DPH=((Y2-ET*EL*WW/RC)**2)/(PH*RVP)
      DET=2.0D+0*EP*EL*WW*(Y2-ET*EL*WW/RC)/(RVP*RC)
      GO TO  900
C
C LIQUID LAYER
C
  450 CONTINUE
      A(4,3)=FRQ2*RH-WN2*EL*(1.0D+0-EP*ET)/R2
      A(4,1)=-2.0D+0*EL*(1.0D+0-EP*ET)/R2
      A(4,2)=-EP/RC
      Y3=(A(4,1)*Y1+A(4,2)*Y2)/A(4,3)
      Y33=Y3*Y3
      WW    =2.0D+0*Y1-WN2*Y3
      DRO=-FRQ2*RH*(Y11+WN2*Y33)
     1   +EL*(1.0D+0-EP*ET)*WW*WW/R2
     2   +Y22/(PH*EL)
      DVP=2.0D+0*(Y22/(PH*EL)+EL*(1.0D+0-EP*ET)*WW*WW/R2)
      DPH=((Y2-ET*EL*WW/RC)**2)/(PH*EL)
      DET=2.0D+0*EP*WW*(Y2-ET*EL*WW/RC)/RC
      GO TO  900
C
C
  900 CONTINUE
      YN(1,N)=CU*DRO
      YN(2,N)=CU*DVP
      YN(3,N)=CU*DVS
      YN(4,N)=CU*DXC
      YN(5,N)=CU*DPH
      YN(6,N)=CU*DET
      GO TO  1000
C
C EXIT
C
 1000 CONTINUE
      RETURN
      END
