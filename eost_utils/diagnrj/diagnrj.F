
C==========================================================================
                           PROGRAM DIAGNRJ
 
C  DIAGRAMME D'ENERGIE (PERIODE - VITESSE DE GROUPE), EXTRACTION DE LA
C  VITESSE DE GROUPE APPARENTE LE LONG DU MAXIMUM D'ENERGIE, DISPERSION
C  RESIDUELLE ET CORRECTION D'AMPLITUDE EVENTUELS, CALCUL DE VITESSE DE
C  PHASE PAR INTEGRATION ET POSSIBILITE DE FILTRAGE VARIABLE DU SIGNAL.

C  VERSION SUN COULEUR OU NOIR/BLANC & SOURIS : Hugues DUFUMIER Sep.1993
C  Cette derniere version permet d'inclure les corrections de source:
C  histoire (vue centroide L.P.) et radiation (phase du mecanisme).
C==========================================================================
c  Entrees:
c   - FISIG: fichier SAC contenant:
c     OBS: signal SAC (si possible tourne, et assez long).
c     ! Si le signal n'est pas corrige de l'instrument et de la source,
c     les vitesses de groupe et de phase obtenues seront apparentes !
c     NOBS,TAUO: nombre de points et pas du signal a etudier
c     TE,T0: heures de debut du signal et origine du seisme le meme jour
c     DKM: distance epicentrale en Km
c     Si ces donnees ne figurent pas dans le header, elles sont demandees.
c   - Les limites desirees du diagramme.
c   - La demi-duree de source pour correction du delai de groupe origine.
c   - Quelques options mais une part decisionnelle reduite.
c  Sorties:
c   - VGR(FRQ): courbe de vitesse de groupe le long du max d'energie
c   - VPH(FRQ): courbe de vitesse de phase integree (en option)
c   - FISIG.filt: signal apres filtrage (variable ou classique).
C  Convention de TF directe en + (Aki) : F(w)= S(f(t).exp(+iwt).dt)
C==========================================================================
C Fichiers utiles: tcoul,tgris,Uprem.Rayl,Uprem.Love,Cprem.Rayl,Cprem.Love.
C Compilation: "make diagnrj" (linker tsg.inc.o, Suncore couleur & PS, SAC)

      REAL OBS(16384)
      COMPLEX XC(16384),REP(16384)
      CHARACTER FISIG*40,FICH*40,STA*4,CMP*1,IVPH*1
      COMMON/STACMP/STA,CMP
c     liens avec dispersion
      REAL U(150),FRQ(54),PRD(54),VPH(54),VGR(54),X(99),V(99)
c     corrections de source
      COMMON/CORSOURS/TGM(54),TC
      PI=4.*ATAN(1.)

      PRINT*,' '
      PRINT*,'DIAGRAMME D''ENERGIE ET FILTRAGE VARIABLE H.Dufumier'
      PRINT*,'----------------------------------------------------'
      PRINT*,'A faire tourner en graphique (de preference couleur)'
      PRINT*,' '
      PRINT*,'Ce programme permet de visualiser l''energie d''un signal'
      PRINT*,'sur un diagramme en vitesse de groupe et periode (en fait'
      PRINT*,'en delai de groupe et frequence). En isolant les paquets'
      PRINT*,'d''interet par un filtrage variable, on obtient un signal'
      PRINT*,'propre ne contenant plus que l''onde desiree. Si de plus,'
      PRINT*,'le signal est corrige de la source et de l''instrument,'
      PRINT*,'il est possible de piquer le long du maximum d''energie'
      PRINT*,'la vitesse de groupe de l''onde, puis de l''integrer pour'
      PRINT*,'obtenir la vitesse de phase de propagation. Sinon les vi-'
      PRINT*,'tesses obtenues seront des vitesses apparentes incluant'
      PRINT*,'les delais de groupe a la source et/ou instrumentaux. Un'
      PRINT*,'moyen simplifie de s''affranchir des effets de source a'
      PRINT*,'longue-periode est de prendre le temps origine et la dis-'
      PRINT*,'tance epicentrale centroides. Pour le reste, le diagramme'
      PRINT*,'d''energie donnera les meilleurs resultats si le signal'
      PRINT*,'est deconvolue de l''instrument, ramene en acceleration,'
      PRINT*,'non filtre, et tourne suivant la polarisation de l''onde.'
      PRINT*,' '

C  CHOIX DIVERS

      PRINT*,'Periodes inferieure et superieure du diagramme (s) ?'
      READ(*,*) TDEB,TFIN
      if(tdeb.ge.tfin) then
        tfin0=tfin
        tfin=tdeb
        tdeb=tfin0
        endif

      PRINT*,'Vitesses de groupe inferieure & superieure (km/s) ?'
      READ(*,*) UINF0,USUP0
      if(uinf0.ge.usup0) then
        usup=usup0
        usup0=uinf0
        uinf0=usup
        endif

      PRINT*,'Voudrez-vous calculer la vitesse de phase par integration'
      PRINT*,'de la vitesse de groupe piquee dans le diagramme [o] ?'
      READ(*,'(A)') IVPH
      IF(IVPH.EQ.' ') IVPH='o'

C  BOUCLE PRINCIPALE
C  -----------------
 
  100 IF(NDIM.EQ.0) THEN
      PRINT*,'Fichier SIGNAL SAC a traiter ?'
      ELSE
      PRINT*,'Autre SIGNAL SAC a traiter (ou Return pour terminer) ?'
      ENDIF
      READ(*,'(A)') FISIG
      IF(FISIG.EQ.' ') GOTO 999

C  LECTURE DU SIGNAL SAC

      NDIM=16384
      CALL RIDSAC(FISIG,OBS,NOBS,TAUO,DKM,TE,T0,NDIM)

C  MOYENNE & DERIVE NULLES, APODISATION

      CALL DTRD(OBS,NOBS,OBS,1)
      CALL TAPER(OBS,NOBS,0.01,0.01,OBS)

C  CALCUL PARAMETRES POUR TF
      
      NP=INT(ALOG(FLOAT(NOBS-1))/ALOG(2.))+1
      LX=2**NP
      NF2=LX/2+1
      DF=1./(LX*TAUO)
      WRITE(*,1001) LX,TAUO,DF
 1001 FORMAT('Calcul TF sur',I5,' pts.  Dt=',F6.3,'s.  Df=',F8.6,' Hz')

C  SELECTION DES PERIODES utiles: de TDEB a TFIN avec une marge autour
c  les frequences seront en ordre croissant, les periodes decroissant

      IF(TDEB.LT.2.*TAUO) PRINT*,'! Periode de debut inferieure a Nyquis 
     *t --> on prend a partir de',2.*TAUO,'s !'
      TDEB=AMAX1(TDEB,2.*TAUO)
      IF(TFIN.GT.1./DF) PRINT*,'! Periode Max doit etre < Duree signal 
     *--> on prend Pmax=',1./DF
      TFIN=AMIN1(TFIN,1./DF)

      WRITE(*,1005) 'Spectre :',NF2-1,.5/TAUO,DF
 1005 FORMAT(A,I5,' frequences de 0 a',F7.3,'Hz par pas de',F8.5,'Hz')
c     limites: points les plus proches dans TDEB-TFIN + marges
      ifin0=int(1./TDEB/DF)
      ideb0=int(1./TFIN/DF)+1
   15 nfu0=ifin0-ideb0+1
c     on prend 1 freq. sur MUL pour avoir au max. 52 frequences dedans
      MUL=int(float(nfu0-1)/51.)+1
      NFU=int(float(nfu0-1)/float(MUL))+3
c     indices et periodes de marge
      ideb=ideb0-MUL
      ifin=ideb+(NFU-1)*MUL
      PMIN=1./float(ifin)/DF
      PMAX=1./float(ideb)/DF
      if(PMIN.lt.2.*TAUO) ifin0=ifin0-1
      if(PMAX.lt.0..or.PMAX.gt.600.) ideb0=ideb0+1
      if(PMIN.lt.2.*TAUO.or.PMAX.lt.0..or.PMAX.gt.600.) goto 15
      pasf=float(MUL)*DF
      DO 20 I=1,NFU
      VGR(I)=0.
      FRQ(I)=1./PMAX+FLOAT(I-1)*pasf
   20 PRD(I)=1./FRQ(I)
      WRITE(*,'(A,I3,A,2F7.2,A,I2,A)')'Selection',NFU,' periodes dans l'
     &'intervalle',PMIN,PMAX,' s (1 pt/',MUL,'):'
      WRITE(*,'(12F7.2)') (PRD(I),I=NFU,1,-1)

C  ENTREES POUR CORRECTIONS DE SOURCE (HISTOIRE ET RADIATION)

      PRINT*,' '
      PRINT*,'Si vous n''avez pas deconvolue le signal de la fonction'
      PRINT*,'source, ni utilise le t0 et la distance epicentrale cen-'
      PRINT*,'troides, entrez la DEMI-DUREE DE SOURCE pour corriger le'
      PRINT*,'diagramme du delai de groupe a la source (ou 0 sinon).'
      READ(*,*) HALFDUR
c     temps origine centroide
      TC=T0+HALFDUR

   35 PRINT*,'Si vous voulez corriger le diagramme de la radiation de'
      PRINT*,'la source (phase induite par le mecanisme), entrez le'
      PRINT*,'fichier SAC contenant cette phase (sinon Return).'
      READ(*,'(A)') FICH
      DO 36 I=1,NFU
   36 TGM(I)=0.
      IF(FICH.NE.' ') THEN
       CALL RSAC1(FICH,V,NV,X(1),DX,99,IERR)
       DO 37 I=1,NV
       X(I)=X(1)+FLOAT(I-1)*DX
   37  V(I)=V(I)/2./PI
       IF(IERR.NE.0) PRINT*,'Fichier Phase du Meca incorrect. Autre ?'
       IF(IERR.NE.0) GOTO 35
       CALL INTERPOL(V,NV,X,.001,TGM,FRQ,NFU,1)
       CALL WSAC2('delaigr.rad',tgm,nfu,frq,ierr)
       IF(IERR.EQ.0) PRINT*,'Delai de groupe associe dans delaigr.rad'
       IF(IERR.NE.0) PRINT*,'! Probleme dans ecriture delai de groupe'
       ENDIF

C  CALCUL DES 150 VITESSES DU DIAGRAMME

      UINF=UINF0
      USUP=USUP0
      TN=TE+FLOAT(NOBS-1)*TAUO
      TUINF=DKM/UINF+TC
      IF(TN.LT.TUINF) THEN     
        UINF=DKM/(TN-TC)
        PRINT*,'! Signal trop court pour intervalle de Vit.Groupe desire
     & --> on prend Uinf=',UINF
        IF(UINF.GT.USUP) PRINT*,'! Diagramme impossible: Uinf > Usup. !'
        UINF=AMIN1(UINF,USUP)
        TUINF=TN
        ENDIF
      TUSUP=DKM/USUP+TC
      IF(TE.GT.TUSUP) THEN
        USUP=DKM/(TE-TC)
        PRINT*,'! Signal debutant trop tard pour intervalle de Vit.Group
     &e desire --> on prend Usup=',USUP
        TUSUP=TE
        ENDIF
c     on fixe ici le pas en U pour avoir maxi 120 points verticalement
      UPAS=(USUP-UINF)/119.
      UPAS=AMAX1(UPAS,(USUP**2)*TAUO/DKM)
      KU=INT((USUP-UINF)/UPAS)+1
      DO 10 I=1,KU
   10 U(I)=UINF+FLOAT(I-1)*UPAS

C  TRANSFORMEE DE FOURIER DIRECTE temps --> spectre d'energie frequence

      DO 30 I=1,LX
   30 XC(I)=CMPLX(OBS(I)*TAUO,0.)
      CALL FOURT(XC,NP,+1)

C  DECONVOLUTION REPONSE INSTRUMENT --> SPECTRE ACCELERATION EN M

      CALL INSPZ(REP,LX,DF)
      DO 31 I=2,NF2
   31 XC(I)=XC(I)/REP(I)
      XC(1)=CMPLX(0.,0.)

C  CHOIX COURBE DE BASE POUR POINTAGE

      PRINT*,' '
      PRINT*,'TRACE DU DIAGRAMME D''ENERGIE A ',STA,CMP
      PRINT*,'Pour une recherche automatique de la vitesse de groupe,'
   40 PRINT*,'entrez une COURBE DE BASE SAC (par exemple Uprem.Rayl [R]'
      PRINT*,'ou Uprem.Love [L]) sinon RETURN POUR POINTER a la souris.'
      READ(*,'(A)') FICH
      IF(FICH.NE.' ') THEN
       IF(FICH.EQ.'R') FICH='Uprem.Rayl'
       IF(FICH.EQ.'L') FICH='Uprem.Love'
       CALL RSAC2(FICH,V,NV,X,99,IERR)
       IF(IERR.NE.0) GOTO 40
       CALL INTERPOL(V,NV,X,.001,VGR,FRQ,NFU,0)
       ENDIF

C  DIAGRAMME D'ENERGIE (trace dans l'unite du signal)

      CALL FDRNRJ(XC,NOBS,TAUO,DKM,TE,T0,FRQ,NFU,U,KU,VGR,IT1,IT2)

C  SAUVEGARDE SAC DE LA DISPERSION POINTEE

c     on vire la longue-periode de marge de VGR(FRQ)
      I1=MAX0(2,IT1)
      NFU=IT2-I1+1
      PRINT*,' '

c     sauvegarde vitesse de groupe
      L=LNBLNK(STA)
      FICH='Unrj.'//STA(:L)//'.'//CMP
      CALL WSAC2(FICH,VGR(I1),NFU,PRD(I1),IERR)
      IF(IERR.EQ.0) PRINT*,'Sauvegarde SAC vitesse de groupe dans ',FICH

c     integration et sauvegarde vitesse de phase
      IF(IVPH.EQ.'o') THEN
      CALL VGVPH(CMP,VGR(I1),FRQ(I1),NFU,VPH(I1))
      FICH='Cnrj.'//STA(:L)//'.'//CMP
      CALL WSAC2(FICH,VPH(I1),NFU,PRD(I1),IERR)
      IF(IERR.EQ.0) PRINT*,'Sauvegarde SAC vitesse de phase dans ',FICH
      ENDIF
      PRINT*,'Utilisez la macro SAC "dispnrj" pour tracer ces vitesses.'

C  SAUVEGARDE DU SIGNAL FILTRE

      DO 50 I=1,LX
   50 OBS(I)=REAL(XC(I))
      FISIG=FISIG(:LNBLNK(FISIG))//'.filt'
      CALL SAUVSAC(FISIG,OBS,NOBS,TAUO,TE,T0,DKM)
      PRINT*,'Signal filtre sauve dans ',FISIG
      PRINT*,'Comparez-le au signal brut avec la macro SAC "sigfilt".'
      PRINT*,' '

C  NOUVEAU SIGNAL
C  --------------

      GOTO 100

  999 CALL findes
      PRINT*,'_________________________________________________________'
      PRINT*,' '
      PRINT*,'! Pensez a detruire les fichiers PostScript apres usage !'
      PRINT*,'_________________________________________________________'
      END

C=======================================================================

      SUBROUTINE RIDSAC(name,X,NPT,dt,DKM,TE,T0,ndim)

C=======================================================================
c     LECTURE FICHIER SAC ET INFOS HEADER UTILES. APPELLE LIBRAIRIE SAC.
      CHARACTER name*30,event*16,sta*4,type*8,comp*8,cmp*1
      COMMON/STACMP/sta,cmp
      COMMON/DATE/iy,id
      REAL X(*)

    1 call rsac1(name,x,NPT,beg,dt,ndim,nerr)
      if(nerr.ne.0) then
        PRINT*,'Fichier incorrect. Fichier signal SAC a lire ?'
        READ(*,'(a)') name
        goto 1
        endif
      if(NPT.eq.ndim) print*,'! On a coupe a ',ndim,' points ! Decimez l
     &e signal ou coupez-le avant traitement !'
      l=lnblnk(name)
      call getfhv('DIST',DKM,nerr)
      if(DKM.eq.-12345..or.nerr.ne.0) then
        PRINT*,'Distance epicentrale inconnue, entrez-la en km.'
        READ(*,*) DKM
        endif
      call getnhv('NZYEAR',iy,nerr)
      call getnhv('NZJDAY',id,nerr)
      call getnhv('NZHOUR',ih,nerr)
      ihbeg=iho
      call getnhv('NZMIN',im,nerr)
      call getnhv('NZSEC',is,nerr)
      call getnhv('NZMSEC',ims,nerr)
      sec=float(is)+float(ims)/1000
      call getihv('IZTYPE',type,nerr)
      call getfhv('O',T0,nerr)
      call getkhv('KEVNM',event,nerr)
      call getkhv('KSTNM',sta,nerr)
      if(sta(1:4).eq.'-123'.or.nerr.ne.0) then
        PRINT*,'Station inconnue, entrez son code [BID].'
        READ(*,'(a)') sta
        if(sta.eq.' ') sta='BID'
        endif
      comp=' '
      call getkhv('KCMPNM',comp,nerr)
      if(comp(1:6).ne.'-12345'.and.comp.ne.' ') then
        cmp=comp(lnblnk(comp):lnblnk(comp))
        else
        cmp=name(l:l)
        endif
      print*, event,'   Date:',iy,id
      print*,'Station: ',sta,' Composante: ',cmp,'  Distance (km):',DKM
      print*,'Debut:',beg,'s apres ',ih,im,sec,'   Nbre Pts,Pas:',NPT,dt
      TE=3600.*ih+60.*im+sec+beg
      if(T0.eq.-12345..and.type.ne.'IO') then
        PRINT*,'Heure Origine du Seisme non trouvee dans header.'
        PRINT*,'Entrez-la en heure min sec du meme jour.'
        READ(*,*) IH0,IM0,S0
        T0=3600.*IH0+60.*IM0+S0
        else if(type.eq.'IO') then
        T0=3600.*ih+60.*im+sec
        else
        T0=3600.*ih+60.*im+sec+T0
        endif
      ih0=int(T0/3600.)
      im0=int((T0-3600.*float(ih0))/60.)
      s0=T0-3600.*float(ih0)-60.*float(im0)
      dt0=T0-TE
      print*,'Heure origine:',ih0,im0,s0,' (T0-Te =',dt0,'s).'
      RETURN
      END

C=======================================================================

      SUBROUTINE SAUVSAC(FICH,X,NX,PX,TE,T0,DKM)

C=======================================================================
C  Sauvegarde d'un signal reel au format SAC binaire.
C  Header minimum comprenant les heures de debut et le temps origine.

      CHARACTER FICH*40,STA*4,CMP*1
      COMMON/STACMP/STA,CMP
      COMMON/DATE/iy,id
      REAL X(*)

      IH=INT(TE/3600.)
      IM=INT((TE-3600.*FLOAT(IH))/60.)
      SEC=TE-3600.*FLOAT(IH)-60.*FLOAT(IM)
      IS=SEC
      IMS=1000*(SEC-IS)
      DT0=T0-TE
      CALL NEWHDR
      CALL SETNHV('NPTS',NX,IERR)
      CALL SETLHV('LEVEN',.TRUE.,IERR)
      CALL SETFHV('B',0.,IERR)
      CALL SETFHV('E',FLOAT(NX-1)*PX,IERR)
      CALL SETFHV('DELTA',PX,IERR)
      CALL SETIHV('IFTYPE','ITIME',IERR)
      CALL SETIHV('IZTYPE','IB',IERR)
      CALL SETNHV('NZYEAR',iy,IERR)
      CALL SETNHV('NZJDAY',id,IERR)
      CALL SETNHV('NZHOUR',IH,IERR)
      CALL SETNHV('NZMIN',IM,IERR)
      CALL SETNHV('NZSEC',IS,IERR)
      CALL SETNHV('NZMSEC',IMS,IERR)
      CALL SETFHV('O',DT0,IERR)
      CALL SETFHV('DIST',DKM,IERR)
      CALL SETKHV('KSTNM',STA,IERR)
      CALL SETKHV('KCMPNM',CMP//'.filt',IERR)
      CALL SETKHV('KEVNM',FICH(1:16),IERR)
      CALL WSAC0(FICH,X,X,IERR)
      IF(IERR.NE.0) PRINT*,'Probleme SAC a l''ecriture du signal filtre'
      RETURN
      END

C=======================================================================
 
      SUBROUTINE FDRNRJ(XC,NPT,PAS,DKM,TE,T0,FRQ,KT,U,KU,UM,IT1,IT2)

C=======================================================================
C  DIAGRAMME D'ENERGIE (PERIODE - VITESSE DE GROUPE), EXTRACTION DE LA 
C  VITESSE DE GROUPE APPARENTE LE LONG DU MAXIMUM D'ENERGIE EN AUTOMA-
C  TIQUE SI ON A UNE COURBE DE BASE (UM) ET/OU A LA SOURIS, DISPERSION
C  RESIDUELLE ET CORRECTION D'AMPLITUDE EVENTUELS, DESSIN DES LIMITES
C  DU FILTRAGE VARIABLE ET FILTRAGE VARIABLE OU CLASSIQUE DU SIGNAL.
C
C  ADAPTATION DES VERSIONS CARA 1987-1989. REFERENCES: CARA 1973,76,78.
C  VERSION SUN COULEUR ET SOURIS POUR DIAGNRJ: Hugues DUFUMIER Mai 1993
C=======================================================================
c Entrees:
c     XC: spectre (si possible tourne et corrige de l'instrument)
c     NPT,PAS: nombre de points et pas du signal a etudier
c     TE: heure de debut d'enregistrement du signal
c     T0: heure origine du seisme (en sec. le meme jour)
c     DKM: distance epicentrale en Km
c     FRQ,KT: frequences (croissantes) du diagramme et nombre de freq.
c     U,KU: vitesses de groupe et nombre de vit.de groupe du diagramme.
c     UM(FRQ): courbe de vitesse de groupe de base pour pointe du max.
c Common Corrections de source:
c     TGM(FRQ): delai de groupe de la radiation du mecanisme
c     TC: heure origine centroide (en sec. le meme jour)
c Sorties:
c     UM(FRQ): courbe de vitesse de groupe le long du max d'energie.
c     IT1,IT2: indices limites de selection de la courbe UM(FRQ).
c     XC est le signal filtre en sortie (prendre sa partie reelle).
C=======================================================================

c     spectres
      COMPLEX XC(16384),XD(16384),YC(16384),CI2
c     filtres
      REAL CF(32768),SF(32768),FO(501)
c     diagrammes: contenu et axes
      REAL AMP(150,54),PHASE(150,54),FRQ(54),U(150)
c     courbes intermediaires et sorties
      REAL UM(54),AM(54),TAU(54),T1(999),T2(999)
c     corrections de source
      COMMON/CORSOURS/TGM(54),TC
      INTEGER NBAND(54)
      PI=4.*ATAN(1.)
      NDIM=16384
      anorf=sqrt(pi/(6.*alog(10.)))

C  PREMIER PASSAGE: IBAND=0 FILTRAGE CLASSIQUE
      IBAND=0

C  RECALCUL DES PARAMETRES DE LA TRANSFORMEE DE FOURIER
      IP1=INT(ALOG(FLOAT(NPT-1))/ALOG(2.))+1
      NP=2**IP1
      DNU=1./FLOAT(NP)/PAS

      UPAS=(U(KU)-U(1))/FLOAT(KU-1)
      PADIA=DKM*UPAS/(U(KU)**2)
      N3=1./(DNU*PADIA)
      IP3=1+ALOG(FLOAT(N3-1))/ALOG(2.)
      NP3=2**IP3

C  CALCUL DU FILTRE GAUSSIEN FO
      LBAND=500
      MIBAND=251
      BETA=-13.812
      DO  5 I=1,MIBAND
      FO(I)=EXP(BETA*((FLOAT(I-MIBAND)/FLOAT(LBAND))**2))
    5 FO(LBAND-I+2)=FO(I)

C  CHOIX LARGEUR FILTRAGE

   10 DO 15 I=1,NP
   15 XD(I)=XC(I)

      PRINT*,'Entrez la largeur de bande relative a -30db (1 est usuel)'
      READ(*,*) PAR
      WRITE(*,1001)'Largeur du filtre a -30 Db :',PAR,'*frequence'
 1001 FORMAT(A,F4.2,A)
      IT1=1
      IT2=KT

C  PARTIE PRINCIPALE: FILTRAGE ET OPTIONS
  100 DO 102 I=IT1,IT2
      IR=NINT(FRQ(I)*PAR/DNU+1.)/2
  102 NBAND(I)=2*IR+1

C  CORRECTION D'AMPLITUDE
      IF(IBAND.EQ.1) CALL SZA(AM,FRQ,DNU,IT1,KTT,XD,NP)

C  DISPERSION RESIDUELLE
      IF(IABS(IBAND).EQ.1)
     &CALL SZD(UM,FRQ,DNU,DKM,IT1,KTT,CF,SF,TAU,NBAND)

C  FILTRAGE MULTIPLE
      PRINT*,'Construction du Diagramme: Filtrage multiple. Patience...'
      DO 105 K=1,KT
      DO 105 I=1,KU
  105 AMP(I,K)=0.

C     GRANDE BOUCLE DE QUADRILLAGE
      DO 120 K=1,KT
      DOM=PI*FLOAT(NBAND(K)-1)*DNU
      N=NBAND(K)
      IF(N.GT.NDIM) PRINT*,'! DIM. TABLEAUX DEPASSEES ! LARGEUR DE BANDE
     * TROP GRANDE !'
      NN=N-1
      NF1=NINT(FRQ(K)/DNU+1.-FLOAT(N-1)/2.)
      IF(NF1.GE.1) GOTO 114
      N=NINT(FRQ(K)/DNU)
      N=2*N+1
      IF(N.GT.NDIM) PRINT*,'! DIM. TABLEAUX DEPASSEES ! PAS EN FREQUENCE
     * TROP PETIT !'
      NN=N-1
      NF1=1
  114 NC=NF1+NN/2
c     application du filtre
      DO 122 J=1,N
      AJ=FLOAT(J-1)*500./FLOAT(N-1)+1.
      J1=AJ
      AJ1=J1
      J2=J1+1
      IF(J2.EQ.502) THEN
       FIL=FO(501)
      ELSE
       FO1=FO(J1)
       FIL=FO1+(FO(J2)-FO1)*(AJ-AJ1)
       I=NF1+J-1
      ENDIF
      IF(IBAND.EQ.0.OR.K.LT.IT1.OR.K.GT.IT2) THEN
       YC(J)=FIL*XD(I)
      ELSE
       FLR=FIL*(CF(I)*CF(NC)+SF(I)*SF(NC))
       FLI=FIL*(SF(NC)*CF(I)-SF(I)*CF(NC))
       YC(J)=CMPLX(FLR,-FLI)*XD(I)
      ENDIF
  122 CONTINUE
      IP2=1+ALOG(FLOAT(N-1))/ALOG(2.)
      NP2=2**IP2
      N2=N+1
      IF(NP2.LT.NP3) THEN
        NP2=NP3
        IP2=IP3
        ENDIF
      DO 124 I=N2,NP2
  124 YC(I)=0.
c     signal temporel reconstitue par periode
      CALL FOURT(YC,IP2,-1)
      PJ=PAS*FLOAT(NP)/FLOAT(NP2)

C     CALCUL AMPLITUDES ET PHASES SUR LA GRILLE DU DIAGRAMME
      DO 120 I=1,KU
c     recherche du point arrivant a la vitesse U dans le signal YC
c     en referentiel corrige des delais de groupe de source
      TPS=DKM/U(I)-(TE-TC-TGM(K))
      J=NINT(TPS/PJ+1.)
      DTP=TPS-FLOAT(J-1)*PJ
      IF(IABS(IBAND).EQ.1.AND.K.GE.IT1.AND.K.LE.IT2) 
     & J=NINT((TPS-TAU(K))/PJ+1.)
  125 IF(J.LT.1) J=J+NP2
      IF(J.GT.NP2) J=J-NP2
      IF(J.LT.1.OR.J.GT.NP2) GOTO 125
      TPSJ=FLOAT(J-1)*PJ
      TPS=TPSJ+DTP
      J1=J-1
      J2=J
      J3=J+1
      IF(J1.LT.1) J1=NP2
      IF(J3.GT.NP2) J3=1
c     phase station et dephasage
      P1=ATAN2(AIMAG(YC(J1)),REAL(YC(J1)))+DOM*FLOAT(J1-1)*PJ
      P2=ATAN2(AIMAG(YC(J2)),REAL(YC(J2)))+DOM*FLOAT(J2-1)*PJ
      P3=ATAN2(AIMAG(YC(J3)),REAL(YC(J3)))+DOM*FLOAT(J3-1)*PJ
      PHASE(I,K)=P2
      DP=(P3-P1)/PJ/2.
c     amplitude: moyenne sur j-1,j,j+1
      A1=CABS(YC(J1))
      A2=CABS(YC(J2))
      A3=CABS(YC(J3))
      X1=TPSJ-PJ
      X2=TPSJ
      X3=TPSJ+PJ
      Z1=(TPS-X3)*(TPS-X2)/((X1-X3)*(X1-X2))
      Z2=(TPS-X1)*(TPS-X3)/((X2-X1)*(X2-X3))
      Z3=(TPS-X1)*(TPS-X2)/((X3-X1)*(X3-X2))
      AMREL=A1*Z1+A2*Z2+A3*Z3
      AMP(I,K)=AMREL/(FRQ(K)*PAR*anorf)
c     commentaire mis pour faire aussi figurer la corr.d'ampl. sur diagramme
c     IF(IBAND.EQ.1) AMP(I,K)=AMP(I,K)*AM(K)
      IF(IBAND.NE.0.AND.(K.LT.IT1.OR.K.GT.IT2)) AMP(I,K)=0.
  120 CONTINUE

C  DESSIN COULEUR DU DIAGRAMME D'AMPLITUDES

c     UM = courbe de base eventuelle en entree au 1er passage
c          courbe de vitesse le long du maximum d'energie ensuite
  150 CALL SZGRAPH(AMP,U,FRQ,KU,KT,IT1,KTT,UM,IBAND)

C  RECHERCHE DU MAXIMUM D'ENERGIE AUTOUR DU TEMPS DE GROUPE TAU DE BASE
      IF(IABS(IBAND).NE.1) THEN
        IT2=IT1+KTT-1
        DO 160 I=IT1,IT2
  160   TAU(I)=DKM/UM(I)
        ENDIF

C  LARGEUR DE RECHERCHE: 0.4 km/s EN AUTOMATIQUE AU PREMIER PASSAGE,
C  0.2 km/s ENSUITE OU SI PIQUAGE A LA SOURIS, MODIFIABLE.
      IF(IBAND.EQ.0) DU=0.4
      IF(IABS(IBAND).EQ.1.OR.IBAND.EQ.2) DU=0.2
      CALL SZMAX(AMP,PHASE,U,FRQ,KU,IT1,KTT,TAU,DKM,DU,UM,AM)

C  DESSIN DES LIMITES DU FILTRAGE VARIABLE
      CALL SZOND(FRQ(IT1),UM(IT1),KTT,DNU,DKM,T1,T2)

C  SECOND PASSAGE EVENTUEL AVEC OPTIONS
      PRINT*,'---> ?'
      IF(IABS(IBAND).EQ.1) THEN
      PRINT*,' 0: Revenir au premier filtrage ou modif. largeur bande;'
      ELSE
      PRINT*,' 0: Changer la largeur de bande du filtrage multiple;'
      PRINT*,'-1: Dispersion residuelle; 1: Disp.res.+ Correct.amplit.;'
      PRINT*,' 2: Refaire le pointage du max. par un piquage souris;'
      ENDIF
      PRINT*,' 3: Appliquer Filtrage Variable; 4: Fin pour ce diagramme'
      READ(*,*) IBAND
      IF(IBAND.EQ.0) GOTO 10
      IF(IABS(IBAND).EQ.1) GOTO 100
      IF(IBAND.EQ.2) GOTO 150

C  FILTRAGE VARIABLE

      IF(IBAND.EQ.3) THEN
c      calcul des limites du filtrage variable
       do 200 i=1,KT
       do 200 j=1,KU
       if(i.lt.IT1.or.i.gt.IT2) AMP(j,i)=0.
  200  continue
       do 202 i=IT1,IT2
       it=1+NINT(FRQ(i)/DNU)-NINT(FRQ(IT1)/DNU)
       tm=(T1(it)+T2(it))/2.
       dt=(T2(it)-T1(it))/2.
       do 202 j=1,KU
       dst=DKM/U(j)-tm
       if(abs(dst).lt.dt) ap=1.-(dst/dt)**4
       if(abs(dst).ge.dt) ap=0.
       AMP(j,i)=ap*AMP(j,i)
  202  continue
c      visualisation graphique du filtrage
       PRINT*,'FILTRAGE VARIABLE en dispersion.'
       CALL SZGRAPH(AMP,U,FRQ,KU,KT,IT1,KTT,UM,3)
c      filtrage variable -> XC = sortie en temps du signal filtre
       DTES=TE-T0
       CALL NPAFNRJ(XC,NP,PAS,DTES,FRQ(IT1),FRQ(IT2),T1,T2)
      ELSE
c      sinon filtrage classique en periode du signal
       WRITE(*,2000) 1./FRQ(IT2),1./FRQ(IT1)
 2000  FORMAT('FILTRAGE CLASSIQUE PASSE-BANDE de',F7.2,' a',F7.2,'s')
       PRINT*,'(filtre Butterworth d''ordre 4 dephasant)'
       FH=FRQ(IT1)
       FB=FRQ(IT2)
       do 210 i=1,KT
       F=FRQ(i)
       CI2=CMPLX(0.,-SQRT(2.)*F)
       do 210 j=1,KU
c      Module du Butterworth d'ordre 4
       AMP(j,i)=AMP(j,i)/CABS((1.+CI2/FB-(F/FB)**2)**4)
  210  AMP(j,i)=AMP(j,i)*CABS(((F/FH)**2/(1.+CI2/FH-(F/FH)**2))**4)
c      visualisation graphique du filtrage
       CALL SZGRAPH(AMP,U,FRQ,KU,KT,IT1,KTT,UM,-3)
c      retour en temps
       NF2=NP/2+1
       XC(1)=CMPLX(0.,0.)
       DO 220 I=2,NF2
       F=FLOAT(I-1)*DNU
c      filtre passe-bande Butterworth d'ordre 4
       CI2=CMPLX(0.,-SQRT(2.)*F)
       XC(I)=XC(I)/(1.+CI2/FB-(F/FB)**2)**4
       XC(I)=XC(I)*(-(F/FH)**2/(1.+CI2/FH-(F/FH)**2))**4
       XC(I)=XC(I)*DNU
  220  XC(NP-I+2)=CONJG(XC(I))
       XC(NF2)=REAL(XC(NF2))
       CALL FOURT(XC,IP1,-1)
      ENDIF

      RETURN
      END

c-----------------------------------------------------------------------
      SUBROUTINE SZA(AM,F,DNU,IT1,KTT,XD,NP)
c-----------------------------------------------------------------------
C     CORRECTION D'AMPLITUDE, CELLE-CI ETANT CONSERVEE DANS AM
      COMPLEX XD(16384)
      REAL F(54),AM(54),DU(54),C0(55),C1(55),C2(55),C3(55)
      REAL*8 R(440)

      BRUIT=70.
      PRINT*,'Correction d''amplitude jusqu''a -70 Db'
      BR=10.**((99.-BRUIT)/20.)
      IT2=IT1+KTT-1
      F2=F(IT2)
      F1=F(IT1)
c     amplitudes extremes, marge de lissage
      AMAX=AM(IT1)
      DO 1 I=IT1,IT2
    1 IF(AM(I).GT.AMAX) AMAX=AM(I)
c     A1=AMAX/BR est le bruit estime
      A1=AMAX/BR
      DO 2 I=1,KTT
    2 DU(I)=AMAX/500.
      S=FLOAT(KTT)
      SO=S/2.
      ILI=0
c     lissage des amplitudes
   20 CALL LISSE(KTT,F(IT1),AM(IT1),DU,S,C0,C1,C2,C3,R)
      IF(C0(1).NE.0.) GOTO 30
      WRITE(*,*)'! DIFFICULTE DANS LISSAGE AMPLITUDES !'
      ILI=ILI+1
      S=SO*FLOAT(ILI)
      IF(ILI.LT.10) GOTO 20
      WRITE(*,*)'! IMPOSSIBLE A LEVER AVEC FLUCTUATIONS DE S !'
c     correction d'amplitudes sur le spectre
   30 DO 50 I=1,NP
      FO=FLOAT(I-1)*DNU
      IF(FO.GT.F1.AND.FO.LT.F2) A=SMOO(KTT,F(IT1),C0,C1,C2,C3,FO)
      IF(FO.LE.F1) A=SMOO(KTT,F(IT1),C0,C1,C2,C3,F1)
      IF(FO.GE.F2) A=SMOO(KTT,F(IT1),C0,C1,C2,C3,F2)
      IF(A.LT.A1) A=A1
      IF(A.GT.AMAX) A=AMAX
   50 XD(I)=XD(I)/A
c     restitution des amplitudes lissees
      DO 70 I=IT1,IT2
      FO=F(I)
      AM(I)=SMOO(KTT,F(IT1),C0,C1,C2,C3,FO)
      IF(AM(I).GT.AMAX) AM(I)=AMAX
      IF(AM(I).LT.A1) AM(I)=A1
   70 CONTINUE
      RETURN
      END

c-----------------------------------------------------------------------
      SUBROUTINE SZD(UM,F,DNU,DKM,IT1,KTT,CF,SF,TAU,NBAND)
c-----------------------------------------------------------------------
C     DISPERSION RESIDUELLE
c     entree: vit.gr. UM; sortie: tps de groupe lisse TAU (de IT1 a IT2)
      REAL UM(54),F(54),DU(54),TAU(54),OM(54),CF(*),SF(*)
      INTEGER NBAND(54)
      REAL C0(55),C1(55),C2(55),C3(55)
      REAL*8 R(440)

      PI2=8.*ATAN(1.)
      IT2=IT1+KTT-1
      DO 1 I=1,KTT
      OM(I)=F(I+IT1-1)*PI2
      TAU(I)=DKM/UM(I+IT1-1)
    1 DU(I)=0.01*DKM/(UM(I+IT1-1)**2)
      S=KTT
      SO=S/25.
      ILI=0
c     lissage du temps de groupe
   40 CALL LISSE(KTT,OM,TAU,DU,S,C0,C1,C2,C3,R)
      IF(C0(1).NE.0.) GOTO 30
      WRITE(*,*)'! DIFFICULTE DANS LISSAGE TEMPS DE GROUPE !'
      ILI=ILI+1
      S=SO*FLOAT(ILI)+FLOAT(KTT)-10.
      IF(ILI.LT.10) GOTO 40
      WRITE(*,*)'! IMPOSSIBLE A LEVER AVEC FLUCTUATIONS DE S !'
c     dispersion residuelle
   30 I1=OM(1)/DNU/PI2-NBAND(IT1)/2+1
      I2=OM(KTT)/DNU/PI2+NBAND(IT2)/2+2
      IF(I1.LT.1) I1=1
      CF(I1)=1.
      SF(I1)=0.
      I11=I1+1
      FO=0.
      DO 2 I=I11,I2
      P=FLOAT(I-1)*DNU*PI2
      T=OM(1)-P
      IF(T) 11,12,13
   11 DO 7 J=2,KTT
      TT=OM(J)-P
      IF(TT.GT.0.) GOTO 14
    7 CONTINUE
      K=KTT+1
      GOTO 15
   13 K=1
   15 C2(K)=0.
      C3(K)=0.
      GOTO 20
   12 K=2
      GOTO 20
   14 K=J
   20 KP=K-1
      IF(KP.EQ.0) KP=1
      D=P-OM(KP)
      FP=C0(K)*D+C1(K)*D*D/2.+C2(K)*D*D*D/3.+C3(K)*D*D*D*D/4.
      IF(K.EQ.1.OR.K.EQ.KTT+1) GOTO 16
      D=OM(K)-OM(KP)
      ZF=C0(K)*D+C1(K)*D*D/2.+C2(K)*D*D*D/3.+C3(K)*D*D*D*D/4.
      FP=FP-ZF
   16 IF(I.EQ.I11) KM=K
      IF(KM.EQ.K.OR.KM.EQ.KTT) GOTO 8
      L=KM
      NS=K-KM
      DO 10 J=1,NS
      LL=L
      L=L+1
      D=OM(L)-OM(LL)
   10 FO=FO+C0(L)*D+C1(L)*D*D/2.+C2(L)*D*D*D/3.+C3(L)*D*D*D*D/4.
      KM=L
    8 FT=FO+FP
      CF(I)=COS(FT)
    2 SF(I)=-SIN(FT)
c     restitution des temps de groupe lisses
      DO 90 I=1,IT2
      OMO=F(I)*PI2
      TAU(I)=SMOO(KTT,OM,C0,C1,C2,C3,OMO)
      IF(I.LT.IT1) TAU(I)=0.
   90 CONTINUE
      RETURN
      END

c-----------------------------------------------------------------------
      SUBROUTINE SZGRAPH(AMP,U,FRQ,KU,KT,IT1,KTT,UM,IBAND)
c-----------------------------------------------------------------------
c     dessin couleur et PostScript du diagramme + trace courbe de base
c     ou piquage souris pour recherche du max. UM.
c     IBAND=0 filtrage classique (+ courbe de base pour recherche max.).
c     IBAND=2 piquage souris sur meme graphe pour recherche max + precise.
c     IBAND=-1 dispersion residuelle + recherche max automatique.
c     IBAND=+1 dispers. res. & corr.amplitude + recherche max automatique.
c     IBAND=+/-3 filtrage variable/passe-bande sur dernier graphe.
!include "../convsuncore.h"
      COMMON/DIA/fu1,fku,h,x1,xkt,w
      REAL AMP(150,54),U(150),FRQ(54),UM(54)
      CHARACTER ligne*80

      if(IBAND.eq.2) goto 50

C  DESSIN DU DIAGRAMME
      if(IBAND.eq.0) then
       upas=U(2)-U(1)
       ul=upas/2.
       fu1=U(1)-ul
       fku=U(KU)+ul
       h=fku-fu1
       df=FRQ(2)-FRQ(1)
       fl=500.*df
       x1=-1000.*FRQ(KT)-fl
       xkt=-1000.*FRQ(1)+fl
       w=xkt-x1
       endif

c     ouverture,legendes...
      if(iabs(IBAND).le.1) then
       CALL OPENDIAG(U(1),U(KU),FRQ,KT)
       iseg=1
      else
       call DelRetainSegment(10)
       call DelRetainSegment(11)
       call DelRetainSegment(12)
      endif

c     legende type de diagramme
      iseg=iseg+1
      call CreateRetainSeg(iseg)
      call SetCharSize(w/100.,.016*h)
      call MoveAbs2(x1-.11*w,fu1-.08*h)
      if(IBAND.eq.0) call Text('premier filtrage classique')
      if(IBAND.eq.-1) call Text('dispersion residuelle')
      ligne(1:80)='disp. residuelle & correction d''amplitude'
      if(IBAND.eq.1) call Text(ligne(1:lnblnk(ligne)))
      call CloseRetainSeg(iseg)

c     recherche de l'amplitude max
      amax=0.
      do 25 i=1,kt
      do 25 j=1,ku
      if(AMP(j,i).gt.amax) amax=AMP(j,i)
   25 continue
      A1=0.434294
      AO=20.*A1*alog(amax)

c     coloriage diagramme
      iseg=iseg+1
      call CreateRetainSeg(iseg)
      do 30 k1=kt,1,-1
      do 30 k2=1,ku
      if(AMP(k2,k1).gt.0.) i=nint(AO-20.*A1*ALOG(AMP(k2,k1)))
      if(i.gt.48.or.AMP(k2,k1).le.0.) i=48
      i=max0(i,1)
      call SetFillIndex(i)
      gf=-1000.*FRQ(k1)
      CALL rectangle(gf-fl,gf+fl,U(k2)-ul,U(k2)+ul,0)
   30 continue
      call CloseRetainSeg(iseg)

c     legende filtrage
      if(IABS(IBAND).eq.3) then
       call CreateRetainSeg(10)
       call SetTextIndex(55)
       call MoveAbs2(xkt-.2*w,fu1-.08*h)
       if(IBAND.eq.3) call Text('+ filtrage variable')
       if(IBAND.eq.-3) call Text('+ filtrage passe-bande')
       call CloseRetainSeg(10)
      endif

C  TRACE DE LA COURBE DE VIT.DE GROUPE DE BASE
      KTT=0
      call CreateRetainSeg(11)
      call SetLineIndex(30)
      CALL TraceMoi(UM,FRQ,KT,IT1,IT2)

c     teste si courbe de base existe
      if(IT2.ne.0) then
c      complement legendes
       call SetTextIndex(30)
       call MoveAbs2(x1+.31*w,fu1-.08*h)
       call Text('+ vit. de groupe')
       endif
      call CloseRetainSeg(11)

      if(IT2.ne.0.and.IABS(IBAND).le.1) then
       call CreateRetainSeg(12)
       call SetTextIndex(53)
       call MoveAbs2(x1+.5*w,fu1-.08*h)
       call Text('+ recherche max. automatique')
       call CloseRetainSeg(12)
       endif
      if(IT2.ne.0) goto 90

C  POINTE D'APRES PIQUAGE SOURIS
c     initialisation piquage
   50 IBAND=2
      PRINT*,'Piquage souris sur bouton gauche; Sortie sur bouton droit'
      call DelRetainSegment(12)
      call CreateRetainSeg(12)
      call SetTextIndex(53)
      call MoveAbs2(x1+.5*w,fu1-.08*h)
      call Text('+ recherche max. par piquage souris')
      call CloseRetainSeg(12)
      IT1=0
      IT2=0
      uwm=0.
      xw0=x1
      xw=0.
c     test piquage d'un point a la souris
   55 call AwtButtonGetLoc2(20000000,1,ibout,xb,yb)
      if(ibout.eq.0) PRINT*,'Alors, ca vient ce Piquage ?'
      if(ibout.eq.0) goto 55
      if(ibout.ne.1.and.ibout.ne.3) PRINT*,'! Bouton incorrect !'
c     teste sortie possible
      if(ibout.eq.3.and.xw.eq.0.) PRINT*,'! Vous n''avez rien pique !'
      if(ibout.eq.3.and.xw.eq.0.) goto 55
      if(ibout.eq.3.and.xw.le.(x1+xkt)/2.) then
        PRINT*,'! Terminez de piquer plus a droite (2eme moitie) !'
        goto 55
        endif
      if(ibout.eq.3) goto 90
      if(ibout.ne.1) goto 55
c     point utile en coordonnees world
      call MapNdcToWorld2(xb,yb,xw,uw)
      if(xw.lt.x1.or.xw.gt.xkt.or.uw.lt.fu1.or.uw.gt.fku) then
        PRINT*,'! Piquage hors cadre (sortie sur bouton de droite) !'
        goto 55
        endif
      if(xw0.eq.x1.and.xw.ge.(x1+xkt)/2.) then
        PRINT*,'! Commencez a piquer plus a gauche (1ere moitie) !'
        goto 55
        endif
      if(xw.lt.xw0) PRINT*,'! Piquez de gauche a droite. Recommencez !'
      if(xw.lt.xw0) goto 55
c     piquage correct: dessin d'une croix
      xw0=xw
      ikt=1+int((xkt-xw)/2./fl)
      if(IT2.eq.0) IT2=ikt
      PRINT*,'Point pique: Periode',-1000./xw,'   Vit.Groupe',uw
      call CreateTempSeg()
      call SetLineIndex(40)
      call MoveAbs2(xw-fl/2.,uw)
      call LineRel2(fl,0.)
      call MoveAbs2(xw,uw-ul/2.)
      call LineRel2(0.,ul)
      call CloseTempSeg()
      if(IT1.eq.0) goto 60
c     interpolation lineaire entre 2 piquages successifs
      do 63 ik=IT1,ikt,-1
   63 UM(ik)=uwm+(uw-uwm)*float(ik-IT1)/float(ikt-IT1)
c     passage au point suivant
   60 IT1=ikt
      uwm=uw
      goto 55

C  SORTIE COURBE DE BASE POUR RECHERCHE DU MAX
c     sortie: KTT nombre de points entre le 1er et le dernier pique
   90 KTT=IT2-IT1+1
      RETURN
      END

c-----------------------------------------------------------------------
      SUBROUTINE TraceMoi(Y,FRQ,N,IT1,IT2)
c-----------------------------------------------------------------------
c     TRACE LA COURBE Y(FRQ) SUR LE GRAPHE.
c     RESSORT BORNES IT1 ET IT2 ENTRE LESQUELLES LA COURBE EST NON NULLE
!include "../convsuncore.h"
      REAL Y(*),FRQ(*)
      IT1=0
      IT2=0
      DO 10 I=1,N
      if(Y(i).eq.0.) goto 10
      IT2=i
      if(IT1.eq.0) call MoveAbs2(-1000.*FRQ(i),Y(i))
      if(IT1.eq.0) IT1=i
      call LineAbs2(-1000.*FRQ(i),Y(i))
   10 continue
      RETURN
      END

c-----------------------------------------------------------------------
      SUBROUTINE SZMAX(AMP,PHASE,U,FRQ,KU,IT1,KTT,TAU,DKM,DU,UM,AM)
c-----------------------------------------------------------------------
C     RECHERCHE MAX. D'ENERGIE AU VOISINAGE DE LA COURBE DE TEMPS DE
C     GROUPE TAU ENTREE. ZONE DE RECHERCHE EN VITESSE DE GROUPE:
C     U-DU a U+DU AVEC DU=0.4 Km/s AU PREMIER PASSAGE OU EN AUTOMATIQUE
C     ET 0.2 km/s ENSUITE OU 0.15 km/s AUTOUR PIQUAGE SOURIS. MODIFIABLE.
C     COURBES SORTIES LE LONG DU MAX: VIT.GR.UM(FRQ), AMPLITUDE AM(FRQ).
!include "../convsuncore.h"
      REAL AMP(150,54),PHASE(150,54),U(150),FRQ(54)
      REAL TAU(54),UM(54),AM(54),PM(54)
      parameter (iseg0=33)
      if(iseg.eq.0) iseg=iseg0-1

c     courbe TAU donnee entre indices IT1 et IT2
      IT2=IT1+KTT-1
      fl=250.*(FRQ(2)-FRQ(1))
      UPAS=U(2)-U(1)
      ul=UPAS/3.

  100 WRITE(*,'(A,F5.2)')'Recherche Max entre U-DU et U+DU pour DU=',DU
      iseg=iseg+1
      call CreateRetainSeg(iseg)

      DO 1 K2=IT2,IT1,-1
      UO=DKM/TAU(K2)
      UINF=UO-DU
      USUP=UO+DU
      KINF=1+INT((UINF-U(1))/UPAS)
      KSUP=KU-INT((U(KU)-USUP)/UPAS)
c     IF(KINF.LT.1.OR.KSUP.GT.KU) PRINT*,'BANDE DE RECHERCHE HORS CADRE'
      KINF=MAX0(KINF,1)
      KSUP=MIN0(KSUP,KU)
c     recherche de l'indice K1M du max entre KINF et KSUP
      AMAX=AMP(KINF,K2)
      K1M=KINF
      DO 5 K1=KINF,KSUP
      IF(AMP(K1,K2).LE.AMAX) GOTO 5
      AMAX=AMP(K1,K2)
      K1M=K1
    5 CONTINUE
      PM(K2)=PHASE(K1M,K2)
      IF(K1M.LE.KINF.OR.K1M.GE.KSUP) THEN
       WRITE(*,1000) 1./FRQ(K2)
 1000  FORMAT('MAX. ENERGIE EN BOUT INTERVALLE DE RECHERCHE A',F7.2,'s')
       AM(K2)=AMAX
       IF(K1M.LE.KINF) UM(K2)=UO-DU
       IF(K1M.GE.KSUP) UM(K2)=UO+DU
       call SetLineIndex(54)
       call SetFillIndex(5)
      else
       A1=AMP(K1M-1,K2)
       A2=AMP(K1M,K2)
       A3=AMP(K1M+1,K2)
       X1=DKM/U(K1M-1)
       X2=DKM/U(K1M)
       X3=DKM/U(K1M+1)
       Q1=A1/(X1-X3)/(X1-X2)
       Q2=A2/(X2-X3)/(X2-X1)
       Q3=A3/(X3-X2)/(X3-X1)
       XX=((X2+X3)*Q1+(X1+X3)*Q2+(X1+X2)*Q3)/(Q1+Q2+Q3)/2.
       AM(K2)=(XX-X3)*(XX-X2)*Q1+(XX-X3)*(XX-X1)*Q2+(XX-X2)*(XX-X1)*Q3
       UM(K2)=DKM/XX
       call SetLineIndex(53)
       call SetFillIndex(15)
      endif
c     trace du max (carres vert:ok, orange:bout) sur le diagramme
      gf=-1000.*FRQ(K2)
      CALL rectangle(gf-fl,gf+fl,UM(K2)-ul,UM(K2)+ul,1)
    1 CONTINUE

      DO 150 I=1,54
      IF(I.LT.IT1.OR.I.GT.IT2) UM(I)=0.
  150 CONTINUE
      PRINT*,'Courbe en vitesse de groupe apparente:',KTT,' points'
      write(*,'(18F5.2)') (UM(i),i=IT1,IT2)
      PRINT*,'Verifiez le bon pointe du max.: vert OK, rouge AU BOUT..'
      call CloseRetainSeg(iseg)

      PRINT*,'Clickez au CENTRE pour VALIDER, a GAUCHE pour chercher MOI
     &NS LARGE, a DROITE pour chercher PLUS LARGE'
  200 call AwaitAnyButton(2E7,ibout)
      if(ibout.eq.0) PRINT*,'Alors, ca vient ce Clickage ?'
      if(ibout.eq.0) goto 200
      if(ibout.eq.1) DU=DU-.1
      if(ibout.eq.3) DU=DU+.1
      IF(DU.LT.0.) PRINT*,'! On ne peut plus reduire la recherche !'
      IF(DU.GT.0.75) PRINT*,'! On ne peut augmenter indefiniment !'
      IF(DU.LT.0..OR.DU.GT.0.75) goto 250
      if(ibout.eq.1.or.ibout.eq.3) goto 100

  250 do 300 i=iseg0,iseg
  300 call DelRetainSegment(i)
      iseg=0
      RETURN
      END

c-----------------------------------------------------------------------
      subroutine OPENDIAG(UINF,USUP,FRQ,KT)
c-----------------------------------------------------------------------
      include "./usercore77.h"
!include "../convsuncore.h"
      character ligne*60,STA*4,CMP*1
      real FRQ(54)
      COMMON/STACMP/STA,CMP
      COMMON/DIA/fu1,fku,h,x1,xkt,w

c     ouverture fenetre
      if(iseg.eq.0) then
       iseg=1
c      on essaye en couleur sinon passera en gris
       icol=1
       CALL initdes(x1-.12*w,xkt+.08*w,fu1-.1*h,fku+.15*h,9,1,0,icol)
       call SetCharPrecision(CHARACTER)
      else
       call DelAllRetainSegs()
       call CreateRetainSeg(1)
      endif

c     cadre
      call SetFillIndex(51)
      CALL rectangle(x1-.11*w,xkt+.07*w,fu1-.05*h,fku+.1*h,0)
      call SetFont(0)
      call SetTextIndex(35)
      call SetCharSize(w/60.,.04*h)
      call MoveAbs2(x1+.02*w,fku+.05*h)
      ligne='Period-Group Velocity  ENERGY DIAGRAM - '//STA//' '//CMP
      call Text(ligne(1:lnblnk(ligne)))
      call SetFont(4)
 
c     legendes
      call SetTextIndex(38) 
      call SetCharSize(w/100.,.015*h) 
      do 20 U=UINF,USUP+h/50.,h/25.
      write(ligne,'(f5.2)') U
      call MoveAbs2(x1-.08*w,U) 
   20 call Text(ligne(1:5)) 
      call MoveAbs2(x1-.095*w,fku+.04*h) 
      call Text('U (km/s)') 
      ligne='     ' 
      do 21 i=1,KT,2 
      write(ligne,'(i4)') nint(1./FRQ(i)) 
      call MoveAbs2(-1000.*FRQ(i)-.7,fu1-.03*h) 
   21 call Text(ligne(1:4)) 
      call MoveAbs2(xkt+.01*w,fu1) 
      call Text('T (s)')
      call CloseRetainSeg(1)

      return
      end

C=======================================================================

      SUBROUTINE SZOND(FRQ,VGR,NFU,DF,DKM,T1,T2)

C=======================================================================
C     Calcul des limites en temps de groupe T1,T2 du filtrage variable.
C     Graphisme: dessine la vit.de groupe lissee & la bande du filtrage.
C     Entrees: VGR(FRQ) courbe de vitesse de groupe de base (NFU points)
C      DF pas en freq. de calcul de T1,T2 (1./duree du signal a filtrer)
C      DKM distance epicentrale en Kms
C     Sorties: T1 et T2 temps de groupe extremes tous les DF.
!include "../convsuncore.h"

      REAL FRQ(*),VGR(*),FQC(999),TG(999),T1(999),T2(999)

c     EL definit la marge de lissage (variable)
      EL=-0.1
c     ER definit la largeur du filtrage (de 8. a 12.)
      ER=10.

c     calcul des frequences par pas de DF
      NQ1=1+NINT(FRQ(1)/DF)
      NQ2=1+NINT(FRQ(NFU)/DF)
      IF(NQ2.GT.999) PRINT*,'! Signal trop long pour filtrage. On est ob
     &lige de tronquer le filtre a Haute-Frequence !'
      NQ2=MIN0(NQ2,999)
      N=NQ2-NQ1+1
      DO 10 I=1,999
   10 FQC(I)=0.
      DO 20 I=NQ1,NQ2
      J=I-NQ1+1
   20 FQC(J)=FLOAT(I-1)*DF

c     lissage moyen du temps de groupe par spline
   25 CALL INTERPOL(VGR,NFU,FRQ,EL,TG,FQC,N,2)

c     limites du filtrage variable
      DO 30 I=1,N
      FR=FQC(I)
      ST=DKM*TG(I)
c     largeur de filtrage analytique inspiree de Dziewonski (68)
      AL=AMAX1(ER/(FR**0.7),ER*10.)
      IF(ST.LT.AL) PRINT*,'! Signal debutant trop tard ou Filtre trop la
     &rge pour periode ',1./FR,'. Filtre tronque !'
      AL=AMIN1(ST,AL)
c     indices de debut et fin des ondelettes
      T1(I)=ST-AL
      T2(I)=ST+AL
   30 CONTINUE

c     dessin de la bande de filtrage
      CALL LIMFVR(T1,T2,DKM,FQC,N)

c     modification eventuelle du lissage
      PRINT*,'Clickez au CENTRE pour VALIDER, a GAUCHE pour MOINS LISSER
     &, a DROITE pour LISSER DAVANTAGE'  
   50 call AwaitAnyButton(2E7,ibout) 
      if(ibout.eq.0) PRINT*,'Alors, ca vient ce Clickage ?' 
      if(ibout.eq.0) goto 50
      if(ibout.eq.1) EL=EL/2. 
      if(ibout.eq.3) EL=EL*2.
      IF(EL.LT.-1.) PRINT*,'! On ne peut physiquement lisser plus !'
      IF(EL.GT.-0.01) PRINT*,'! Ca ne sert a rien de moins lisser !'
      IF(EL.LT.-0.01.AND.EL.GT.-1..AND.ibout.ne.2) goto 25

c     recalcul de vit.groupe lissee
      DO 60 I=1,NFU
      IT=1+NINT(FRQ(I)/DF)-(NQ1-1)
   60 VGR(I)=2.*DKM/(T1(IT)+T2(IT))

      RETURN
      END

c-----------------------------------------------------------------------
      subroutine LIMFVR(T1,T2,DKM,FRQ,KT)
c-----------------------------------------------------------------------
!include "../convsuncore.h"
      COMMON/DIA/fu1,fku,h,x1,xkt,w
      real frq(*),T1(*),T2(*)
      call DelRetainSegment(10)
      call CreateRetainSeg(10)
      call SetTextIndex(55)
      call MoveAbs2(xkt-.1*w,fu1-.08*h)
      call Text('+ limites filtrage')
      call SetLineIndex(55)
      do 1 i=2,kt
      tm=(T1(i)+T2(i))/2.
      tm1=(T1(i-1)+T2(i-1))/2.
      if(DKM/T1(I-1).lt.fu1.and.DKM/T1(I).lt.fu1) GOTO 2
      if(DKM/T1(I-1).gt.fku.and.DKM/T1(I).gt.fku) GOTO 2
      call MoveAbs2(-1000.*FRQ(i-1),DKM/T1(i-1))
      call LineAbs2(-1000.*FRQ(i),DKM/T1(i))
   2  if(DKM/T2(I-1).lt.fu1.and.DKM/T2(I).lt.fu1) GOTO 3
      if(DKM/T2(I-1).gt.fku.and.DKM/T2(I).gt.fku) GOTO 3
      call MoveAbs2(-1000.*FRQ(i-1),DKM/T2(i-1))
      call LineAbs2(-1000.*FRQ(i),DKM/T2(i))
   3  if(DKM/tm1.lt.fu1.and.DKM/tm.lt.fu1) GOTO 1
      if(DKM/tm1.gt.fku.and.DKM/tm.gt.fku) GOTO 1
      call MoveAbs2(-1000.*FRQ(i-1),DKM/tm1)
      call LineAbs2(-1000.*FRQ(i),DKM/tm)
   1  continue
      call CloseRetainSeg(10)
      return
      end

C=======================================================================

      SUBROUTINE NPAFNRJ(XC,LX,TAUO,DTES,F1,F2,T1,T2)

C=======================================================================
c     Filtrage d'un signal entre les courbes de temps de groupe T1 et T2
c     Entrees: XC spectre du signal + caracteristiques + heure origine + 
c            courbes limites T1 et T2 de F1 a F2 par pas de DF=1/LX/TAUO
c     Sortie:  XC= signal temps complexe filtre.

      COMPLEX XC(*),XI,XDT,CI
      REAL YF(20000),T1(999),T2(999)
      CI=CMPLX(0.,1.)
      PI=4.*ATAN(1.)
      DO 23 I=1,20000
  23  YF(I)=0.
      DF=1./LX/TAUO
      NF2=LX/2+1

c     calcul des indices limites des frequences
      NQ1=1+NINT(F1/DF)
      NQ2=1+NINT(F2/DF)
      IF(NQ2.GT.NF2) PRINT*,'! Attention, frequence Max > Nyquist !'
      NQ2=MIN0(NQ2,NF2)
      N=NQ2-NQ1+1
      IF(N.GT.999) PRINT*,'! Signal trop long pour filtrage. On est obli
     &ge de tronquer le filtre a Haute-Frequence !'
      N=MIN0(N,999)

c     TF-1 numerique incluant le produit par ondelettes 
      DO 26 J=NQ1,NQ1+N-1
      K=J-NQ1+1
      I1=INT(T1(K)/TAUO)+1
      I2=INT(T2(K)/TAUO)+2
      TIC=(T1(K)+T2(K))/2./TAUO+1.
      DT=(T2(K)-T1(K))/2.
      OM=2.*PI*FLOAT(J-1)/FLOAT(LX)
      XI=CEXP(CI*FLOAT(I1-2)*OM)
      XDT=CEXP(CI*OM*DTES/TAUO)
      DO 25 I=I1,I2
      DST=(FLOAT(I)-TIC)*TAUO
c     ap donne la forme de l'ondelette
      AP=2.*(1.-(DST/DT)**4)
      IF(DST.GE.DT) AP=0.
      XI=XI*CEXP(CI*OM)
c     TF-1: on se place / au t0 du seisme
  25  YF(I)=YF(I)+AP*DF*REAL(XC(J)*CONJG(XI)*XDT)
  26  CONTINUE

c     XC: signal temps complexe ramene au temps d'enregistrement
      NST=NINT(DTES/TAUO)
      IF(NST.GE.0) THEN
      DO 30 I=1,LX
  30  XC(I)=CMPLX(YF(I+NST),0.)
      ELSE
      DO 34 I=LX,1-NST,-1
  34  XC(I)=CMPLX(YF(I+NST),0.)
      DO 35 I=1,-NST
  35  XC(I)=CMPLX(0.,0.)
      ENDIF

      END

C==========================================================================

      SUBROUTINE VGVPH(CMP,U,FRQ,N,C)

C==========================================================================
c     Integre la vitesse de groupe U pour obtenir la vitesse de phase C.
c     Quelques valeurs de la vitesse de phase a longue-periode (100-600s)
c     sont tabulees pour Love et Rayleigh: cmp doit etre T pour Love, ou
c     Z,V ou L pour les Rayleigh; sinon demandera la valeur de C a Pmax.
c     Trace aussi la vitesse sur le diagramme.
!include "../convsuncore.h"
      COMMON/DIA/fu1,fku,h,x1,xkt,w

      character fich*40,cmp*1
      real u(54),c(54),frq(54),v(99),x(99),cl(6),cr(6)
c     valeurs des vitesses de phase initiales a 100,200,...600s
      data cr/4.06,4.55,5.24,5.90,6.38,6.65/
      data cl/4.63,4.90,5.17,5.57,5.83,6.01/

      pmax=amax1(1./frq(1),1./frq(n))
      fmin=1./pmax

c     valeur de la vitesse de phase a periode max
      print*,'CALCUL VITESSE DE PHASE: integration de la Vit.de groupe:'
   1  write(*,'(a,f5.0,a)')'Entrez fichier SAC de vitesse de phase pour 
     &caler la courbe a',pmax,'s'
      print*,'(par exemple Cprem.Rayl [R] ou Cprem.Love [L])'
      if(pmax.ge.100.and.pmax.le.600.and.(cmp.eq.'Z'.or.cmp.eq.'z'.or.cm
     &p.eq.'V'.or.cmp.eq.'L'.or.cmp.eq.'l'.or.cmp.eq.'T'.or.cmp.eq.'t'))
     &print*,'ou Return pour utiliser les valeurs par defaut.'
      READ(*,'(a)') fich
      if(fich.eq.' ') then
        i=min0(int(pmax/100.),5)
        if(cmp.eq.'T'.or.cmp.eq.'t') then
        cbeg=cl(i)+(cl(i+1)-cl(i))*(pmax-100.*i)/100.
        else
        cbeg=cr(i)+(cr(i+1)-cr(i))*(pmax-100.*i)/100.
        endif
      else
        if(fich.eq.'R') fich='Cprem.Rayl'
        if(fich.eq.'L') fich='Cprem.Love'
        call RSAC2(fich,v,nv,x,99,IERR)
        if(IERR.NE.0) goto 1
        CALL INTERPOL(v,nv,x,.001,cbeg,fmin,1,0)
      endif
      write(*,'(a,f5.0,a,f6.3)')'Vitesse de phase a',pmax,'s:',cbeg
      cte=1./pmax/cbeg

c     integration en frequence de 1/u
      CALL INTERPOL(u,n,frq,.005,c,frq,n,3)
c     calcul de la vitesse de phase integree
      do 9 i=1,n
    9 c(i)=frq(i)/(c(i)+cte)

c     trace sur le diagramme
      call CreateTempSeg()
      call SetLineIndex(53)
      ibeg=0
      ifin=0
      CALL TraceMoi(c,frq,n,ibeg,ifin)
      if(ifin.ne.0) then
      call SetTextIndex(53)
      call MoveAbs2(x1+.5*w,fu1-.08*h)
      call Text('+ vitesse de phase')
      endif
      call CloseTempSeg()

      return
      end

C==========================================================================

      SUBROUTINE INTERPOL(V,NF,X,DELTA,VOUT,FRQ,NFU,IOP)

C==========================================================================
C     Interpole la courbe V(X) (NF<=99 points). Determine si les X sont des
C     periodes ou frequences et passe ensuite en frequences F croissantes.
C     Sort la courbe VOUT aux NFU freq. FRQ par lissage d'ajustement DELTA:
C     ! Si DELTA < 0, prend un ajustement variable -SQRT(FRQ)*DELTA.
C     Si IOP=0 ressort le lissage de la fonction V (Ex: vitesse de groupe).
C     Si IOP=1 ressort la derivee de la fonction V (Ex: derivee de phase).
C     Si IOP=2 ressort le lissage de la fonction 1/V (Ex: delai de groupe).
C     Si IOP=3 ressort l'integrale de la fonction 1/V (Ex: pour Vit.phase).
C     Les courbes entrees ne sont pas modifiees.

      REAL V(*),X(*),F(99),Y(99),DY(99),A(100),B(100),C(100),D(100)
      REAL*8 R(800)
      REAL FRQ(*),VOUT(*)

c     IFP=1: X=frequences, IFP=2: X=periodes -> F=frequences croissantes
      IF(X(1).LE.1.) IFP=1
      IF(X(1).GT.1.) IFP=2
      IF(X(NF).LE.1.) IFP0=1
      IF(X(NF).GT.1.) IFP0=2
      IF(IFP.NE.IFP0) THEN
        PRINT*,'INTERPOL: tableau en frequences (1) ou periodes (2) ?'
        READ(*,*) IFP
        ENDIF
      IF(IFP.EQ.1) F0=X(1)
      IF(IFP.EQ.1) FN=X(NF)
      IF(IFP.EQ.2) F0=1./X(1)
      IF(IFP.EQ.2) FN=1./X(NF)
      IRENV=0
      IF(F0.GT.FN) IRENV=1

c     on se limite a la partie non nulle de la courbe
      I0=1
      J0=0
      DO 10 I=1,NF
      IF(IRENV.EQ.0) J=I
      IF(IRENV.EQ.1) J=NF-I+1
      IF(V(J).EQ.0.) GOTO 10
      JM=I
      IF(J0.EQ.0) J0=I
      IF(IFP.EQ.1) F(I)=X(J)
      IF(IFP.EQ.2) F(I)=1./X(J)
      IF(IOP.GE.2) THEN
        Y(I)=1./V(J)
        DY(I)=DELTA/(V(J)**2)
        ELSE
        Y(I)=V(J)
        DY(I)=DELTA
        ENDIF
c     si delta negatif, calcule un ajustement variable en racine(f)
      IF(DELTA.LT.0.) DY(I)=-DY(I)*SQRT(F(I))
c     on limite la courbe de depart autour de l'intervalle de sortie
      IF(F(I).LT.FRQ(1)) I0=I
      IF(F(I).LT.FRQ(NFU)) IM=I+1
   10 CONTINUE
      I0=MAX0(I0,J0,1)
      IM=MIN0(IM,JM,NF)
      N=IM-I0+1
      DO 20 I=1,N
      F(I)=F(I+I0-1)
      Y(I)=Y(I+I0-1)
   20 DY(I)=DY(I+I0-1)

      S0=FLOAT(N)/5.
      S=FLOAT(N)
      ITER=0
c     lissage de Y a DY pres sur N points -> coefficients A,B,C,D
   30 CALL LISSE(N,F,Y,DY,S,A,B,C,D,R)
      IF(A(1).EQ.0.) THEN
        PRINT*,'! Lissage difficile. Marge trop etroite ? Reessai... !'
        ITER=ITER+1
        S=S0*FLOAT(ITER)
        IF(ITER.LT.10) GOTO 30
        PRINT*,'! LISSAGE IMPOSSIBLE ! On sort comme on peut !'
        PRINT*,'! LISSEZ DAVANTAGE LA COURBE AVANT INTERPOLATION !'
        ENDIF

c     restitution aux frequences FRQ
      DO 50 I=1,NFU
      FQC=FRQ(I)
      IF(IOP.EQ.0.OR.IOP.EQ.2) THEN
      IF(FQC.LT.F(1)) PRINT*,'! Courbe de dispersion trop courte a Longu
     &e Periode ! Extrapolation B.F. (',1./fqc,'s) !'
      IF(FQC.GT.F(N)) PRINT*,'! Courbe de dispersion trop courte a Court
     &e Periode ! Extrapolation H.F. (',1./fqc,'s) !'
      VOUT(I)=SMOO(N,F,A,B,C,D,FQC)
      ENDIF
      IF(IOP.EQ.1) THEN
      IF(FQC.LT.F(1)) PRINT*,'! Extrapolation au 1er ordre de la derivee
     & a Longue Periode (',1./fqc,'s) !'
      IF(FQC.GT.F(N)) PRINT*,'! Extrapolation au 1er ordre de la derivee
     & a Courte Periode (',1./fqc,'s) !'
      VOUT(I)=SMDER(N,F,A,B,C,D,FQC)
      ENDIF
      IF(IOP.EQ.3) THEN
      IF(FQC.LT.F(1).OR.FQC.GT.F(N)) THEN
       PRINT*,'! Integration impossible hors intervalle a',1./fqc,'s.'
       PRINT*,'->  Vitesse de phase nulle !'
       VOUT(I)=0.
       ELSE
       VOUT(I)=SZSOML(N,F,A,B,C,D,FQC)
       ENDIF
      ENDIF
   50 CONTINUE

      RETURN
      END
C=======================================================================

      SUBROUTINE INSPZ(XC,LX,DF)
 
C=======================================================================
C   CALCULE LA FONCTION DE TRANSFERT XC A PARTIR DE SES POLES & ZEROS.
C   Pour un fichier de poles et zeros de type SAC.
C   RESSORT UNE FCT DE TRANSFERT / A L'ACCELERATION DU SOL EN M.S-2 SI
C   REPONSE EN DIGITS/M(/S(/S)).

      CHARACTER DVA*6
      CHARACTER TYPE*1
      COMPLEX XC(*),POLE(30),ZERO(30),S,Z
      COMMON/PZ/NZE,ZERO,NPO,POLE,A0
      PII=6.2831852
      NF2=LX/2+1

      A0=0.
      DVA='      '

      CALL DATIPG(DVA)

      IF(DVA(1:2).EQ.'M ') TYPE='D'
      IF(DVA(1:3).EQ.'M/S'.OR.DVA(1:4).EQ.'MS-1') TYPE='V'
      IF(DVA(1:6).EQ.'M/S**2'.OR.DVA(1:5).EQ.'M/S/S'.OR.DVA(1:4).EQ.
     &'MS-2') TYPE='A'
      IF(TYPE.NE.'A'.AND.TYPE.NE.'V'.AND.TYPE.NE.'D') THEN
        WRITE(6,*)'! TYPE DE FONCTION (/D,V,A du sol) INCONNU !'
        IF(NZE.EQ.1) TYPE='A'
        IF(NZE.EQ.2) TYPE='V'
        IF(NZE.EQ.3) TYPE='D'
        IF(NZE.GE.1.AND.NZE.LE.3) THEN
          WRITE(6,*)'Le nombre de zeros semble indiquer qu''il s''agit d
     *u type',TYPE,'Confirmez: D,V ou A'
          ELSE
          WRITE(6,*)'Je ne peux deduire le type du nombre de zeros. Entr
     *ez le type: D,V ou A'
          ENDIF
        READ(5,'(A1)') TYPE
        ENDIF
      WRITE(6,*)'Type fonction de transfert:',TYPE
 
c     convention de TF avec + (Aki)
      DO 100 I=2,NF2
        S=CMPLX(0.,PII*FLOAT(I-1)*DF)
c       correction de gain en cm
        XC(I)=CMPLX(A0,0.)
c       fct transf en deplacement
        DO 10 J=1,NZE
  10    XC(I)=XC(I)*(S-ZERO(J))
        DO 20 J=1,NPO
  20    XC(I)=XC(I)/(S-POLE(J))
        XC(I)=CONJG(XC(I))
        Z=XC(I)
c       fct transf en vitesse > fct en acceleration
        IF(TYPE.EQ.'V') XC(I)=-XC(I)/S
c       fct transf en deplacement > fct en acceleration
        IF(TYPE.EQ.'D') XC(I)=XC(I)/S**2
 100  CONTINUE

      DO 300 I=NF2+1,LX
 300  XC(I)=(0.,0.)

      RETURN
      END


C=======================================================================
      SUBROUTINE DATIPG(DVA)
 
C   LECTURE (UNITE 20) D'UN FICHIER INDIVIDUEL DE POLES ET ZEROS POUR
C   UNE COMPOSANTE DONNEE, AU FORMAT TYPE IPG (COMPATIBLE SAC).
C   FORMAT TYPE: zeros    NZE
C                Re(zeros) Im(zeros) NZEx
C                poles    NPO
C                Re(poles) Im(poles) NPOx
C                constant A0 et unite
C   SI LE TYPE DE FCT DE TRANSFERT (DONNE PAR L'UNITE DE A0) N'APPARAIT
C   PAS, ON LE DEMANDE (D,V ou A) OU ON LE DEDUIT DU NOMBRE DE ZEROS NULS
C   ( 1-> A, 2-> V, 3-> D ).
C   SI TOUTE LA DERNIERE LIGNE N'APPARAIT PAS, ON PREND 1 DU/M(/S(/S))
 
      CHARACTER FICH*40,LIGNE*8,DVA1*40,DVA*6
      COMPLEX POLE(30),ZERO(30)
      COMMON/PZ/NZE,ZERO,NPO,POLE,A0
      DATA IU/20/
 
      WRITE(6,*)'Nom du Fichier Poles/Zeros type SAC ?'
      READ(5,'(A40)') FICH
  5   OPEN(IU,FILE=FICH,STATUS='OLD',IOSTAT=ICO)
      IF(ICO.NE.0) THEN
        WRITE(6,*)'Fichier ',FICH(:LNBLNK(FICH)),' non trouve. Autre ?'
        READ(5,'(A40)') FICH
        GOTO 5
        ENDIF
      REWIND(IU)

      READ(IU,*) LIGNE,NZE
      IF(LIGNE(1:5).NE.'ZEROS'.AND.LIGNE(1:5).NE.'zeros')
     *WRITE(6,*)'Erreur format sur ligne des Zeros'
      DO 17 I=1,NZE
      READ(IU,*) ZR,ZI
  17  ZERO(I)=CMPLX(ZR,ZI)

      READ(IU,*) LIGNE,NPO
      IF(LIGNE(1:5).NE.'POLES'.AND.LIGNE(1:5).NE.'poles')
     *WRITE(6,*)'Erreur format sur ligne des Poles'
      DO 18 I=1,NPO
      READ(IU,*) PR,PI
  18  POLE(I)=CMPLX(PR,PI)
      READ(IU,'(A8,A)',END=999) LIGNE,DVA1
      IF(LIGNE.NE.'CONSTANT'.AND.LIGNE.NE.'constant')
     *WRITE(6,*)'Erreur format sur ligne de la Constante'
      READ(DVA1,*) A0
      DO 20 I=1,LNBLNK(DVA1)-2
      IF(DVA1(I:I+2).EQ.'DU/'.OR.DVA1(I:I+2).EQ.'du/') GOTO 21
  20  CONTINUE
  21  DVA=DVA1(I+3:LNBLNK(DVA1))
      WRITE(6,*)'Sensibilite globale:',A0,'units/'//DVA

  999 CLOSE(IU)
      RETURN
      END

