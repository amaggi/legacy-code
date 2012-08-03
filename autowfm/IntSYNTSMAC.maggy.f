      parameter(NPROFMAX=50,NBTRAMAX=6000,NSTAMAX=200)
      parameter(NPERTMAX=50,NPTSMAX=8000,NCOUPMAX=20)
      real slat,slon,elat,elon,depth,sig,vscoup,profmin,profmax

      real v(NPROFMAX,NPTSMAX),c(NPROFMAX,90,180)
      real a1(NPROFMAX,90,180),a2(NPROFMAX,90,180)
      real cinlatc(NPROFMAX,NCOUPMAX,180),cingc(NPROFMAX,NPTSMAX)
      real vout(NPROFMAX),prof(NPROFMAX),vsout(NPROFMAX)
      real coupin(NPROFMAX),dcoupin(NPROFMAX)
      real lentsom(NPROFMAX),vsin(NPROFMAX),dvsin(NPROFMAX)
      real voutsom(NPROFMAX)
      real aa(NPROFMAX),bb(NPROFMAX)
      real cc(NPROFMAX),dd(NPROFMAX)
      real cref(NPROFMAX),crefi(NPROFMAX),dcrefi(NPROFMAX)
      real crefpass,crefbid,crefint(NPROFMAX)
      real*8 rr(8*(NPROFMAX+1))

      real latc(NCOUPMAX),lat1(NCOUPMAX),lon1(NCOUPMAX)
      real lat2(NCOUPMAX),lon2(NCOUPMAX)
      real latmod(NPERTMAX),lonmod(NPERTMAX),depthmod(NPERTMAX) 
      real late(NBTRAMAX),lone(NBTRAMAX) 
      real lats(NSTAMAX),lons(NSTAMAX) 
      real err(NBTRAMAX,NPROFMAX) 
      real zsmac(23),zone(4),len(NPTSMAX)
      data zsmac/50,60,70,80,90,100,115,130,150,175,200,
     *250,300,350,390,430,460,510,530,560,610,640,680/
      integer icout(NCOUPMAX),icout2(NCOUPMAX)
      character*85 toto
      character*4 stat(NSTAMAX),sta2
      character*2 nocoup
      character*10 sta(NBTRAMAX),stabid
      character*80 nomtra(NBTRAMAX),filein,filecout
      character*80 nomtrajet,filecout2

      common/c1/nbtra,nsta,ncouch

c----------------------------------------------------------
c    14 JUIN 2002 : version adaptee pour effectuer un test synthetique
c    pour Maggy    
c    25 JUIN 2001 : CODE ADAPTE A PARTIR DE IntsurtrajetsVs1anis.f
c    en vue du calcul d'un test synthetique complet
c    incluant waveform modelling et inversion 3D.
c
c    Programme de calcul de models synthetiques dans 3SMAC
c    Il calcul l'integrale des lenteurs sur un trajet donne.
c    les coordonnees des noms de trajets  et des epicentres
c    et stations sont dans un fichier format Montagner.
c
c    Cette version lit 
c        -3SMAC 
c        - un fichier intomodesVs 
c          pour la couverture en trajet que l'on veut mettre dans le test.
c    Cette version permet 
c        - d'ajouter une anisotropie synthetique a 3SMAC.
c        - d'ajouter un dVs a 3SMAC a une longitude et latitude donnee, etalees
c          sur 2 bloces en longitude et un en latitude
c    
c    Les vitesses integrees le long des trajets dans 3SMAC sont interpollees
c    tout les 25 km->700km.
c    Cette version sort en plus une coupe a lat=17 deg du model 3SMAC +dVS
c----------------------------------------------------------

      ncouch=32
      nvar=2
      isor=12
      pi=4.*atan(1.)
      degtorad=pi/180.0
      radtodeg=180.0/pi      

      write(*,*)'Lecture du fichier intomodesVs'
      call lecintomo(sta,stat,lats,lons,nomtra,late,lone,err)

      print *,'entrer latsup,latinf,lonouest,lonest'
      print *,'pour le plot de la carte GMT'
      print *,'(51 -41 -3 121)'
      print *,'(46 -71 -128 51) pour l\'Amerique du Sud'
      read(*,*)zone(1),zone(2),zone(3),zone(4)
      write(*,*)zone(1),zone(2),zone(3),zone(4)

      print *,'Entrez corr, facteur de correction pour le plot de'
      print *,'l\'anisotropie(corr=5->1cm = 5% d\'aniso sur les cartes)'
      read(*,*)corr

c  Lecture de 3SMAC pour toute les profondeurs comprises
c  entre 50 km et 660km (couches 12 a 34).
c     
      write(*,*)'Entrer 1 pour ajouter une perturbation de Vs a 3SMAC'
      read(*,*)ipert
      if(ipert.eq.1)then
        write(*,*)'Lat Lon(0->359) des 2 perturbations?'
        write(*,*)'(e.g. 16 40 16 44) '
        write(*,*)'(e.g. -20 292 -24 292) South America'
        read(*,*)latpert1,lonpert1,latpert2,lonpert2
        write(*,*)'Prof min et max de la perturbation? '
        write(*,*)'(e.g. 300 460) '
        write(*,*)'(e.g. 50 300)  South America'
        read(*,*)profmin,profmax
      endif

       print *,'Number of cross-sections within the starting model '
       print *,'and lat of the cross-sections '
       write(*,*)'(e.g. 2,15,17)'
       read(*,*)nbcoup,(latc(ij),ij=1,nbcoup)
       write(*,*)'Also', nbcoup,' cross section along great circles'
       do 4 i=1,nbcoup
       print *,'Coord of 2 points to define the minor arc along which'
       write(*,*)'the cross-section ',i,' is plotted '
       write(*,*)'(lon1,lat,lon2,lat2)'
       write(*,*)'(e.g. 25,-15,52,25)'
       write(*,*)'(e.g. 290,-30,310,10)'
       read(*,*)lon1(i),lat1(i),lon2(i),lat2(i)
4      write(*,*)lon1(i),lat1(i),lon2(i),lat2(i)

       open(9,file='gc.cross.xy')
       write(9,'(1a)')'>'
       do 5 jc=1,nbcoup
         write(9,*)lat1(jc),lon1(jc)
         write(9,*)lat2(jc),lon2(jc)
         write(9,'(1a)')'>'
        icout(jc)=24+jc
        icout2(jc)=44+jc
        write(nocoup,'(i2.2)')jc
        filecout= 'VSSMACcoupgc'//nocoup//'.xyz'
        filecout2= 'VSSMACcouplat'//nocoup//'.xyz'
        open(icout(jc),file=filecout)
        open(icout2(jc),file=filecout2)
5      continue
       close(9)

c------------------------------------------------------------------

      write(*,*)'Lecture de 3SMAC'
      write(*,*)'NB : les discontinuites sont lissees'
      ipass=0 
      comptpert=0
      do 10 jp=12,34
       write(filein,'("VS.",i2)')jp
       call lecdata(c,filein,jp,zsmac)

c---AJOUT DE PERTURBATIONS ISOTROPES ET/OU ANISOTROPES DANS 3SMAC---
c  a1 a2 pour avoir 2% d'anisotropie avec c0 =4.5 km/s
      cref(jp)=0
      do 6 i=1,90
      do 6 j=1,180
        lat= 90-(i+(i-1))
        lon=j+(j-1)
        cref(jp)=cref(jp)+c(jp,i,j)
      if((ipert.eq.1.).and.
     *(((abs(lat-latpert1).le.1).and.(abs(lon-lonpert1).le.1)).or.
     *((abs(lat-latpert2).le.1).and.(abs(lon-lonpert2).le.1))).and.
     *(zsmac(jp-11).le.profmax.and.zsmac(jp-11).ge.profmin)) then
       ilatcoup=i
       write(*,*)'Lat',lat,' Lon',lon,' CP ',c(jp,i,j),' Z',zsmac(jp-11) 
       c(jp,i,j)=c(jp,i,j)-0.05*c(jp,i,j)
       if(ipass.eq.0) then
         comptpert=comptpert+1
         latmod(comptpert)=lat
         lonmod(comptpert)=lon
         depthmod(comptpert)=zsmac(jp-11)
       endif
      endif
c          a1(jp,i,j)= 0.   
c          a2(jp,i,j)= 0.045
c       if(  lon.le.40.or.
c    *       lon.gt.85) then
c    *      (lon.gt.60.and.lon.le.80).or.
c    *      (lon.gt.90.and.lon.le.100)   )then
c   fast axes a 45N, 2% d'aniso crete crete
c          a1(jp,i,j)= 0.   
c          a2(jp,i,j)= 0.045
c       else
c   fast axes a 135N, 2% d'aniso crete crete
c          a1(jp,i,j)= 0.   
c          a2(jp,i,j)= -0.045
c       endif
c   fast axie a horizontal 2% d'aniso crete crete
c      a1(jp,i,j)= -0.045   
c      a2(jp,i,j)= 0.0
c   fast axes a 22.5N, 2% d'aniso crete crete
c      a1(jp,i,j)= 0.031819805
c      a2(jp,i,j)= 0.031819805
       a1(jp,i,j)= 0.0   
       a2(jp,i,j)= 0.0
6     continue
      cref(jp)=cref(jp)/90/180
c
c---Ecriture du model en INPUT 3SMAC en format xyz--------------------
c   Ecritures des cartes a differentes profondeurs
      crefpass=cref(jp)
      call writedata(c,a1,a2,filein,jp,corr,zone,zsmac,crefpass)
      if(comptpert.ne.0)ipass=1
      crefi(jp-11)=cref(jp)*100.
      dcrefi(jp-11)=crefi(jp-11)/10000.

10    continue

         sig=23.
         call lisse(23,zsmac,crefi,dcrefi,sig,aa,bb,cc,dd,rr)
         prof(1)=75
c        On interpolle jusqu'a 675 km de profondeur
c        car 3SMAC va jusqu'a 680 km
         write(*,*)''
         write(*,*)'cref interpolle aux prof. des coupes'
         do 16 ii=1,25
         if(ii.ge.2)prof(ii)=prof(ii-1)+25
         if(prof(ii).gt.700)prof(ii)=prof(ii-1)+50
         depth=prof(ii)
         crefbid=smoo(23,zsmac,aa,bb,cc,dd,depth)
         crefint(ii)=crefbid/100.
         write(*,*)prof(ii),crefint(ii)
16    continue

c--------COUPES DANS LE MODELE INITIAL-----------------------------

      do 200 jc=1,nbcoup
       alat1=lat1(jc)
       alon1=lon1(jc)
       alat2=lat2(jc)
       alon2=lon2(jc)
       call mygrt(alat1,alon1,alat2,alon2,disd,az12,az21,gc)
       delta=disd
       azes=az12
       if(azes.lt.0.) azes=azes+360.
25     if(azes.le.360.) go to 30 
       azes=azes-360.
       go to 25 
30     continue
       npas= int(5.*delta+1.)
       dpas= delta/float(npas)
       npat= npas+1
       if(npat.ge.NPTSMAX) then
         write(*,*) 'NPAT ',npat
       endif

c       Boucle sur chaque point de l'azimuth epicentre-station
        do 32  jcle=1,npat
          dstcle=dpas*float(jcle-1)*111.195
          alat1=lat1(jc)
          alon1=lon1(jc)
          call mygds(alat1,alon1,azes,dstcle,blat,blon)
          i=int(((90-blat+1)/2)+0.5)
          if(blon.lt.0.)blon=360.+blon
          j=int(((blon+1)/2)+0.5)
          do 35 jp=12,34
           jt=jp-11
           coupin(jt)=c(jp,i,j)*100.
35         dcoupin(jt)=coupin(jt)/10000.
         sig=23.
         call lisse(23,zsmac,coupin,dcoupin,sig,aa,bb,cc,dd,rr)
         prof(1)=75
c        On interpolle jusqu'a 675 km de profondeur
c        car 3SMAC va jusqu'a 680 km
         do 40 ii=1,25
         if(ii.ge.2)prof(ii)=prof(ii-1)+25
         if(prof(ii).gt.700)prof(ii)=prof(ii-1)+50
         depth=prof(ii)
         vscoup=smoo(23,zsmac,aa,bb,cc,dd,depth)
         vscoup=vscoup/100.
        cingc(ii,jcle)=100*(vscoup-crefint(ii))/crefint(ii)
40    continue
32    continue

      do 45  ii=1,25 
      do 45  jcle=1,npat 
       write(icout(jc),1200)-prof(ii),jcle,vscoup,cingc(ii,jcle)
45    continue

200   continue

c  Coupes a latitudes constantes

      do 300 ic=1,nbcoup
      ilatc=int(latc(ic))
      ilatc=int(((90-ilatc+1)/2)+0.5)

      do 50 jlon=1,180
          do 60 jp=12,34
           jt=jp-11
           coupin(jt)=c(jp,ilatc,jlon)*100.
           dcoupin(jt)=coupin(jt)/10000.
60    continue
         sig=23.
         call lisse(23,zsmac,coupin,dcoupin,sig,aa,bb,cc,dd,rr)
         prof(1)=75
c        On interpolle jusqu'a 675 km de profondeur
c        car 3SMAC va jusqu'a 680 km
         do 70 ii=1,25
         if(ii.ge.2)prof(ii)=prof(ii-1)+25
         if(prof(ii).gt.700)prof(ii)=prof(ii-1)+50
         depth=prof(ii)
         vscoup=smoo(23,zsmac,aa,bb,cc,dd,depth)
         vscoup=vscoup/100.
        cinlatc(ii,ilatc,jlon)=100*(vscoup-crefint(ii))/crefint(ii)
70    continue
50    continue

      do 80 ii=1,25
      do 80 jlon=1,180
          lon=jlon+(jlon-1)
          write(icout2(ic),1200)-prof(ii),lon,cinlatc(ii,ilatc,jlon)
80     continue
300    continue


1200  format(f8.2,2x,i5,2x,f8.3)


c    On ferme les fichiers contenant les coupes
      do 90 jc=1,nbcoup
       close(icout(jc))
       close(icout2(jc))
90    continue

c   Ecris les lat lon prof ou dVs a ete apporte
      write(*,*)'dVs SMAC apporte a'
      do 100 ij=1,comptpert
      write(*,*)'lat ',latmod(ij),' lon ',lonmod(ij),
     *' depth ',depthmod(ij)
100   continue

      write(*,*)''
      write(*,*)'Initialisation du model de depart achevee'
      write(*,*)''

c--------FIN DES COUPES-----------------------------

       open(isor,status='new',file='intomo3SMAC')
       write(isor,*)nsta
       do 120 j=1,nsta
       write(isor,1001)stat(j),lats(j),lons(j)
120    continue
      
      do 125 jp=1,NPROFMAX
       voutsom(jp)=0.
125    lentsom(jp)=0.

      do 130 itra=1,nbtra
      nomtrajet=nomtra(itra)
      istart=index(nomtrajet,'19')
      is2=index(nomtrajet,'2000')
      if(is2.lt.istart.and.is2.ge.1)istart=is2
      if(istart.eq.0)istart=is2
      nomtrajet=nomtrajet(istart:lnblnk(nomtrajet))
      elat=late(itra)
      elon=lone(itra)
      stabid=sta(itra)
      sta2=stabid(1:4)

      if (itra.lt.10) then
      write(isor,1004) sta2,".z0",itra,nomtrajet
      elseif (itra.lt.100) then
      write(isor,1005) sta2,".z",itra,nomtrajet
      elseif (itra.lt.1000) then
      write(isor,1006) sta2,".z",itra,nomtrajet
      elseif (itra.lt.10000) then
      write(isor,1007) sta2,".z",itra,nomtrajet
      endif
      write(isor,*)elat,elon

1004  format(a,a,i1,5x,a)
1005  format(a,a,i2,5x,a)
1006  format(a,a,i3,5x,a)
1007  format(a,a,i4,5x,a)

      call findcoord(stat,sta2,slat,slon,lats,lons)
      if(itra.eq.802) write(*,*)stabid,slat,slon

	call mygrt(elat,elon,slat,slon,disd,az12,az21,gc)
	delta=disd
	azes=az12
      if(azes.lt.0.) azes=azes+360.
140    if(azes.lt.360.) go to 150   
	azes=azes-360.
	go to 140   
150  	continue
c     npas= int(5.*delta+1.)
      npas= int(20.*delta+1.)
	dpas= delta/float(npas)
	npat= npas+1
c	write(*,*) dpas, npas, npat

c  Boucle sur chaque point de l'azimuth epicentre-station

      do 160 jcle=1,npat 
        dstcle=dpas*float(jcle-1)*111.195 
	call mygds(elat,elon,azes,dstcle,blat,blon)
        call mygrt(elat,elon,blat,blon,dstcle,azim,azback,gc)
        i=int(((90-blat+1)/2)+0.5)
        if(blon.lt.0.)blon=360.+blon
        j=int(((blon+1)/2)+0.5)

        do 170 jp=12,34 
          v(jp,jcle)=c(jp,i,j)
          v(jp,jcle)=v(jp,jcle)+a1(jp,i,j)*cos(2*degtorad*azim)+
     *a2(jp,i,j)*sin(2*degtorad*azim)
          if(v(jp,jcle).eq.0) then
            write(*,*)'V NUL TRA',itra
            stop
          endif
170  	continue 
 
160  	continue 

c  Il ne reste plus qu'a integrer sur le trajet.
    
      do 180 jp=12,34
        do 190 jcle=1,npat 
           len(jcle)=1./v(jp,jcle)
190      continue
        call trapez(npat,len,dpas,sint)
        if(delta.eq.0) then 
             write(*,*)'DELTA NUL TRA',itra
             stop
        endif
        vout(jp)=sint/delta
        if(vout(jp).eq.0)then
             write(*,*)'VOUT NUL TRA',itra
             stop
        endif
        vout(jp)=1./vout(jp)
c Modif eric pour MAGGY 14/06/2002
        voutsom(jp)=voutsom(jp)+vout(jp)
180    continue

c -Interpolation des modeles 1D tout les 25 km 
      do 210 ii=12,34
      jj=ii-11
      vsin(jj)=vout(ii)*100.
      dvsin(jj)=vsin(jj)/10000.
210   continue

      sig=23. 
      call lisse(23,zsmac,vsin,dvsin,sig,aa,bb,cc,dd,rr)
      
c     toto='des.S'//nomtrajet
c     write(*,'(a)')toto 
c     open(9,file=toto,status='new') 
c     write(9,'(a)') 'Synthetic within 3SMAC + dvs'
c     write(9,*)' 1'
c     write(9,*)' 25 ',nvar

      do 220 ii=1,25
        depth=prof(ii)
        vsout(ii)=smoo(23,zsmac,aa,bb,cc,dd,depth)
c       write(9,*)depth,vsout(ii)/100., 0.0050
        lentsom(ii)=lentsom(ii)+(100./vsout(ii))
 220   continue
      close(9)


c       write(isor,1003) (vsout(ii)/100.,ii=1,25)
c       write(isor,1003) (err(itra,ii),ii=1,25)

c modif pour MAGGY  on ecrit dans IntSYNTSMAC.f VS 3SMAC a 50km a toutes les profondeurs
c (ca sera le modele initial a toutes les profondeurs).
c par contre, pour l'erreur a posteriori, on met celle des donnees reeles a chaque profondeurs
c depuis 75 km -> la 25eme couche.

        write(isor,1003) (vout(12),ii=1,25)
        write(isor,1003) (err(itra,ii),ii=1,25)

130    continue

c  Affichage de la moyenne des lenteurs a 100 km utiliser comme modele
c  a priori

c modif affichage pour Maggy 14/06/2002     
        write(*,*)' VS a priori pour MAGGY issue de 3SMAC 50 km (VS.12)'
      do 230 ii=1,25
c       write(*,*) prof(ii),1./(lentsom(ii)/nbtra),' .05 0. 0. 0.005'
        write(*,*) prof(ii),(voutsom(12)/nbtra),' .05 0. 0. 0.005'
 230   continue

 
      write(*,*) 'ATTENTION : toutes les couches Vs de 3SMAC sont'
      write(*,*) 'prises en compte dans intomo3SMAC. Intsurtrajetanis.f'
      write(*,*) 'permet de prendre en compte une seule couche '
      
      close (isor)
1001  format(a4,2x,2f9.4)
1003  format(50f8.4)
      end

c @@@@@@@@@@@@@@@@@@ Liste des Subroutines @@@@@@@@@@@@@@@@@@@@
      SUBROUTINE MYGRT(ALAT1,ALON1,ALAT2,ALON2,DISD,AZ12,AZ21,GC)

c Great circle program. Given epicenter and station, finds distance,
c take-off and back- azimuths, and length of great circle

      PI = 4.*atan(1.0)
      ATH=6378.388
      BTH=6356.912
      RAD = PI/180.
      H = 1. - BTH*BTH/(ATH*ATH)
      P = H/(1. - H)
      GR = ALON1*RAD
      TR = ALAT1*RAD
      SINTR = SIN(TR)
      COSTR = COS(TR)
      IF (SINTR .EQ. 0.) SINTR = .000001
      IF (COSTR .EQ. 0.) COSTR = .000001
      R1 = ATH/SQRT(1. - H*SINTR*SINTR)
      Z1 = R1*(1. - H)*SINTR
      G = ALON2*RAD
      T = ALAT2*RAD
      IF (T .EQ. 0.) T = .00001
      SINT = SIN(T)
      COST = COS(T)
      R2 = ATH/SQRT(1. - H*SINT*SINT)
      DG = G - GR
      COSDG = COS(DG)
      COSDG = COS(DG)
      SINDG = SIN(DG)
      DGR = GR - G
      DT = T - TR
      Q = SINT*COSTR/((1. + P)*COST*SINTR) + H*R1*COSTR/(R2*COST)
      X = R2*COST*COSDG
      Y = R2*COST*SINDG
      Z = R2*(1. - H)*SINT
      AZ12 = ATAN2(SINDG,(Q - COSDG)*SINTR)
      Q = SINTR*COST/(COSTR*SINT*(1. + P)) + H*R2*COST/(R1*COSTR)
      AZ21 = ATAN2(SIN(DGR),SINT*(Q - COS(DGR)))
      COS12 = COS(AZ12)
      CTA2 = COSTR*COSTR*COS12*COS12
      P0 = P*(CTA2 + SINTR*SINTR)
      B0 = (R1/(1. + P0))*SQRT(1. + P*CTA2)
      E0 = P0/(1. + P0)
      GC = 2.*PI*B0*SQRT(1. + P0)*(1. - E0*(.25 + E0*(3./64.
     *                                          + 5.*E0/256.)))
      C0 = 1. + P0*(.25 - P0*(3./64 - 5.*P0/256.))
      C2 = P0*(-.125 + P0*(1./32. - 15.*P0/1024.))
      C4 = (-1./256. + 3.*P0/1024.)*P0*P0
      U0 = ATAN2(SINTR,COSTR*COS12*SQRT(1. + P0))
      U = ATAN2(R1*SINTR + (1. + P0)*(Z - Z1),(X*COS12 - Y*SINTR*
     *                                         SIN(AZ12))*SQRT(1. + P0))
      DISD = U - U0
      IF (U .LT. U0) DISD = PI + PI + DISD
      DIST = B0*(C0*( DISD ) + C2*(SIN(U + U) - SIN(U0 + U0))
     *                       + C4*(SIN(4.*U) - SIN(4.*U0)))
      DISD = DISD/RAD
      AZ12 = AZ12/RAD
      AZ21 = AZ21/RAD
      IF (AZ12 .LT. 0.) AZ12 = 360. + AZ12
      IF (AZ21 .LT. 0.) AZ21 = 360. + AZ21
        if(disd.gt.355.) disd=360.-disd
      RETURN
      END

c-------------------------------------------------------------

      SUBROUTINE MYGDS(ALAT,ALON,AZ,DIST,BLAT,BLON)

c Geodesic program. Given epicenter (alat,alon), azimuth of travel and distance
c traveled, finds arrival point (blat,blon);  dist in km

      ALT=ALAT
      ALN=ALON
      IF(ALAT.EQ.90.)ALAT=ALAT-0.000001
      IF(AZ.EQ.90.0.OR.AZ.EQ.270.0)AZ=AZ-0.000001
      A0=6378.388
      B0=6356.912
      F=(A0-B0)/A0
      PI=4.*ATAN(1.0)
      RAD=PI/180.0
      ALAT=ALAT*RAD
      ALON=ALON*RAD
      AZ=AZ*RAD
      E2=1.-(B0**2)/(A0**2)
      EPS=E2/(1.-E2)
      SIN1=SIN(ALAT)
      V=A0/(SQRT(1.-E2*(SIN1**2)))
      SINA=SIN(AZ)
      COSA=COS(AZ)
      COS1=COS(ALAT)
      TC2=(COSA*COSA)*(COS1*COS1)
      C2=TC2+(SIN1*SIN1)
      EPS0=C2*EPS
      TB1=SQRT(1.+EPS*TC2)
      B1=(V*TB1)/(1.+EPS0)
      TAN1=TAN(ALAT)
      IF(COSA.EQ.0.0)COSA=0.00001
      S1=SQRT(1.+EPS0)
      G0=1.-(EPS0/4.)+((7.*EPS0*EPS0)/64.)-((15.*(EPS0**3))/256.)
      G2=(EPS0/8.)-(.0625*EPS0*EPS0)+((145.*(EPS0**3))/2048.)
      G4=((5.*EPS0*EPS0)/256.)-((5.*(EPS0**3))/256.)
      G6=(29.*(EPS0**3))/6144.
      SIGM=(DIST*G0)/B1
      U1P=ATAN2(TAN1,(COSA*S1))
      SIN2P=SIN(2.*U1P)
      SIN4P=SIN(4.*U1P)
      TSS=((EPS0/4.)-((EPS0*EPS0)/8.))
      S12=2.*U1P-TSS*SIN2P-((EPS0*EPS0)/128.)*SIN4P
      SIGP=S12+SIGM
      T1=SIGM+(2.*G2*SIN(SIGM))*COS(SIGP)
      T2=(2.*G4*SIN(2.*SIGM))*COS(2.*SIGP)
      T3=(2.*G6*SIN(3.*SIGM))*COS(3.*SIGP)
      U2P=U1P+T1+T2+T3
      SINU1=TAN1/(SQRT(1.+EPS+(TAN1*TAN1)))
      C=SQRT(C2) 
      SINU2=(((B1*C)/B0)*SIN(U2P))-((EPS-EPS0)/(1.+EPS0))*SINU1
      U2=ASIN(SINU2)
      SINP1=SINU2/(SQRT(1.-E2*(COS(U2)*COS(U2))))
      BLAT=ASIN(SINP1)
      A1=B1*(SQRT(1.+EPS0))
      IF(COS(U2).EQ.0.0)U2=U2-0.00001
      Q1=(A1*COS(U2P))/(A0*COS(U2))
      if(q1.gt.-1..and.q1.lt.1.) Q2=ACOS(Q1)
      if(q1.eq.1) q2=0.
      if(q1.eq.-1.) q2=pi
      if(q1.gt.1.) go to 100
      if(q1.lt.-1.) go to 200
300   X1=SIN1*SINA
      AMU=ATAN2(X1,COSA)
      AZ=AZ/RAD
      U2P=U2P/RAD
      IF(AZ.GT.180.0)Q2=-Q2
      IF(U2P.GT.180.0.OR.U2P.LT.0.0)Q2=-Q2
      DLAMB=Q2-AMU
      BLON=DLAMB+ALON
      BLAT=BLAT/RAD
      BLON=BLON/RAD
      IF(ABS(BLON).GT.180.0)BLON=BLON-SIGN(360.,BLON)
      ALAT=ALT
      ALON=ALN
      RETURN
c100   write(6,*)'Flag in "Q2=ACOS(Q1)",Q1 = ',q1,';Q2 taken as 0'
100   q2=0.
      go to 300
c200   write(6,*)'     Flag in "Q2=ACOS(Q1)", Q1= ',q1,';Q2 taken as PI'
200   q2=pi
      go to 300
      END
c-------------------------------------------------------------
      SUBROUTINE lecdata(c,filein,jp,zsmac)
      parameter(NPROFMAX=50,NBTRAMAX=6000,NSTAMAX=200)
      dimension c(NPROFMAX,90,180),cmoy(NPROFMAX)
      character*80 ligne,filein,toto
      dimension zsmac(1)
c---------------------------------------------------------
c Cette routine lit 3SMAC et sort la vitesse en fonction 
c de la latitude et longitude.
c---------------------------------------------------------

c     Les fichiers 3SMAC sont donnes tout les 2*2 degres
c     avec val(ilat,ilon)
c     colat de 0 a 180 degres en partant du pole Nord
c     lon de 0 a 360 degres  
  
c    ouverture des fichiers
      toto='/util01/sismo/3-SMAC/para/'//filein(:lnblnk(filein))
      open(7,file=toto)
       read(7,'(a)')ligne
       write(*,'(a)')ligne
       write(*,*)'(Valeurs extrapolee a la prof de ',zsmac(jp-11),' km)'
       read(7,'(a)')ligne
       read(7,'(a)')ligne

      cmoy(jp)=0
       do 50 j=1,180 
       read(7,*)(c(jp,i,j),i=1,90)
      cmoy(jp)=cmoy(jp)+c(jp,i,j)
 50   continue
       cmoy(jp)=cmoy(jp)/(90*180)
      close(7)
      end
c-------------------------------------------------------------
      SUBROUTINE writedata(c,a1,a2,filein,jp,corr,zone,zsmac,cref)
      parameter(NPROFMAX=50,NBTRAMAX=6000,NSTAMAX=200)
      real c(NPROFMAX,90,180),a1(NPROFMAX,90,180)
      real a2(NPROFMAX,90,180)
      real zone(4)
      real cref
      dimension zsmac(1)
      integer sauti,sautj
      dimension cp(90,180),alpha(90,180),phiao(90,180)
      character*80 filein,fileout,fileaout,fileaout2

      common/c1/ nbtra,nsta,ncouch
c---------------------------------------------------------
c Ce pgm calcule la perturbation de Vs de 3SMAC par rapport
c a la vitesse de reference choisie.
c Il ecrit les resultats dans des fichiers :
c      -lon lat parametres
c---------------------------------------------------------

c     Les fichiers 3SMAC sont donnes tout les 2*2 degres
c     avec val(ilat,ilon)
c     ilat=colat de 0 a 180 degres en partant du pole Nord
c     ilon=lon de 0 a 360 degres

      pi=4.*atan(1.)
      degtorad=pi/180.0
      radtodeg=180.0/pi      
      iout=14
      iaout=15
      iaout2=16
    
      ilat=90
      ilon=180
      jj=jp-11

      write(*,*)zone(1),zone(2),zone(3),zone(4)
      ulat=zone(1)
      dlat=zone(2)
      wlon=zone(3)
      elon=zone(4)

c     calcul du pourcentage de perturbation par rapport a la valeur
c     moyenne de reference choisie.


        print*,'La valeur de reference pour le plot'
        print*,'du modele initial est :',cref

       do 110 i=1,ilat
       do 110 j=1,ilon
        cp(i,j)=(c(jp,i,j)-cref)/cref*100
 110     continue

       fileout= filein(:lnblnk(filein))//'.xyz'
       fileaout= 'an'//filein(:lnblnk(filein))//'.xyz'
       fileaout2= 'an2'//filein(:lnblnk(filein))//'.xyz'
       open(iout,file=fileout)
       open(iaout,file=fileaout)
       open(iaout2,file=fileaout2)

c     Conversion des indices colonnes en latitude et longitude
       sauti=0    
       sautj=0    

       do 130 i=1,ilat
        sauti=sauti+1
        if (sauti.eq.2)sauti=0
        lat= 90-(i+(i-1))
       do 130 j=1,ilon
        sautj=sautj+1
        if (sautj.eq.2)sautj=0
        lon=j+(j-1)
c  MODIF pour maggy le 14/06/2002----
        if(lon.gt.180)lon=lon-360
c---------------------------------
        if(lat.ge.dlat.and.lat.le.ulat.and.
     *lon.ge.wlon.and.lon.le.elon)then
       alpha(i,j)=sqrt((a1(jp,i,j)**2)+(a2(jp,i,j)**2))
       alpha(i,j)= 100*(alpha(i,j)/c(jp,i,j))
c   Facteur d'echelle tel que sur les cartes corr%aniso=1cm
c   1 unite est representee par 2.54cm (1pouce) par gmt
       alpha(i,j)=alpha(i,j)/(2.54*corr)
       phiao(i,j)= 0.5*radtodeg*atan2(a2(jp,i,j),a1(jp,i,j))
c   Phiao est l'angle par rapport au nord. Dans GMT
c   psxy plot les vecteurs par rapport a l'Est.
c   On convertit
          phiao(i,j)=90.-phiao(i,j)

       write(iout,'(2i5,f8.3)')lat,lon,cp(i,j)
        if(sauti.eq.0.and.sautj.eq.0)then
         write(iaout,'(2i5,2f8.3)')lat,lon,phiao(i,j),alpha(i,j)
         write(iaout2,'(2i5,2f8.3)')lat,lon,(phiao(i,j)+180),alpha(i,j)
        endif
       endif
  
 130    continue

      close(iout)
      close(iaout)
      close(iaout2)
      end

c----------------------------------------------------------
      subroutine trapez(npoint,f,ds,sint)
      dimension f(*)
      sint=0.
      nmax=npoint-1
      do 1 ns=2,nmax
  1   sint=sint+2.*f(ns)
      sint=sint+f(1)+f(npoint)
      sint=sint*ds*0.5
      return
      end
c----------------------------------------------------------
      subroutine lecintomo(sta,stat,lats,lons,nomtra,late,lone,err)
      parameter(NPROFMAX=50,NBTRAMAX=6000,NSTAMAX=200)
      real late(NBTRAMAX),lone(NBTRAMAX)
      real lats(NSTAMAX),lons(NSTAMAX)
      real vbid,err(NBTRAMAX,NPROFMAX)
      character*4  stat(NSTAMAX)
      character*10 sta(NBTRAMAX)
      character*80 nomtra(NBTRAMAX)
      character*80 filenem
      
      common/c1/ nbtra,nsta,ncouch
      
      write(*,*)'On supose qu\'il ya 32 profondeurs dans intomodes'
      write(*,*)'en input'
      write(*,*)'Reinitialiser la variable ncouch si necessaire'
c     On lit le fichier intomdesVs
      write(*,*)''
      write(*,*)'Nom du fichier intomodes?'
      write(*,*)'/eric4/SAMERICA/VS.5850/intomodesVs'
      write(*,*)'en verifiant que le modele contient bien 32 couches '
      write(*,*)'et en supposant que la seconde est bien celle a 100 km'
      read(*,*)filenem
      open(10,status='old',file=filenem)
      read(10,*) nsta,nbtra
      do 10 i=1,nsta
      read(10,1001)stat(i),lats(i),lons(i)
c     write(*,*)stat(i),lats(i),lons(i)
10    continue
      do 20 i=1,nbtra
      read(10,1002) sta(i),nomtra(i)
      read(10,*)late(i),lone(i)
      read(10,*) (vbid,jp=1,ncouch)
      read(10,*) (err(i,jp),jp=1,ncouch)
20    continue
      close(10)
1000  format(a32,f7.3,3x,f7.3,1x,a58,i1)
1001  format(a4,2x,2f9.4)
1002  format(a10,3x,a)
c1002  format(a10,3x,a32)
      end
c----------------------------------------------------------
      subroutine findcoord(stat,sta,slat,slon,lats,lons)
      parameter(NPROFMAX=50,NBTRAMAX=6000,NSTAMAX=200)
      real lats(NSTAMAX),lons(NSTAMAX)
      character*4 stat(NSTAMAX),sta,car
      
      common/c1/nbtra,nsta,ncouch

      
c     On retrouve les coord de la station
      icompt=0
      do 30  j=1, nsta
      car=stat(j)
      if(sta(:lnblnk(sta)).eq.car(:lnblnk(car))) then
      slat=lats(j)
      slon=lons(j)
       else
        icompt=icompt+1
       endif
30    continue
      if(icompt.eq.nsta) then
             write(*,*)car
             stop 'Une station est non trouvee'
      endif
      end
c----------------------------------------------------------
      SUBROUTINE  LISSE(N,X,Y,DY,S,A,B,C,D,R)
C-----------------------------------------------------------------------
C     CETTE SUBROUTINE PERMET LE CALCUL DE LA FONCTION SPLINE D'AJUSTEME
C        D'ORDRE 2 SUR LES POINTS X(I),Y(I),I=1,N .
C        LES X(I) SONT DONNES DANS L'ORDRE CROISSANT  . F EST CHOISIE TELLE
C        QUE SIGMA(((F(X(I))-Y(I))/DY(I))**2) )S .
C     VALEURS DE SORTIES : A,B,C,D (VECTEURS DE DIMENSION N+1)
C        F(T)=A(I)+B(I)*H+C(I)*H**2+D(I)*H**3  AVEC H=T-X(I-1) , X(I-1))
C        F(T)=A(1)+B(1)*(T-X(1)) SI T)X(1)
C        F(T)=A(N+1)+B(N+1)*(T-X(N)) SI T>X(N)
C     R : VECTEUR DE TRAVAIL EN DOUBLE PRECISION DE DIMENSION 8*(N+1)
C
C-----------------------------------------------------------------------
c     DIMENSION X(1),Y(1),DY(1),A(1),B(1),C(1),D(1),R(1)
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
