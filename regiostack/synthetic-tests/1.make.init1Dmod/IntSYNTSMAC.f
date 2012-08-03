      parameter(NPROFMAX=50,NBTRAMAX=40000,NSTAMAX=300)
      parameter(NPERTMAX=20,NPTSMAX=10000,NCOUPMAX=20)
      real slat,slon,elat,elon,depth,sig,vscoup,profmin,profmax

      real v(NPROFMAX,NPTSMAX),c(NPROFMAX,90,180)
      real a1(NPROFMAX,90,180),a2(NPROFMAX,90,180)
      real cinlatc(NPROFMAX,NCOUPMAX,180),cingc(NPROFMAX,NPTSMAX)
      real vout(NPROFMAX),prof(NPROFMAX),vsout(NPROFMAX)
      real coupin(NPROFMAX),dcoupin(NPROFMAX)
      real lentsom(NPROFMAX),vsin(NPROFMAX),dvsin(NPROFMAX)
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
c    25 JUIN 2001 : CODE ADAPTE A PARTIR DE IntsurtrajetsVs1anis.f
c    en vue du calcul d'un test synthetique complet
c    incluant waveform modelling et inversion 3D.
c
c    Programme de calcul de modeles synthetiques dans 3SMAC
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
c    Cette version sort en plus une coupe a lat=17 deg du modele 3SMAC +dVS
c----------------------------------------------------------

      ncouch=34
      nvar=2
      isor=12
      pi=4.*atan(1.)
      degtorad=pi/180.0
      radtodeg=180.0/pi      
      zone(1)=89
      zone(2)=-89
      zone(3)=121
      zone(4)=301


      write(*,*)'Lecture du fichier intomodesVs'
      call lecintomo(sta,stat,lats,lons,nomtra,late,lone,err)

c  Lecture de 3SMAC pour toute les profondeurs comprises
c  entre 50 km et 660km (couches 12 a 34).
c     
      write(*,*)'Entrer 1 pour ajouter une perturbation de Vs a 3SMAC'
      read(*,*)ipert
      if(ipert.eq.1)then
        write(*,*)'Lat Lon(0->359) des 2 perturbations?'
        write(*,*)'(e.g. 16 40 16 44) '
        read(*,*)latpert1,lonpert1,latpert2,lonpert2
        write(*,*)'Prof min et max de la perturbation? '
        write(*,*)'(e.g. 300 460) '
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
       write(*,*)'(e.g. 40,0,40,30)'
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
       call lec3smac(c,filein,jp,zsmac)

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
c---Ecriture du modele en INPUT 3SMAC en format xyz--------------------
c   Ecritures des cartes a differentes profondeurs
      crefpass=cref(jp)
      call writedata(c,a1,a2,filein,jp,zone,zsmac,crefpass)
      if(comptpert.ne.0)ipass=1
      crefi(jp-11)=cref(jp)*100.
      dcrefi(jp-11)=crefi(jp-11)/10000.

10    continue

         sig=23.
         call lisse(23,zsmac,crefi,dcrefi,sig,aa,bb,cc,dd,rr)
c        prof(1)=75
c	 bad bad hack to cheat intomo3SMAC to have values for 50
c        km rather than 75km as the integrated velocities.
         prof(1)=50
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
c        prof(1)=75
c	 bad bad hack to cheat intomo3SMAC to have values for 50
c        km rather than 75km as the integrated velocities.
         prof(1)=50
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
c        prof(1)=75
c	 bad bad hack to cheat intomo3SMAC to have values for 50
c        km rather than 75km as the integrated velocities.
         prof(1)=50

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
      write(*,*)'Initialisation du modele de depart achevee'
      write(*,*)''

c--------FIN DES COUPES-----------------------------

       open(isor,status='new',file='intomo3SMAC')
       write(isor,*)nsta
       do 120 j=1,nsta
       write(isor,1001)stat(j),lats(j),lons(j)
120    continue
      
      do 125 jp=1,NPROFMAX
125    lentsom(jp)=0

      do 130 itra=1,nbtra
      nomtrajet=nomtra(itra)
      istart=index(nomtrajet,'19')
      is2=index(nomtrajet,'20')
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
180    continue

c -Interpolation des modeles 1D tout les 25 km 
      do 210 ii=12,34
      jj=ii-11
      vsin(jj)=vout(ii)*100.
      dvsin(jj)=vsin(jj)/10000.
210   continue

      sig=23. 
      call lisse(23,zsmac,vsin,dvsin,sig,aa,bb,cc,dd,rr)
      
      toto='des.S'//nomtrajet
      write(*,'(a)')toto 
      open(9,file=toto,status='new') 
      write(9,'(a)') 'Synthetic within 3SMAC + dvs'
      write(9,*)' 1'
      write(9,*)' 25 ',nvar

      do 220 ii=1,25
        depth=prof(ii)
        vsout(ii)=smoo(23,zsmac,aa,bb,cc,dd,depth)
        write(9,*)depth,vsout(ii)/100., 0.0050
        lentsom(ii)=lentsom(ii)+(100./vsout(ii))
220   continue
      close(9)

c       Modify here to only output integrated 3SMAC at 50 km depth
c       but keep the a priori errors the same
c       write(isor,1003) (vsout(ii)/100.,ii=1,25)
c       write(isor,1003) (err(itra,ii),ii=1,25)
        write(isor,1003) (vsout(1)/100.,ii=1,25)
        write(isor,1003) (err(itra,ii),ii=1,25)

130    continue

c  Affichage de la moyenne des lenteurs a 100 km utiliser comme modele
c  a priori
     

      do 230 ii=1,25
        write(*,*) prof(ii),1./(lentsom(ii)/nbtra),' .05 0. 0. 0.005'
230   continue

 
      write(*,*) 'ATTENTION : toutes les couches Vs de 3SMAC sont'
      write(*,*) 'prises en compte dans intomo3SMAC. Intsurtrajetanis.f'
      write(*,*) 'permet de prendre en compte une seule couche '
      
      close (isor)
1001  format(a4,2x,2f9.4)
1003  format(50f8.4)
      end


