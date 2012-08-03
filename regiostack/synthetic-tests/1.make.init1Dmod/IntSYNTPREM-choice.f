c sylvana: added perturb. magnitude in input parameters      
c alessia: modified to do checkerboard test

      parameter(NPROFMAX=50,NBTRAMAX=40000,NSTAMAX=300)
      parameter(NPTSMAX=8000,NCOUPMAX=20, NPERTMAX=20)
      real slat,slon,elat,elon,crefpass
      real lat, lon, latkm, lonkm, reflon

      real v(NPROFMAX,NPTSMAX),c(NPROFMAX,90,180)
      real a1(NPROFMAX,90,180),a2(NPROFMAX,90,180)
      real vout(NPROFMAX),prof(NPROFMAX)
      real vs(NPROFMAX),bidlogq
      real cref(NPROFMAX)
      real lentsom(NPROFMAX)
      integer n_cyl
      real cyl_lat(NPERTMAX), cyl_lon(NPERTMAX)
      real cyl_mag(NPERTMAX), cyl_rad(NPERTMAX)
      real cyl_minprof(NPERTMAX), cyl_maxprof(NPERTMAX)
     
      real mag, pert_size, profmin,profmax
      integer anis_type, psign, synt_type

      common latc(NCOUPMAX),lat1(NCOUPMAX),lon1(NCOUPMAX)
      common lat2(NCOUPMAX),lon2(NCOUPMAX) 
      real late(NBTRAMAX),lone(NBTRAMAX) 
      real lats(NSTAMAX),lons(NSTAMAX) 
      real err(NBTRAMAX,NPROFMAX) 
      real zone(4),len(NPTSMAX)
      character*2 nocoup
      integer isor,jc,nbcoup,nbcouplat
      integer icout(NCOUPMAX),icoutlat(NCOUPMAX)
      character*92 toto
      character*4 stat(NSTAMAX),sta2 
      character*10 sta(NBTRAMAX),stabid
      character*80 nomtra(NBTRAMAX),filein
      character*80 nomtrajet
      character*80 filecout,filecoutlat

      common/c1/nbtra,nsta,ncouch

c----------------------------------------------------------
c    25 JUIN 2001 : CODE ADAPTE A PARTIR DE IntsurtrajetsVs1anis.f
c    en vue du calcul d'un test synthetique complet
c    incluant waveform modelling et inversion 3D.
c
c    Programme de calcul de modeles synthetiques dans PREM + perturbation
c    Il calcul l'integrale des lenteurs sur un trajet donne.
c    les coordonnees des noms de trajets  et des epicentres
c    et stations sont dans un fichier format Montagner.
c
c    Cette version lit 
c        - Un model PREM 1D (des.in) extrapolle en chaque points geographique
c          ou aucune perturbation n'est ajoute 
c        - un fichier intomodesVs 
c          pour la couverture en trajet que l'on veut mettre dans le test.
c    Cette version permet 
c        - d'ajouter une anisotropie synthetique a 3SMAC.
c        - d'ajouter un dVs a PREM pour une region et un interval de profondeur donne 
c    
c----------------------------------------------------------
      il1=12
      isor=12
      pi=3.1415926
      degtorad=pi/180.0
      radtodeg=180.0/pi      
      zone(1)=89
      zone(2)=-89
      zone(3)=121
      zone(4)=301

c     center longitude of plot
      reflon=210

c     read intomodes file
      write(*,*)'Lecture du fichier intomodesVs'
      call lecintomo(sta,stat,lats,lons,nomtra,late,lone,err)

c     ask for synthetic type
      write(*,*) 'Type of synthetic test desired:'
      write(*,*) '1: Checkerboard'
      write(*,*) '2: Cylinders'
      read(*,*)synt_type
   
      if (synt_type .eq. 2 ) then
c     -----------------------------------------------------------------
c     cylinders
c     -----------------------------------------------------------------
      print *,'Number of cylindrical perturbations ',
     *         '(max ', NPERTMAX,')'
       read(*,*)n_cyl
       if (n_cyl.gt.NPERTMAX) then
         write(*,*) 'Too many perturbations'
         stop
       endif
       do 2 i=1,n_cyl
         print *,'Lon, lat, radius(km), mindepth, maxdepth, ',
     *           'magnitude (percentage) for cylinder', i
         write(*,*)'(e.g. 25, -15, 200, 50, 410, 0.03)'
         read(*,*)cyl_lon(i),cyl_lat(i),cyl_rad(i),cyl_minprof(i),
     *            cyl_maxprof(i),cyl_mag(i)
2      continue
c     -----------------------------------------------------------------

      else
c     -----------------------------------------------------------------
c     checkerboard
c     -----------------------------------------------------------------
c       ask for checkerboard parameters
        write(*,*) 'Lateral size of chekerboard (km)'
        write(*,*)'(e.g. 1000) '
        read(*,*) pert_size
        write(*,*)'Prof min et max de la perturbation? '
        write(*,*)'(e.g. 300 450) '
        read(*,*) profmin,profmax
        write(*,*)'Magnitude (%) of perturbation? (e.g. 0.05)'
        read(*,*) mag
        write(*,*)'Anisotropic perturbation? (0=no, 1=yes)'
        read(*,*) anis_type
c     -----------------------------------------------------------------
      endif


c     -----------------------------------------------------------------
c     cross-sections    
c     -----------------------------------------------------------------
       print *,'Number of latitude cross-sections (max ', NCOUPMAX,')'
       print *,'and lat of the cross-sections '
       write(*,*)'(e.g. 2,15,17)'
       read(*,*)nbcouplat,(latc(ij),ij=1,nbcouplat)
       print *,'Number of great circle path cross-sections ',
     *         '(max ', NCOUPMAX,')'
       read(*,*)nbcoup
       do 4 i=1,nbcoup
         print *,'Coord of 2 points to define the minor arc along which'
         write(*,*)'the cross-section ',i,' is plotted '
         write(*,*)'(lon1,lat,lon2,lat2)'
         write(*,*)'(e.g. 25,-15,52,25)'
         read(*,*)lon1(i),lat1(i),lon2(i),lat2(i)
4      continue

c      write gmt file of crossections
       open(9,file='gc.cross.xy')
       write(9,'(1a)')'>'
       do 5 jc=1,nbcoup
         write(9,*)lat1(jc),lon1(jc)
         write(9,*)lat2(jc),lon2(jc)
         write(9,'(1a)')'>'
c        open a file for each crossection
         icout(jc)=24+jc
         write(nocoup,'(i2.2)')jc
         filecout= 'VSPREMcoupgc'//nocoup//'.xyz'
         open(icout(jc),file=filecout)
5      continue
       do 6 jc=1,nbcouplat
c        open a file for each latitude cut
         icoutlat(jc)=44+jc
         write(nocoup,'(i2.2)')jc
         filecoutlat= 'VSPREMcouplat'//nocoup//'.xyz'
         open(icoutlat(jc),file=filecoutlat)
6      continue
       close(9)
c     -----------------------------------------------------------------

c-----------------------------------------------
c     ADD PERTURBATIONS
c-----------------------------------------------
      write(*,*)'Lecture de PREM'
      open(il1, file='des.in', status='old')
      read(il1,'(a)') titre
      read(il1,*) nbid
      read(il1,*) ncouch,nvar

c------------------------------------------------------------------
c    start loop over depth layers
c------------------------------------------------------------------
      do 10 jp=1,ncouch
        if (jp.lt.10)then 
          write(filein,'("VS.0",i1)')jp
        else
          write(filein,'("VS.",i2)')jp
        endif
        read(il1,*) prof(jp),vs(jp),bidlogq

        cref(jp)=0
        do 12 i=1,90
          do 12 j=1,180
c           find the latitude and longitude of this point
            lat= 90-(i+(i-1))
            lon=j+(j-1)
c           set the velocity at each point to be PREM for that depth
            c(jp,i,j)=vs(jp)
c           set the anisotropy at each point to be 0 
            a1(jp,i,j)=0.0
            a2(jp,i,j)=0.0
c           apply the relevant perturbation
            if(synt_type .eq. 2) then
c-----------------------------------------------------------------
c             cylinders
c-----------------------------------------------------------------
               c(jp,i,j)=vs_cyl(c(jp,i,j),lat,lon,prof(jp),
     *                          cyl_lat, cyl_lon, cyl_rad,
     *                          cyl_maxprof, cyl_minprof,
     *                          cyl_mag,n_cyl)
            else
c-----------------------------------------------------------------
c             checkerboard
c-----------------------------------------------------------------
              c(jp,i,j)=vschecker(c(jp,i,j),lat,lon,prof(jp),pert_size,
     *                             profmin, profmax,mag)
              if (anis_type .eq. 1) then
                a2(jp,i,j)=a2checker(lat,lon,prof(jp),pert_size,
     *                             profmin, profmax)
              endif
            endif
c-----------------------------------------------------------------
c           add to reference velocity sum
            cref(jp)=cref(jp)+c(jp,i,j)
12      continue
c       reference velocity for depth jp is the mean of the velocities
        cref(jp)=cref(jp)/90/180
c-----------------------------------------------


c----------Ecriture du modele en INPUT en format GMT--------------------
c   Ecritures des cartes a differentes profondeurs
        crefpass=cref(jp)
        call writedata(c,a1,a2,filein,jp,zone,prof,crefpass)

c--------COUPES DANS LE MODELE INITIAL-----------------------------

        do 140 jc=1,nbcoup
          alat1=lat1(jc)
          alon1=lon1(jc)
          alat2=lat2(jc)
          alon2=lon2(jc)
          call mygrt(alat1,alon1,alat2,alon2,disd,az12,az21,gc)
          delta=disd
          azes=az12
          do while (azes.lt.0.) 
            azes=azes+360.
          enddo
          do while (azes.gt.360.) 
            azes=azes-360.
          enddo
          npas= int(5.*delta+1.)
          dpas= delta/float(npas)
          npat= npas+1

c       Boucle sur chaque point de l'azimuth epicentre-station
          do 140 jcle=1,npat
            dstcle=dpas*float(jcle-1)*111.195
            alat1=lat1(jc)
            alon1=lon1(jc)
            call mygds(alat1,alon1,azes,dstcle,blat,blon)
            i=int(((90-blat+1)/2)+0.5)
            if(blon.lt.0.)blon=360.+blon
            j=int(((blon+1)/2)+0.5)
            write(icout(jc),1200)-prof(jp),jcle,
     *           (c(jp,i,j)-cref(jp))/cref(jp)
140     continue
        do 17 ii=1,nbcoup
          ilatc=int(latc(ii))
          ilatc=int(((90-ilatc+1)/2)+0.5)
          do 17 jlon=1,180
            lon=jlon+(jlon-1)
            write(icoutlat(ii),1200)-prof(jp),lon,
     *           (c(jp,ilatc,jlon)-cref(jp))/cref(jp)
17      continue
1200    format(f8.2,2x,f8.2,f8.3)


10    continue
c------------------------------------------------------------------
c    end loop over depth layers
c------------------------------------------------------------------
      close(il1)

c    On ferme les fichiers contenant les coupes
      do 150 jc=1,nbcoup
        close(icout(jc))
        close(icoutlat(jc))
150   continue

      write(*,*)''
      write(*,*)'Initialisation du modele de depart achevee'
      write(*,*)''

c--------FIN DES COUPES-----------------------------

      open(isor,status='new',file='intomoPREM')
      write(isor,*)nsta
      do 20 j=1,nsta
        write(isor,1001)stat(j),lats(j),lons(j)
20    continue
      
      do 25 jp=1,ncouch
25    lentsom(jp)=0

      do 30 itra=1,nbtra
        nomtrajet=nomtra(itra)
        istart=index(nomtra(itra),'19')
        is2=index(nomtra(itra),'200')
        if(is2.lt.istart.and.is2.ge.1)istart=is2
        if(istart.eq.0)istart=is2
        nomtrajet=nomtra(itra)(istart:lnblnk(nomtra(itra)))
        elat=late(itra)
        elon=lone(itra)
c       sanity check
        if(elat.lt.-90.or.elat.gt.90.or.elon.lt.-180.or.elon.gt.180)then
          write(*,*) nomtraj, elat, elon
        endif
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
        do while(azes.lt.0.) 
          azes=azes+360.
        enddo
        do while(azes.gt.360.) 
          azes=azes-360.
        enddo
        npas= int(20.*delta+1.)
        dpas= delta/float(npas)
        npat= npas+1
c       write(*,*) dpas, npas, npat

c  Boucle sur chaque point de l'azimuth epicentre-station

        do 60 jcle=1,npat 
          dstcle=dpas*float(jcle-1)*111.195 
          call mygds(elat,elon,azes,dstcle,blat,blon)
          call mygrt(elat,elon,blat,blon,dstcle,azim,azback,gc)
          i=int(((90-blat+1)/2)+0.5)
          if(blon.lt.0.)blon=360.+blon
          j=int(((blon+1)/2)+0.5)

          do 70 jp=1,ncouch 
            v(jp,jcle)=c(jp,i,j)
            v(jp,jcle)=v(jp,jcle)+a1(jp,i,j)*cos(2*degtorad*azim)+
     *                 a2(jp,i,j)*sin(2*degtorad*azim)
            if(v(jp,jcle).eq.0) then
              write(*,*)'V NUL TRA',itra
              stop
            endif
70        continue 
 
60      continue 

c  Il ne reste plus qu'a integrer sur le trajet.
    
        do 80 jp=1,ncouch
          do 90 jcle=1,npat 
            len(jcle)=1./v(jp,jcle)
90        continue
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
80      continue

        toto='des.S'//nomtrajet(:lnblnk(nomtrajet))
c       write(*,'(a)')toto 
        open(9,file=toto,status='new') 
        write(9,'(a)') 'Synthetic within PREM + dvs'
        write(9,*)' 1'
        write(9,*)ncouch-2,nvar

        do 120 ii=3,ncouch
          write(9,*)prof(ii),vout(ii), 0.0050
          lentsom(ii)=lentsom(ii)+(1./vout(ii))
120     continue
        close(9)


        write(isor,1003) (vout(ii)/100.,ii=3,ncouch)
        write(isor,1003) (err(itra,ii),ii=3,ncouch)

30    continue

c  Affichage de la moyenne des lenteurs a 100 km utiliser comme modele
c  a priori
     

      do 130 ii=3,ncouch
        write(*,*) prof(ii),1./(lentsom(ii)/nbtra),' .05 0. 0. 0.005'
130   continue

      close (isor)
1001  format(a4,2x,2f9.4)
1003  format(50f8.4)
      end


