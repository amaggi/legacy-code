
c---------------------------------------------------------
      subroutine read_nomcoo(elat,elon,slat,slon,ntraj)
      parameter(NBTRAMAX=60000)
      real elat(NBTRAMAX),elon(NBTRAMAX)
      real slat(NBTRAMAX),slon(NBTRAMAX)
      character*80 filename, bid

      write(*,*) 'Filename of .txt event file: '
      read(*,*) filename
      write(*,*) 'Number of paths to read: '
      read(*,*) ntraj

      if (ntraj.gt.NBTRAMAX) goto 1001

      ifile=12
      open(unit=ifile, file=filename)
     
      do 10, i=1,ntraj
        read(ifile,*) bid, elat(i), elon(i), slat(i), slon(i)
10    continue

      close(ifile)


      return
1001  STOP 'Increase NBTRAMAX in subroutine read_nomcoo'
      end
c---------------------------------------------------------
      subroutine readtxt(nlines,line,lat,lon)
      parameter(NMAX=500)
      integer nlines
      character*80 line(NMAX), filename
      real lat(NMAX), lon(NMAX)

      write(*,*) 'Filename of .txt event file: '
      read(*,*) filename
      write(*,*) 'Number of events to read: '
      read(*,*) nlines
      if (nlines.gt.NMAX) goto 1001

      ifile=12

c     read all the lines
      open(unit=ifile, file=filename)
      do 10, i=1,nlines
        read(ifile,'(a)') line(i)
10    continue

      rewind(ifile)
c     read only lats and lons
      do 20, i=1,nlines
        read(ifile,'(28x,f8.3,1x,f9.3,24x)') lat(i),lon(i)
20    continue
      close(ifile)

      return
1001  STOP 'Increase NMAX in subroutine readtxt'
      end

c---------------------------------------------------------

      subroutine readsta(nlines,lat,lon,net,sta)
      parameter(NMAX=500)
      integer nlines
      character*80 filename, bidline
      character*2 net(NMAX)
      character*5 sta(NMAX)
      real lat(NMAX), lon(NMAX)

      write(*,*) 'Filename of station file file: '
      read(*,*) filename
      write(*,*) 'Number of lines to ignore: '
      read(*,*) nskip
      write(*,*) 'Number of stations to read: '
      read(*,*) nlines
      if (nlines.gt.NMAX) goto 1001

      ifile=12

c     read all the lines
      open(unit=ifile, file=filename)
      do 10, i=1,nskip
        read(ifile,'(a)') bidline
10    continue

c     read network, station, lats and lons
      do 20, i=1,nlines
        read(ifile,*) net(i), sta(i), lat(i), lon(i)
20    continue
      close(ifile)

      return
1001  STOP 'Increase NMAX in subroutine readsta'
      end


c---------------------------------------------------------
      SUBROUTINE writedata(c,a1,a2,filein,jp,zone,cref)
      parameter(NPROFMAX=50,NBTRAMAX=60000)
      real c(NPROFMAX,90,180),a1(NPROFMAX,90,180)
      real a2(NPROFMAX,90,180)
      real zone(4)
      real cref
      integer sauti,sautj
      dimension cp(90,180),alpha(90,180),phiao(90,180)
      character*80 filein,fileout,fileaout,fileaout2

      common/c1/ nbtra,ncouch
c---------------------------------------------------------
c Ce pgm calcule la perturbation de Vs de 3SMAC par rapport
c a la vitesse de reference choisie.
c Il ecrit les resultats dans des fichiers :
c      -lon lat parametres
c---------------------------------------------------------

c     Les fichiers PREM sont donnes tout les 2*2 degres
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

c      write(*,*)zone(1),zone(2),zone(3),zone(4)
      ulat=zone(1)
      dlat=zone(2)
      wlon=zone(3)
      elon=zone(4)

c     calcul du pourcentage de perturbation par rapport a la valeur
c     moyenne de reference choisie.

c       print*,'La valeur de reference pour le plot'
c       print*,'du modele initial est :',cref

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
        if(lat.ge.dlat.and.lat.le.ulat.and.
     *lon.ge.wlon.and.lon.le.elon)then
       alpha(i,j)=sqrt((a1(jp,i,j)**2)+(a2(jp,i,j)**2))
       alpha(i,j)= 100*(alpha(i,j)/c(jp,i,j))
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

c----------------------------------------------------------
c     READ GEOMETRY INFORMATION ONLY FROM INTOMODES FILE
c----------------------------------------------------------
      subroutine lecintomogeom(lats,lons,late,lone,nbtra)
      parameter(NBTRAMAX=60000)
      real late(NBTRAMAX),lone(NBTRAMAX)
      real lats(NBTRAMAX),lons(NBTRAMAX)
      real vbid,errbid
      character*80 nomtrabid
      character*80 filenem

c----------------------------------------------------------
c     Number of layers - important!
c----------------------------------------------------------
      ncouch=34
c----------------------------------------------------------
c     Get filename from user
c----------------------------------------------------------
      write(*,*)"On supose qu\'il y a ",ncouch, 
     *           " profondeurs dans intomodes"
      write(*,*)'en input'
      write(*,*)'Reinitialiser la variable ncouch si necessaire'
      write(*,*)'Expect number of paths as first line of intomodes file'
c     On lit le fichier intomdesVs
      write(*,*)''
      write(*,*)'Nom du fichier intomodes?'
      read(*,*)filenem
c----------------------------------------------------------
c     Read the intomodes file into memory
c----------------------------------------------------------
      open(10,status='old',file=filenem)
      read(10,*) nbtra
      write(*,*) nbtra, ' TRAJETS'
      if (nbtra.gt.NBTRAMAX) goto 1001
      do 20 i=1,nbtra
        read(10,1002) nomtrabid
        read(10,*)late(i),lone(i),lats(i),lons(i)
        read(10,*) (vbid,jp=1,ncouch)
        read(10,*) (errbid,jp=1,ncouch)
20    continue
      close(10)
1002  format(a80)
      return
1001  STOP 'Exceeded number of paths: NBTRAMAX'
 
      end
c----------------------------------------------------------
c     READ INTOMODES FILE V3 (no station list)
c----------------------------------------------------------
      subroutine lecintomo(lats,lons,nomtra,late,lone,err)
      parameter(NPROFMAX=50,NBTRAMAX=60000)
      real late(NBTRAMAX),lone(NBTRAMAX)
      real lats(NBTRAMAX),lons(NBTRAMAX)
      real vbid,err(NBTRAMAX,NPROFMAX)
      character*80 nomtra(NBTRAMAX)
      character*80 filenem
      integer i,NBTRAMAX
      
      common/c1/ nbtra,ncouch
      
c----------------------------------------------------------
c     Number of layers - important!
c----------------------------------------------------------
      ncouch=34
c----------------------------------------------------------
c     Get filename from user
c----------------------------------------------------------
      write(*,*)"On supose qu\'il y a ",ncouch, 
     *           " profondeurs dans intomodes"
      write(*,*)'en input'
      write(*,*)'Reinitialiser la variable ncouch si necessaire'
c     On lit le fichier intomdesVs
      write(*,*)''
      write(*,*)'Nom du fichier intomodes?'
      read(*,*)filenem
c----------------------------------------------------------
c     Read the intomodes file into memory
c----------------------------------------------------------
      open(10,status='old',file=filenem)
      read(10,*) nbtra
      write(*,*) nbtra, ' TRAJETS'
      do 20 i=1,nbtra
        read(10,1002) nomtra(i)
        read(10,*)late(i),lone(i),lats(i),lons(i)
        read(10,*) (vbid,jp=1,ncouch)
        read(10,*) (err(i,jp),jp=1,ncouch)
20    continue
      close(10)
1002  format(a80)
      end


c-------------------------------------------------------------
      SUBROUTINE lec3smac(c,filein,jp,zsmac)
      parameter(NPROFMAX=50)
      dimension c(NPROFMAX,90,180),cmoy(NPROFMAX)
      character*80 ligne,filein,toto
      real zsmac
      dimension zsmac(*)
c---------------------------------------------------------
c Cette routine lit 3SMAC et sort la vitesse en fonction 
c de la latitude et longitude.
c---------------------------------------------------------

c     Les fichiers 3SMAC sont donnes tout les 2*2 degres
c     avec val(ilat,ilon)
c     colat de 0 a 180 degres en partant du pole Nord
c     lon de 0 a 360 degres  
  
c    ouverture des fichiers
      toto='/home/alessia/pacific/3-SMAC/para/'//filein(:lnblnk(filein))
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
