c---------------------------------------------------------
      SUBROUTINE writedata(c,a1,a2,filein,jp,zone,cref)
      parameter(NPROFMAX=50,NBTRAMAX=40000,NSTAMAX=300)
      real c(NPROFMAX,90,180),a1(NPROFMAX,90,180)
      real a2(NPROFMAX,90,180)
      real zone(4)
      real cref
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
      subroutine lecintomo(sta,stat,lats,lons,nomtra,late,lone,err)
      parameter(NPROFMAX=50,NBTRAMAX=40000,NSTAMAX=300)
      real late(NBTRAMAX),lone(NBTRAMAX)
      real lats(NSTAMAX),lons(NSTAMAX)
      real vbid,err(NBTRAMAX,NPROFMAX)
      character*4  stat(NSTAMAX)
      character*10 sta(NBTRAMAX)
      character*80 nomtra(NBTRAMAX)
      character*80 filenem
      integer i
      
      common/c1/ nbtra,nsta,ncouch
      
      write(*,*)"On supose qu\'il ya 34 profondeurs dans intomodes"
      write(*,*)'en input'
      write(*,*)'Reinitialiser la variable ncouch si necessaire'
c     On lit le fichier intomdesVs
      write(*,*)''
      write(*,*)'Nom du fichier intomodes?'
      read(*,*)filenem
      open(10,status='old',file=filenem)
      read(10,*) nsta,nbtra
      do 10 i=1,nsta
      read(10,1001)stat(i),lats(i),lons(i)
10    continue
      do 20 i=1,nbtra
      read(10,1002) sta(i),nomtra(i)
c     write(*,*) sta(i), nomtra(i)

      read(10,*)late(i),lone(i)
      ncouch=34
      read(10,*) (vbid,jp=1,ncouch)
      read(10,*) (err(i,jp),jp=1,ncouch)
20    continue
      close(10)
1000  format(a32,f7.3,3x,f7.3,1x,a58,i1)
1001  format(a4,2x,2f9.4)
1002  format(a10,1x,a)
c1002  format(a10,3x,a32)
      end
c----------------------------------------------------------
      subroutine findcoord(stat,sta,slat,slon,lats,lons)
      parameter(NPROFMAX=50,NBTRAMAX=40000,NSTAMAX=300)
      real lats(NSTAMAX),lons(NSTAMAX)
      character*4 stat(NSTAMAX),sta,car
      
      common/c1/nbtra,nsta,ncouch

      
c     On retrouve les coord de la station
      icompt=0
      do 30  j=1, nsta
      car=stat(j)

c      write(*,*) "stat(j) is ", stat(j), " sta is ",sta
c      write(*,*) stat(j)
c      write(*,*) "sta(:lnblnk(sta)) is ",sta(:lnblnk(sta))

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


c--------------------------------------------------------------------------
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
