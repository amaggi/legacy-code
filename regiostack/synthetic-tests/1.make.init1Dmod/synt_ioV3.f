c---------------------------------------------------------
c---------------------------------------------------------
c     I/O routines for synthetic 3D model generation
c---------------------------------------------------------
c---------------------------------------------------------
c $Id: synt_ioV3.f,v 1.3 2004/05/04 08:08:27 alessia Exp $
c---------------------------------------------------------

c---------------------------------------------------------
c     getilin(lat) returns line index for latitude in degrees
c     given a gridstep
c---------------------------------------------------------
      function getilin(lat,gridstep)
      integer ind
      real gridstep,lat
      ind=int((90-lat)/gridstep)+1
      getilin=ind
      end function
      
c---------------------------------------------------------
c     geticol(lon) returns column index for longitude in degrees
c---------------------------------------------------------
      function geticol(lon,gridstep)
      integer ind
      real gridstep,lon
      ind=int(lon/gridstep)+1
      geticol=ind
      end function
      
 
c---------------------------------------------------------
c     getlat(ilin) returns latitude in degrees for line index
c     given a gridstep
c---------------------------------------------------------
      function getlat(ilin,gridstep)
      integer ilin
      real gridstep,lat
      lat= 90-((ilin-1)*gridstep+gridstep/2)
      getlat=lat
      end function
      
c---------------------------------------------------------
c     getlon(icol) returns longitude in degrees for column index
c---------------------------------------------------------
      function getlon(icol,gridstep)
      integer icol
      real gridstep,lon
      lon= (icol-1)*gridstep+gridstep/2
      getlon=lon
      end function
      
c---------------------------------------------------------
c     writedata(c,a1,a2,filein,jp,zone,cref,gstep,nanis)
c     writes the data out in gmt readable format
c---------------------------------------------------------
      SUBROUTINE writedata(c,a1,a2,filein,jp,zone,cref,gstep,nanis)
      include 'param_V3.h'
      real c(NPROFMAX,NLIN,NCOL),a1(NPROFMAX,NLIN,NCOL)
      real a2(NPROFMAX,NLIN,NCOL)
      real zone(4)
      real cref
      integer sauti,sautj
      dimension cp(NLIN,NCOL),alpha(NLIN,NCOL),phiao(NLIN,NCOL)
      character*80 filein,fileout,fileaout,fileaout2

      common/c1/ nbtra,ncouch,nlat,nlon
c---------------------------------------------------------
c Ce pgm calcule la perturbation de Vs de 3SMAC par rapport
c a la vitesse de reference choisie.
c Il ecrit les resultats dans des fichiers :
c      -lon lat parametres
c---------------------------------------------------------

c     Les fichiers PREM sont donnes tout les gstep*gstep degres
c     avec val(ilat,ilon)
c     ilat=colat de 0 a 180 degres en partant du pole Nord
c     ilon=lon de 0 a 360 degres

      pi=4.*atan(1.)
      degtorad=pi/180.0
      radtodeg=180.0/pi      
      iout=14
      iaout=15
      iaout2=16
    
      nlat=int(180/gstep)
      nlon=int(360/gstep)

c     calcul du pourcentage de perturbation par rapport a la valeur
c     moyenne de reference choisie.

      do 110 i=1,nlat
      do 110 j=1,nlon
        cp(i,j)=(c(jp,i,j)-cref)/cref*100
110   continue

      fileout= filein(:lnblnk(filein))//'.xyz'
      fileaout= 'an'//filein(:lnblnk(filein))//'.xyz'
      fileaout2= 'an2'//filein(:lnblnk(filein))//'.xyz'
      open(iout,file=fileout)
      open(iaout,file=fileaout)
      open(iaout2,file=fileaout2)

c     Conversion des indices colonnes en latitude et longitude
      sauti=0    
      sautj=0    

      do 130 i=1,nlat
        sauti=sauti+1
        if (sauti.eq.nanis)sauti=0
c       lat= 90-((i-1)*gstep)
        lat=getlat(i,gstep)
      do 130 j=1,nlon
        sautj=sautj+1
        if (sautj.eq.nanis)sautj=0
c       lon=(j-1)*gstep
        lon=getlon(j,gstep)
        alpha(i,j)=sqrt((a1(jp,i,j)**2)+(a2(jp,i,j)**2))
        alpha(i,j)= 100*(alpha(i,j)/c(jp,i,j))
        phiao(i,j)= 0.5*radtodeg*atan2(a2(jp,i,j),a1(jp,i,j))
c       Phiao est l'angle par rapport au nord. Dans GMT
c       psxy plot les vecteurs par rapport a l'Est.
c       On convertit
        phiao(i,j)=90.-phiao(i,j)

        write(iout,'(2i5,f8.3)')lat,lon,cp(i,j)
        if(sauti.eq.0.and.sautj.eq.0)then
          write(iaout,'(2i5,2f8.3)')lat,lon,phiao(i,j),alpha(i,j)
          write(iaout2,'(2i5,2f8.3)')lat,lon,(phiao(i,j)+180),alpha(i,j)
        endif
  
130   continue

      close(iout)
      close(iaout)
      close(iaout2)
      end



c----------------------------------------------------------

c----------------------------------------------------------
c     READ INTOMODES FILE V3 (no station list)
c----------------------------------------------------------
      subroutine lecintomo(lats,lons,nomtra,late,lone,err)
      include 'param_V3.h'
      real late(NBTRAMAX),lone(NBTRAMAX)
      real lats(NBTRAMAX),lons(NBTRAMAX)
      real vbid,err(NBTRAMAX,NPROFMAX)
      character*80 nomtra(NBTRAMAX)
      character*80 filenem
      integer i
      
      common/c1/ nbtra,ncouch,nlat,nlon
      
c----------------------------------------------------------
c     Number of layers - important!
c----------------------------------------------------------
      if (ncouch.eq.0) then
        ncouch=34
      endif
c----------------------------------------------------------
c     Get filename from user
c----------------------------------------------------------
c     write(*,*) 'Combien de couches dans intomodes? '
c     read(*,*) ncouch
       write(*,*)"On supose qu\'il y a ",ncouch, 
     *           " profondeurs dans intomodes"
       write(*,*)'en input'
c     write(*,*)'Reinitialiser la variable ncouch si necessaire'
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
        if (90.0-lats(i) < 0.01 ) lats(i)=lats(i)-0.01
        if (lats(i)+90.0 < 0.01 ) lats(i)=lats(i)+0.01
        read(10,*) (vbid,jp=1,ncouch)
        read(10,*) (err(i,jp),jp=1,ncouch)
20    continue
      close(10)
1002  format(a80)
      end

      subroutine lecintomogeom(lats,lons,late,lone,nbtra)
      include 'param_V3.h'
      real late(NBTRAMAX),lone(NBTRAMAX)
      real lats(NBTRAMAX),lons(NBTRAMAX)
      real vbid,err
      character*80 nomtra
      character*80 filenem
      integer i, nbtra
      
c----------------------------------------------------------
c     Get filename from user
c----------------------------------------------------------
      write(*,*) 'Combien de couches dans intomodes? '
      read(*,*) ncouch
      write(*,*)"On supose qu\'il y a ",ncouch, 
     *           " profondeurs dans intomodes"
      write(*,*)'en input'
c     write(*,*)'Reinitialiser la variable ncouch si necessaire'
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
        read(10,1002) nomtra
        read(10,*)late(i),lone(i),lats(i),lons(i)
        read(10,*) (vbid,jp=1,ncouch)
        read(10,*) (err,jp=1,ncouch)
20    continue
      close(10)
1002  format(a80)
      end

c-------------------------------------------------------------
c     reads PREM parameters and fills in starting c array
c     note: the action of reading PREM fixes the global 
c           variable ncouch, and fills in prof array
c-------------------------------------------------------------
      subroutine lecPREM(c,a1,a2,gstep,ncouch,prof)
      include 'param_V3.h'
      real c,a1,a2,prof,gstep, vs
      integer ncouch
      dimension c(NPROFMAX,NLIN,NCOL)
      dimension a1(NPROFMAX,NLIN,NCOL)
      dimension a2(NPROFMAX,NLIN,NCOL)
      dimension prof(NPROFMAX), vs(NPROFMAX)

      nlat=int(180/gstep)
      nlon=int(360/gstep)

      il=7
      write(*,*) 'Reading PREM'
      open(il,file='des.in',status='old')
      read(il,'(a)') title
      read(il,*) nbid
      read(il,*) ncouch, nvar
      
      do 10, jp=1,ncouch
        read(il,*) prof(jp),vs(jp),bidlogq
        write(*,*) prof(jp), vs(jp)
        do 20,i=1,nlat 
        do 20,j=1,nlon 
          c(jp,i,j)=vs(jp)
          a1(jp,i,j)=0.0
          a2(jp,i,j)=0.0
20      continue
10    continue
      close(il)

      end

c-------------------------------------------------------------
c     reads 3smac parameters and fills in starting c array
c-------------------------------------------------------------
      subroutine lec3smac(c,a1,a2,gstep,ncouch,prof,i50)
      include 'param_V3.h'
      dimension c(NPROFMAX,NLIN,NCOL)
      dimension a1(NPROFMAX,NLIN,NCOL)
      dimension a2(NPROFMAX,NLIN,NCOL)
      dimension c_smac(NPROFMAX,90,180)
      dimension c_smac_interp(NPROFMAX,90,180)
      character*80 ligne,filein,toto
      real gstep
      real zsmac(23),prof(NPROFMAX)
      real lat, lon, latsmac, lonsmac
      integer nlat, nlon,i50
      logical found

c     work arrays for smoothing
      real ctemp(NPROFMAX), dctemp(NPROFMAX)
      real aa(NPROFMAX),bb(NPROFMAX),cc(NPROFMAX),dd(NPROFMAX)
      real*8 rr(8*(NPROFMAX+1))

      data zsmac/50,60,70,80,90,100,115,130,150,175,200,
     *250,300,350,390,430,460,510,530,560,610,640,680/



c     Les fichiers 3SMAC sont donnes tout les 2*2 degres
c     avec val(ilat,ilon)
c     colat de 0 a 180 degres en partant du pole Nord
c     lon de 0 a 360 degres  

      do 10, jp=12,34
c       ouverture des fichiers
        kp=jp-11
        write(filein,'("VS.",i2)')jp
        toto='/home/alessia/share/3-SMAC/para/'//filein(:lnblnk(filein))
        open(7,file=toto)
        read(7,'(a)')ligne
c       write(*,'(a)')ligne
c       write(*,*)'(Valeurs a la prof de ',zsmac(kp),' km)'
        read(7,'(a)')ligne
        read(7,'(a)')ligne
        do 50 j=1,180 
          read(7,*)(c_smac(kp,i,j),i=1,90)
50      continue
        close(7)
10    continue

c     Set up number of layers and final depth array
      ncouch=25
      prof(1)=40
      prof(2)=50
      do 16 ii=1,ncouch
        if(ii.gt.2)prof(ii)=prof(ii-1)+25
        if(prof(ii).gt.700)prof(ii)=prof(ii-1)+50
16    continue

      sig=0.01
c     depth-interpolate the smac c array to wfm depths
      do 20, i=1,90
      do 20, j=1,180
        do 25, kp=1,23
          ctemp(kp)=c_smac(kp,i,j)
          dctemp(kp)=ctemp(kp)/100.0
25      continue
c       create the arrays to do the interpolation
        call lisse(23,zsmac,ctemp,dctemp,sig,aa,bb,cc,dd,rr)
c       interpolate for each depth
        do 27, jp=1,ncouch
          depth=prof(jp)
          c_smac_interp(jp,i,j)=smoo(23,zsmac,aa,bb,cc,dd,depth)
27      continue
20    continue



c     Now transfer the c_smac_interp array to the c_array, taking into
c     account the fact that the gstep is different.
c     NOTE: am going to assume that 3smac values refer to center of
c           2x2 degree cells.  
c           Will assign 3smac value to a lat lon point if point is
c           closer in lat, lon to center of 3smac cell than half the
c           3smac gridsize

c     find number of lat/lons in c grid
      nlat=int(180/gstep)
      nlon=int(360/gstep)

c     iterate over c grid
      do 60, i=1,nlat
      do 60, j=1,nlon
        found=.false.
        lat=getlat(i,gstep)
        lon=getlon(j,gstep)

c       iterate over smac grid exiting early when you can
c       first the latitude grid - do longitude only if latitude is ok
        do 70, ic=1,90
          latsmac=getlat(ic,2.)
          if (abs(latsmac-lat).lt.1.) then
          do 80, jc=1,180
             lonsmac=getlon(jc,2.)
             if (abs(lonsmac-lon).lt.1.) then
               do 90 jp=1,ncouch
                 if(i50.eq.1) then
c                  want each depth to have 50km 3SMAC velocity
                   c(jp,i,j)=c_smac_interp(2,ic,jc)
                 else
                   c(jp,i,j)=c_smac_interp(jp,ic,jc)
                 endif
c                write(*,*) i,j,jp,c(jp,i,j), c_smac_interp(jp,ic,jc)
                 a1(jp,i,j)=0
                 a2(jp,i,j)=0
90             continue
               found=.true.
               goto 65
             endif
80        continue
          endif
70      continue

c       get here either because we have broken out of the loop or 
c       because we failed to find 3smac value for this i,j pair
65      if (.not.found) stop 'Failed to allocate 3smac value to c array'

60    continue

      end
c-------------------------------------------------------------

c-------------------------------------------------------------
c     subroutine readfin(c,da1,da2,ndep,ili,icol)
c-------------------------------------------------------------
      subroutine readfin(profs,c,da1,da2,ndep,ili,icol)
      include 'param_V3.h'
      real profs,c,da1,da2
      integer ndep,ili,icol
      dimension profs(NPROFMAX),c(NPROFMAX,NLIN,NCOL),
     *          da1(NPROFMAX,NLIN,NCOL),da2(NPROFMAX,NLIN,NCOL)
      
      ivs=11
      ianvs=12
      if (ili >  NLIN) stop 'NLIN not big enough in readfin'
      if (icol > NCOL) stop 'NCOL not big enough in readfin'

      ivs=11
      ianvs=12
      open(unit=ivs,file='vsfin_input')
      open(unit=ianvs,file='anvsfin_input')
      do idep=1,ndep
        read(ivs,*)tt
        read(ianvs,*)tt
        profs(idep)=tt
        do i=1,ili
          read(ivs,'(10(f10.6,1x))')(c(idep,i,j),j=icol/2+1,icol)
          read(ivs,'(10(f10.6,1x))')(c(idep,i,j),j=1,icol/2)
        enddo
        do i=1,ili
          read(ianvs,'(10(f10.6,1x))')(da1(idep,i,j),j=icol/2+1,icol)
          read(ianvs,'(10(f10.6,1x))')(da1(idep,i,j),j=1,icol/2)
        enddo
        do  i=1,ili
          read(ianvs,'(10(f10.6,1x))')(da2(idep,i,j),j=icol/2+1,icol)
          read(ianvs,'(10(f10.6,1x))')(da2(idep,i,j),j=1,icol/2)
        enddo
      enddo
      close(ivs)
      close(ianvs)
 

      end

c-------------------------------------------------------------
c     subroutine writefin(c,da1,da2,ndep,ili,icol)
c-------------------------------------------------------------
      subroutine writefin(profs,c,da1,da2,ndep,ili,icol)
      include 'param_V3.h'
      real profs,c,da1,da2
      integer ndep,ili,icol
      dimension profs(NPROFMAX),c(NPROFMAX,NLIN,NCOL),
     *          da1(NPROFMAX,NLIN,NCOL),da2(NPROFMAX,NLIN,NCOL)
      
      ivs=11
      ianvs=12
      open(unit=ivs,file='vsfin')
      open(unit=ianvs,file='anvsfin')
      do 10,idep=1,ndep
        write(ivs,*)profs(idep) 
        write(ianvs,*)profs(idep) 
        do 30 i=1,ili
          write(ivs,'(10(f10.6,1x))')(c(idep,i,j),j=icol/2+1,icol)
30        write(ivs,'(10(f10.6,1x))')(c(idep,i,j),j=1,icol/2)
        do 50 i=1,ili
          write(ianvs,'(10(f10.6,1x))')(da1(idep,i,j),j=icol/2+1,icol)
50        write(ianvs,'(10(f10.6,1x))')(da1(idep,i,j),j=1,icol/2)
        do 60 i=1,ili
          write(ianvs,'(10(f10.6,1x))')(da2(idep,i,j),j=icol/2+1,icol)
60        write(ianvs,'(10(f10.6,1x))')(da2(idep,i,j),j=1,icol/2)
10    continue
      close(ivs)
      close(ianvs)
      end
 

c-------------------------------------------------------------
c     subroutine write_coverage
c-------------------------------------------------------------
      subroutine write_coverage(lats,lons,late,lone,nbtra)
      include 'param_V3.h'
      real late(NBTRAMAX),lone(NBTRAMAX)
      real lats(NBTRAMAX),lons(NBTRAMAX)
      character*80 filenem
      integer i
 
      write(*,*)'Name of output file for coverage (gmt format)?'
      read(*,*)filenem
      open(10,status='unknown',file=filenem)
      do i=1,nbtra
        write(10,101) late(i), lone(i)
        write(10,101) lats(i), lons(i)
        write(10,*) '>'
      enddo
      close(10)

101   format(2(1x,f10.4))
      end subroutine

