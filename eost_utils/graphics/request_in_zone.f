c
c     calculates the path density from the intomodes file
c
      parameter(NMAX=500,NBTRAMAX=60000)
      real slat,slon,elat,elon
      real lat, lon, latkm, lonkm, reflon

      real late(NMAX),lone(NMAX) 
      real lats(NMAX),lons(NMAX) 

      real nomlate(NBTRAMAX),nomlone(NBTRAMAX) 
      real nomlats(NBTRAMAX),nomlons(NBTRAMAX) 
 
      character*2 net(NMAX)
      character*5 sta(NMAX)
      real zone, kmdeg 

      integer nev,nsta,nnomcoo
      character*80 evline(NMAX), filename

      logical inzone, havepath


      pi=4.*atan(1.)
      degtorad=pi/180.0
      radtodeg=180.0/pi      
      kmdeg=111.195

c     read txt file
      call readtxt(nev,evline,late,lone)
      call readsta(nsta,lats,lons,net,sta)
      call read_nomcoo(nomlate,nomlone,nomlats,nomlons,nnomcoo)

      write(*,*) 'Finished setup'


c     ask for path step in degrees
c     write(*,*) 'Path step in degrees : '
c     read (*,*) pathstep
      pathstep=1
      distance=1
      zone=400

      write(*,*) 'Center point lat, lon:'
      read(*,*) clat, clon

      write(*,*) 'Output file:'
      read(*,*) filename

      isor=12
      open(isor,file=filename)

      write(isor,*) '#Rays within zone : '
      write(isor,*) '#Center : ',clat,' - ',clon
      write(isor,*) '#Radius: ',zone

c     start loop over stations

      do 30 is=1,nsta
        slat=lats(is)
        slon=lons(is)
c       only bother if the station itself is outside of the box
        call mygrt(clat,clon,slat,slon,disd,az12,az21,gc)
        if(disd*kmdeg.le.2*zone) goto 30
        write(*,'(a,1x,a)') net(is), sta(is)
        do 35 ie=1,nev
          inzone=.false.
          elat=late(ie)
          elon=lone(ie)

          call mygrt(elat,elon,slat,slon,disd,az12,az21,gc)
          azes=az12
          do while(azes.lt.0.) 
            azes=azes+360.
          enddo
          do while(azes.gt.360.) 
            azes=azes-360.
          enddo
          nsteps= (disd/pathstep+1)

c  Boucle sur chaque point de l'azimuth epicentre-station

          do 60 jcle=1,nsteps 
            dstcle=pathstep*jcle*111.195 
            call mygds(elat,elon,azes,dstcle,blat,blon)
            call mygrt(blat,blon,clat,clon,disd,az12,az21,gc)
            if(disd*kmdeg.lt.zone) inzone=.true.
60        continue 
          if (inzone) then
c           check we don't have this path yet
            havepath=.false.
            do 80, i=1,nnomcoo
              elatnom=nomlate(i)
              elonnom=nomlone(i)
              slatnom=nomlats(i)
              slonnom=nomlons(i)
              call mygrt(elat,elon,elatnom,elonnom,disde,az12,az21,gc)
              call mygrt(slat,slon,slatnom,slonnom,disds,az12,az21,gc)
              if(disde.lt.distance.and.disds.lt.distance) then
                havepath=.true.
                goto 81
              endif
80          continue
81          if(.not.havepath) then 
              write(isor,'(a,1x,a,1x,a,1x,f8.3,1x,f9.3)')
     *                   evline(ie),net(is),sta(is),lats(is),lons(is)
            endif
          endif
35      continue
30    continue

      close(isor)

      return

      end


