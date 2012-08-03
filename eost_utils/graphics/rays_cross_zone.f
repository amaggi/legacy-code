c
c     calculates the path density from the intomodes file
c
      parameter(NBTRAMAX=60000)
      real slat,slon,elat,elon
      real lat, lon, latkm, lonkm, reflon

      real late(NBTRAMAX),lone(NBTRAMAX) 
      real lats(NBTRAMAX),lons(NBTRAMAX) 
 
      real minlat, maxlat, minlon, maxlon, clat, clon
      real zone, kmdeg

      logical inzone(NBTRAMAX)

      integer nbtra


      pi=4.*atan(1.)
      degtorad=pi/180.0
      radtodeg=180.0/pi      

c     read intomodes file
      write(*,*)'Lecture du fichier intomodesVs'
      call lecintomogeom(lats,lons,late,lone,nbtra)

c     ask for path step in degrees
      write(*,*) 'Path step in degrees : '
      read (*,*) pathstep

c     zone is correlation length
      zone=400
      kmdeg=111.195

c     ask for center (for azimuth calculation)
      write(*,*) 'clat, clon : '
      read(*,*) clat, clon 
   

      isor=12
      open(isor,file='rays_in_zone.dat')

c     start loop over paths

      do 30 itra=1,nbtra
        inzone(itra)=.false.
        elat=late(itra)
        elon=lone(itra)
        slat=lats(itra)
        slon=lons(itra)

        call mygrt(clat,clon,slat,slon,disd,az12,az21,gc)
c       only bother if the station itself is outside of the zone
        if(disd*kmdeg.gt.2*zone) then

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
            if(disd*kmdeg.le.zone) inzone(itra)=.true.
60        continue 
        endif
30    continue

c START AGAIN HERE - OUTPUT TO FILE

      write(isor,*) '#Rays crossing zone : '
      write(isor,*) '#Centered : ',clat,' - ',clon
      write(isor,*) '#Radius: ',zone

      do 90, itra=1,nbtra
        if(inzone(itra)) then 
c         calculate the event azimuth as seen from the center
          elat=late(itra)
          elon=lone(itra)
          call mygrt(elat,elon,clat,clon,disd,az12,az21,gc)
          azes=az21
          write(isor,*) lats(itra), lons(itra), 
     *                                 late(itra), lone(itra), azes
        endif
90    continue

      close(isor)

      return

      end


