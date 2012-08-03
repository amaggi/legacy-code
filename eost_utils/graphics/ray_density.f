c
c     calculates the path density from the intomodes file
c
      parameter(NBTRAMAX=60000)
      real slat,slon,elat,elon
      real lat, lon, latkm, lonkm, reflon

      real late(NBTRAMAX),lone(NBTRAMAX) 
      real lats(NBTRAMAX),lons(NBTRAMAX) 

      integer ncount(360,720)

      real gridstep

      integer nbtra


      pi=4.*atan(1.)
      degtorad=pi/180.0
      radtodeg=180.0/pi      

c     read intomodes file
      write(*,*)'Lecture du fichier intomodesVs'
      call lecintomogeom(lats,lons,late,lone,nbtra)

c     ask for grid step in degrees
      write(*,*) 'Grid step in degrees : '
      read (*,*) gridstep

      pathstep=gridstep/2.0

      nlat=int(180/gridstep)
      nlon=int(360/gridstep)

      isor=12
      open(isor,status='new',file='raydensity')

c     start loop over paths

      nmin=0
      nmax=0
      do 30 itra=1,nbtra
        elat=late(itra)
        elon=lone(itra)
        slat=lats(itra)
        slon=lons(itra)

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

        ilatsave=-999
        ilonsave=-999
        do 60 jcle=1,nsteps 
          dstcle=pathstep*jcle*111.195 
          call mygds(elat,elon,azes,dstcle,blat,blon)
c         the array is in colatitude
          ilat=int((90-blat)/gridstep)
          if(blon.lt.0.)blon=360.+blon
          ilon=int(blon/gridstep)
   
          if(ilat.gt.nlat) goto 1001
          if(ilon.gt.nlon) goto 1002

          if (ilat.ne.ilatsave.or.ilon.ne.ilonsave) then
            ilatsave=ilat
            ilonsave=ilon
            ncount(ilat,ilon) = ncount(ilat,ilon) + 1
          endif

60      continue 
30    continue

c START AGAIN HERE - OUTPUT TO FILE

      do 90, ilat=1,nlat
        do 91, ilon=1,nlon
          lat=90-ilat*gridstep
          lon=ilon*gridstep
          write(isor,*) lat,lon,ncount(ilat,ilon)
91      continue
90    continue

      close(isor)

      return

1001  STOP 'Latitude index exceeds nlat'
1002  STOP 'Longitude index exceeds nlon'

      end


