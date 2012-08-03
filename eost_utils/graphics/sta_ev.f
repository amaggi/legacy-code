c
c     calculates the path density from the intomodes file
c
      parameter(NBTRAMAX=60000)
      real slat,slon,elat,elon
      real lat, lon, latkm, lonkm, reflon

      real late(NBTRAMAX),lone(NBTRAMAX) 
      real lats(NBTRAMAX),lons(NBTRAMAX) 

      integer nbtra



c     read intomodes file
      write(*,*)'Lecture du fichier intomodesVs'
      call lecintomogeom(lats,lons,late,lone,nbtra)

      isor=12
      open(isor,file='pseudo-nomcoo')

c     start loop over paths

      write(isor,*) '#Stla      Stlo       Evla       Evlo'
      do 30 i=1,nbtra
        write(isor,*) lats(i), lons(i), late(i), lone(i)
30    continue

      close(isor)

      return

      end


