      parameter(NPROFMAX=34,NBTRAMAX=150000)

      real vin(NBTRAMAX,NPROFMAX),errin(NBTRAMAX,NPROFMAX)
      real sigmoy(NBTRAMAX,NPROFMAX)
      real colate(NBTRAMAX),late(NBTRAMAX),lone(NBTRAMAX)
      integer tracer(NBTRAMAX),iclass(NBTRAMAX)

      real vsom(NPROFMAX),errsom(NPROFMAX),errmax(NPROFMAX)
      real vmoy(NPROFMAX),sigma(NPROFMAX)
      real vkksom(NPROFMAX)

      real colats(NBTRAMAX),lats(NBTRAMAX),lons(NBTRAMAX)
      logical mask(NBTRAMAX)
      real colatsom,latsom,lonsom,colatstasom,latstasom,lonstasom
      integer compti,comptimin,comptimax,kk
      real sdmimj, sdmimjstat, distdeg
      character*4 toto 
      character*12 sta
      character*80 nomtra(NBTRAMAX),nameout
      character*80 filenem

c     20/11/2003
c     Version du code qui genere un fichier intomodesVS.V3
c     avec un format simplifie(coordonnee de la station
c     donnee pour chaque evenement)      
c     16/09/2003
c     Ce programme demarre d'une liste de trajets epicentres-stations :
c     Il prend le 1er trajet de la liste A puis dresse la liste des autres
c     trajets dont l'epicentre se trouve dans un rayon de 100 km du premier.
c     De cette liste il extrait la liste B des trajets enregistres a la meme 
c     station.
c     Puis il fait la meme chose en continuant a partir du premier trajet 
c     suivant qui n'appartient pas a la liste B.
c     
c      
      ncouch=34
      isor=12
      pi=4.*atan(1.)
      degtorad=pi/180.0
      radtodeg=180.0/pi
 
c------------------------------------------------------------------
c     Get the name of the intomodes file and the radius of the cluster
c------------------------------------------------------------------
      write(*,*)'Lecture du fichier intomodesVs'
      write(*,*)'On suppose qu il ya 34 profondeurs dans intomodes'
      write(*,*)'en input'
      write(*,*)'Reinitialiser la variable ncouch si necessaire'
      write(*,*)''
      write(*,*)'Nom du fichier intomodes?'
      read(*,*)filenem
      write(*,*)'Radius of the cluster in degrees?  (e.g. 2 )'
      read(*,*)distdeg

c------------------------------------------------------------------
c     Open the input file for reading, and the others for writing
c------------------------------------------------------------------
      open(10,status='old',file=filenem)
      open(12,file='intomodesVs.group')
      open(13,file='gcin.xy')
      write(13,'(1a)')'>'
      open(14,file='gcgroup.xy')
      write(14,'(1a)')'>'

c------------------------------------------------------------------
c     Read number of paths from input intomodesVs
c------------------------------------------------------------------
      read(10,*) nbtra
      write(12,*) nbtra

c------------------------------------------------------------------
c     START INITIALISATION LOOP OVER PATHS
c------------------------------------------------------------------
      write(*,*) 'Starting initialisation'
      do 20 i=1,nbtra
c------------------------------------------------------------------
c       Initalise mask to true, and sigma to 0 for cluster around
c       this path, then read path parameters.
c------------------------------------------------------------------
        mask(i)=.TRUE.
        iclass(i)=0
        do 211 jp=1,ncouch
          sigmoy(i,jp)=0
211     continue
        read(10,1002) sta,nomtra(i)
        read(10,*)late(i),lone(i),lats(i),lons(i)
c------------------------------------------------------------------
c       Fixup longitudes of both events and stations
c------------------------------------------------------------------
        if(lone(i).gt.180)then
          lone(i)=lone(i)-360
        else if(lone(i).le.-180)then
          lone(i)=lone(i)+360
        endif
        if(lons(i).gt.180)then
          lons(i)=lons(i)-360
        else if(lons(i).le.-180)then
          lons(i)=lons(i)+360
        endif
c------------------------------------------------------------------
c       Fudge station location for the south pole
c------------------------------------------------------------------
        if((lats(i).eq.-90).and.lons(i).eq.0.) then
          lats(i)=-89.9999
          lons(i)=0.0001
c------------------------------------------------------------------
c       Write to full path coverage GMT file
c------------------------------------------------------------------
        endif
        write(13,*)late(i),lone(i)
        write(13,*)lats(i),lons(i)
c------------------------------------------------------------------
c       Convert latitudes to colatitudes, and all coordinates
c       to radians.
c------------------------------------------------------------------
        colats(i)=(90-lats(i))*degtorad
        colate(i)=(90-late(i))*degtorad
        late(i)=late(i)*degtorad
        lats(i)=lats(i)*degtorad
        lone(i)=lone(i)*degtorad
        lons(i)=lons(i)*degtorad
c------------------------------------------------------------------
c       Read in Vs and error for path i
c------------------------------------------------------------------
        read(10,*) (vin(i,jp),jp=1,ncouch)
        read(10,*) (errin(i,jp),jp=1,ncouch)
20    continue
c------------------------------------------------------------------
c     END INITIALISATION LOOP OVER PATHS
c------------------------------------------------------------------

c------------------------------------------------------------------
c     kk counts the cluster index 
c------------------------------------------------------------------
      do 22 jp=1,ncouch
        vkksom(jp)=0.
22    continue

      kk=0
      comptimin=999999
      comptimax=1
c------------------------------------------------------------------
c     START EXTERNAL LOOP OVER PATHS
c------------------------------------------------------------------
      write(*,*) 'Starting to loop over paths...'
      do 25 i=1,nbtra
c------------------------------------------------------------------
c       if this path has already been assigned to a group, skip it
        if(mask(i).eqv..FALSE.) go to 25
c------------------------------------------------------------------
        compti=1
        latsom=late(i)*radtodeg
        colatsom=colate(i)*radtodeg
        lonsom=lone(i)*radtodeg
        latstasom=lats(i)*radtodeg
        colatstasom=colats(i)*radtodeg
        lonstasom=lons(i)*radtodeg

        tracer(compti)=i
        do 251 jp=1,ncouch
          vsom(jp)=vin(i,jp)
          errsom(jp)=errin(i,jp)
          errmax(jp)=errin(i,jp)
          sigma(jp)=0.
 251    continue

c------------------------------------------------------------------
c       Create a new cluster with index kk and assign this path to it
c------------------------------------------------------------------
        kk=kk+1
        mask(i)=.FALSE.

c------------------------------------------------------------------
c       START INTERNAL LOOP OVER REMAINING PATHS
c------------------------------------------------------------------
        do 26 j=i+1,nbtra
          if(mask(j).eqv..TRUE.)then
            sdmimj=cosdel(colate(i),lone(i),colate(j),lone(j))
            sdmimjstat=cosdel(colats(i),lons(i),colats(j),lons(j))
c------------------------------------------------------------------
c           if the j'th path has both ends within distdeg of the
c           i'th path, then assign to this cluster
c------------------------------------------------------------------
            if(acos(sdmimj)*radtodeg.lt.distdeg .and. 
     *         acos(sdmimjstat)*radtodeg.lt.distdeg )then
              mask(j)=.FALSE.
              compti=compti+1
              tracer(compti)=j
c------------------------------------------------------------------
c             sum latitudes
c------------------------------------------------------------------
              latsom=latsom+late(j)*radtodeg
              colatsom=colatsom+colate(j)*radtodeg
              latstasom=latstasom+lats(j)*radtodeg
              colatstasom=colatstasom+colats(j)*radtodeg
c------------------------------------------------------------------
c             sum longitudes - taking care over 180 degrees
c------------------------------------------------------------------
              if(abs(lone(j)*radtodeg-lone(i)*radtodeg).gt.180)then
                if(lone(j).lt.0)then
                  lone(j)=(lone(j)*radtodeg+360)*degtorad
                else if(lone(j).ge.0)then 
                  lone(j)=(lone(j)*radtodeg-360)*degtorad
                endif
              endif
              lonsom=lonsom+lone(j)*radtodeg
              if(abs(lons(j)*radtodeg-lons(i)*radtodeg).gt.180)then
                if(lons(j).lt.0)then
                  lons(j)=(lons(j)*radtodeg+360)*degtorad
                else if(lons(j).ge.0)then 
                  lons(j)=(lons(j)*radtodeg-360)*degtorad
                endif
              endif
              lonstasom=lonstasom+lons(j)*radtodeg

c------------------------------------------------------------------
c             sum vs and errors from path j
c------------------------------------------------------------------
              do 261 jp=1,ncouch
                vsom(jp)=vsom(jp)+vin(j,jp)
                errsom(jp)=errsom(jp)+errin(j,jp)
                if(errin(j,jp).gt.errmax(jp))errmax(jp)=errin(j,jp)
261           continue
            endif
          endif
26      continue
c------------------------------------------------------------------
c       END INTERNAL LOOP OVER REMAINING PATHS
c------------------------------------------------------------------

        iclass(compti)=iclass(compti)+1
        if(compti.le.comptimin)comptimin=compti
        if(compti.ge.comptimax)comptimax=compti
c------------------------------------------------------------------
c       Calculate mean and standard deviation of velocities
c------------------------------------------------------------------
        do 252 jp=1,ncouch
          vmoy(jp)=vsom(jp)/compti
          do 253 it=1,compti
            sigma(jp)=sigma(jp) + (vin(tracer(it),jp)-vmoy(jp))**2
253       continue
          sigma(jp)=sqrt(sigma(jp)/compti)
c------------------------------------------------------------------
c         sigmoy does NOT contain sigma_m, but the sum of all sigmas
c         for clusters with compti paths in them
c------------------------------------------------------------------
          sigmoy(compti,jp)=sigmoy(compti,jp)+ sigma(jp)
c------------------------------------------------------------------
c         turn sigma into a sigma_m by dividing by sqrt(compti)
c         use max of this value and the max error on a single path
c------------------------------------------------------------------
          sigma(jp)=sigma(jp)/sqrt(float(compti))
          sigma(jp)=max(errmax(jp),sigma(jp))
c------------------------------------------------------------------
c         vkksom holds the sum of the mean velocities at this depths
c         for all clusters
c------------------------------------------------------------------
          vkksom(jp)=vkksom(jp)+ vmoy(jp)
252     continue

c------------------------------------------------------------------
c       output an identifying line for this cluster
c------------------------------------------------------------------
        if (kk.lt.10) then
          write(12,1003) "z",kk," nbre de trajets:",compti
        elseif (kk.lt.100) then
          write(12,1004) "z",kk," nbre de trajets:",compti
        elseif (kk.lt.1000) then
          write(12,1005) "z",kk," nbre de trajets:",compti
        elseif (kk.lt.10000) then
          write(12,1006) "z",kk," nbre de trajets:",compti
        elseif (kk.lt.100000) then
          write(12,1007) "z",kk," nbre de trajets:",compti
        elseif (kk.lt.1000000) then
          write(12,1008) "z",kk," nbre de trajets:",compti
        endif

c------------------------------------------------------------------
c       calculate centroid lat and lon for both ends of the cluster
c------------------------------------------------------------------
        latsom=latsom/compti
        colatsom=colatsom/compti
        lonsom=lonsom/compti

        latstasom=latstasom/compti
        colatstasom=colatstasom/compti
        lonstasom=lonstasom/compti

c------------------------------------------------------------------
c       SANITY CHECK - if the centroids are too far from the base
c       (i) around which the cluster was made, there is a bug !
c------------------------------------------------------------------
        sdmimj=cosdel(colatsom*degtorad,lonsom*degtorad,
     *                colate(i),lone(i))
        if(acos(sdmimj)*radtodeg.gt.distdeg)then
          write(*,*)'Path ',i,' sdmimj ',acos(sdmimj)*radtodeg
          write(*,*)'centroid event at ',latsom,lonsom
          write(*,*)'i-th event at ',late(i)*radtodeg,lone(i)*radtodeg
        endif
        sdmimj=cosdel(colatstasom*degtorad,lonstasom*degtorad,
     *                colats(i),lons(i))
        if(acos(sdmimj)*radtodeg.gt.distdeg)then
          write(*,*)'Path ',i,' sdmimj ',acos(sdmimj)*radtodeg
          write(*,*)'centroid station at ',latstasom,lonstasom
          write(*,*)'i-th station at ',lats(i)*radtodeg,lons(i)*radtodeg
        endif

c------------------------------------------------------------------
c       Write to clustered coverage GMT file
c------------------------------------------------------------------
        write(14,*)latsom,lonsom
        write(14,*)latstasom,lonstasom
        write(14,'(1a)')'>'
c------------------------------------------------------------------
c       Write to output intomodes file
c------------------------------------------------------------------
        write(12,'(4(2x,f9.4))')latsom,lonsom,latstasom,lonstasom
        write(12,'(34(f8.4,2x))')(vmoy(jp),jp=1,ncouch)
        write(12,'(34(f8.4,2x))')(sigma(jp),jp=1,ncouch)

25    continue
c------------------------------------------------------------------
c     END EXTERNAL LOOP OVER PATHS
c------------------------------------------------------------------
      
   
c------------------------------------------------------------------
c     Write output statistics
c------------------------------------------------------------------
      write(*,*) kk ,' trajets retenus'
      write(*,*) '     Groupe minimum ',comptimin
      write(*,*) '     Groupe maximum ',comptimax

      write(*,*)'----------------------------------'

c------------------------------------------------------------------
c     START OUTPUT LOOP
c------------------------------------------------------------------
      prof=30.
      do 27 icouch=1,ncouch
        vkksom(icouch)=vkksom(icouch)/kk
        if(prof.lt.41.) then
          prof=prof+10
        else if(prof.lt.700.) then
          prof=prof+25
        else
          prof=prof+50
        endif
c------------------------------------------------------------------
c       write a priori model to screen
c------------------------------------------------------------------
        write(*,1100)prof,vkksom(icouch)

        in=14+icouch
        write(toto,'(i4.4)') int(prof)
        nameout= 'sigmoycouch'//toto//'.xy'
        open(in,file=nameout)

c------------------------------------------------------------------
c       write out mean sigma per number of paths in a cluster
c------------------------------------------------------------------
        do 270 ic = comptimin,comptimax
          if(iclass(ic).gt.1) then
            sigmoy(ic,icouch)=sigmoy(ic,icouch)/iclass(ic)
            write(in,*) ic,sigmoy(ic,icouch),iclass(ic)
          endif
270     continue
        close(in)  
27    continue
c------------------------------------------------------------------
c     END OUTPUT LOOP
c------------------------------------------------------------------

      inn=14+ncouch+1
      open(inn,file='nclusterpergroup.xy')
      do 28 ic = comptimin,comptimax
        if(iclass(ic).ne.0) write(inn,*) ic,iclass(ic)
28    continue
      close(inn)

      close(10)
      close(12)
      close(13)
      close(14)
1000  format(a32,f7.3,3x,f7.3,1x,a58,i1)
1001  format(a4,2x,2f9.4)
1002  format(a12,1x,a40)
1003  format(a,i1,5x,a,i5)
1004  format(a,i2,5x,a,i5)
1005  format(a,i3,5x,a,i5)
1006  format(a,i4,5x,a,i5)
1007  format(a,i5,5x,a,i5)
1008  format(a,i6,5x,a,i5)
1100  format(f7.2,1x,f7.3,' .05 0. 0. 0.005')

      end
c----------------------------------------------------------
      function cosdel(tetam1,phim1,tetam2,phim2)
      cosdel=0.
      c1=cos(tetam1)
      c2=cos(tetam2)
      s1=sin(tetam1)
      s2=sin(tetam2)
      cosdel=c1*c2+(s1*s2)*cos(phim1-phim2)
      return
      end
