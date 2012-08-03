      parameter(NPROFMAX=34,NBTRAMAX=150000)

      real vin(NPROFMAX),errin(NPROFMAX),errstat(NPROFMAX)
      integer iseuil(NPROFMAX)
      real late(NBTRAMAX),lone(NBTRAMAX)
      real lats(NBTRAMAX),lons(NBTRAMAX)

      character*30 sta
      character*80 filenem

c     20/11/2003
c     Version du code qui genere un fichier intomodesVS.V3
c     avec un format simplifie(coordonnee de la station
c     donnee pour chaque evenement)
c     16/09/2003
c     Ce programme reconstruit un fichier intomodesVs ou les erreurs
c     sont determinees de facon statistiques.
c     
c      
      ncouch=34
      isor=12
      pi=4.*atan(1.)
      degtorad=pi/180.0
      radtodeg=180.0/pi
 
c------------------------------------------------------------------
c     Open and read Errstat
c------------------------------------------------------------------
      open(13,file='Errstat')
      read(13,*)ilayer
      do 5 j=1,ilayer
        read(13,*) prof,errstat(j),iseuil(j)
5     continue
      close(13)

c------------------------------------------------------------------
c     Open and read grouped intomodes file
c------------------------------------------------------------------
      write(*,*)'Lecture du fichier intomodesVs'
      write(*,*)'On suppose qu il ya 34 profondeurs dans intomodes'
      write(*,*)'en input'
      write(*,*)'Reinitialiser la variable ncouch si necessaire'
      write(*,*)''
      write(*,*)'Nom du fichier intomodes?'
      read(*,*)filenem
      open(10,status='old',file=filenem)
      open(12,file='intomodesVs.group2')

      read(10,*) nbtra
      write(12,*)nbtra

c------------------------------------------------------------------
c     START LOOP OVER PATHS
c------------------------------------------------------------------
      do 20 i=1,nbtra
c------------------------------------------------------------------
c       Copy the cluster identification line to output
c       and store the number of paths in this cluster
c------------------------------------------------------------------
        if(i.lt.10)then
          read(10,1003) sta,inbre
          write(*,1003) sta,inbre
          write(12,1003) sta,inbre
        elseif (i.lt.100) then
          read(10,1004) sta,inbre
          write(12,1004) sta,inbre
        elseif (i.lt.1000) then
          read(10,1005) sta,inbre
          write(12,1005) sta,inbre
        elseif (i.lt.10000) then
          read(10,1006) sta,inbre
          write(12,1006) sta,inbre
        elseif (i.lt.100000) then
          read(10,1007) sta,inbre
          write(12,1007) sta,inbre
        elseif (i.lt.1000000) then
          read(10,1008) sta,inbre
          write(12,1008) sta,inbre
        endif
c------------------------------------------------------------------
c       Read lat lon information, Vs and Errors
c------------------------------------------------------------------
        read(10,*)late(i),lone(i),lats(i),lons(i)
        read(10,*) (vin(jp),jp=1,ncouch)
        read(10,*) (errin(jp),jp=1,ncouch)

c------------------------------------------------------------------
c       If there are fewer paths than required for good stats in 
c       this cluster then apply error from Errstat to this cluster
c       NB. Error is applied as sigma_m!
c------------------------------------------------------------------
        do 212 jp=1,ncouch
          if(inbre.le.iseuil(jp)) then
            errin(jp)=errstat(jp)/sqrt(float(inbre))
          endif
212     continue

        write(12,'(4(2x,f9.4))')late(i),lone(i),lats(i),lons(i)
        write(12,'(34(f8.4,2x))')(vin(jp),jp=1,ncouch)
        write(12,'(34(f8.4,2x))')(errin(jp),jp=1,ncouch)
20    continue
c------------------------------------------------------------------
c     END LOOP OVER PATHS
c------------------------------------------------------------------

      write(*,*)'----------------------------------'

      close(10)
      close(12)
      close(13)
      close(14)
1000  format(a32,f7.3,3x,f7.3,1x,a58,i1)
1001  format(a4,2x,2f9.4)
1002  format(a,a,a,i6)
1003  format(a25,i5)
1004  format(a26,i5)
1005  format(a27,i5)
1006  format(a28,i5)
1007  format(a29,i5)
1008  format(a30,i5)
1100  format(f7.2,1x,f7.3,' .05 0. 0. 0.005')

      end
