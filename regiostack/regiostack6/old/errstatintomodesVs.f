      parameter(NPROFMAX=40,NBTRAMAX=150000,NSTAMAX=1500)

      real vin(NPROFMAX),errin(NPROFMAX),errstat(NPROFMAX)
      integer iseuil(NPROFMAX)
      real late(NBTRAMAX),lone(NBTRAMAX)

      real lats(NSTAMAX),lons(NSTAMAX)
      character*4  stat(NSTAMAX)
      character*4 sta2(NBTRAMAX)
      character*33 sta
      character*80 nomtra(NBTRAMAX)
      character*80 filenem,ligne
     
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
 
      open(13,file='Errstat')
      read(13,*)ilayer
      do 5 j=1,ilayer
      read(13,*) prof,errstat(j),iseuil(j)
5     continue
      close(13)

      write(*,*)'Lecture du fichier intomodesVs'
      write(*,*)'On suppose qu il ya 34 profondeurs dans intomodes'
      write(*,*)'en input'
      write(*,*)'Reinitialiser la variable ncouch si necessaire'
      write(*,*)''
      write(*,*)'Nom du fichier intomodes?'
      read(*,*)filenem
      open(10,status='old',file=filenem)
      open(12,file='intomodesVs.group2')

      read(10,*) nsta,nbtra
      write(12,*) nsta,nbtra
      do 10 i=1,nsta
      read(10,1001)stat(i),lats(i),lons(i)
      write(12,1001)stat(i),lats(i),lons(i)
10    continue

      do 20 i=1,nbtra
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
c     write(*,*)inbre
      sta2(i)=sta(:4)
      read(10,*)late(i),lone(i)
      read(10,*) (vin(jp),jp=1,ncouch)
      read(10,*) (errin(jp),jp=1,ncouch)

      do 212 jp=1,ncouch
      if(inbre.le.iseuil(jp)) then
         errin(jp)=errstat(jp)/sqrt(float(inbre))
      endif
212   continue

       write(12,*)late(i),lone(i)
       write(12,'(34(f8.4,2x))')(vin(jp),jp=1,ncouch)
       write(12,'(34(f8.4,2x))')(errin(jp),jp=1,ncouch)
20    continue

      write(*,*)'----------------------------------'

      close(10)
      close(12)
      close(13)
      close(14)
1000  format(a32,f7.3,3x,f7.3,1x,a58,i1)
1001  format(a4,2x,2f9.4)
1002  format(a,a,a,i6)
1003  format(a29,i5)
1004  format(a30,i5)
1005  format(a31,i5)
1006  format(a32,i5)
1007  format(a33,i5)
1008  format(a34,i5)
1100  format(f7.2,1x,f7.3,' .05 0. 0. 0.005')

      end
