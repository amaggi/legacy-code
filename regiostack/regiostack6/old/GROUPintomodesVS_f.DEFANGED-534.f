      parameter(NPROFMAX=32,NBTRAMAX=150000,NSTAMAX=1500)

      real vin(NBTRAMAX,NPROFMAX),errin(NBTRAMAX,NPROFMAX)
      real sigmoy(NBTRAMAX,NPROFMAX)
      real colate(NBTRAMAX),late(NBTRAMAX),lone(NBTRAMAX)
      real latsta(NBTRAMAX),lonsta(NBTRAMAX)
      integer tracer(NBTRAMAX),iclass(NBTRAMAX)

      real vsom(NPROFMAX),errsom(NPROFMAX),errmax(NPROFMAX)
      real vmoy(NPROFMAX),sigma(NPROFMAX)
      real vkksom(NPROFMAX)

      real lats(NSTAMAX),lons(NSTAMAX)
      logical mask(NBTRAMAX)
      real latsom,lonsom
      integer compti,comptimin,comptimax,kk
      character*4 toto 
      character*4  stat(NSTAMAX)
      character*4 sta2(NBTRAMAX)
      character*12 sta
      character*80 nomtra(NBTRAMAX),nameout
      character*80 filenem
     
c     16/09/2003
c     Ce programme demarre d'une liste de trajets epicentres-stations :
c     Il prend le 1er trajet de la liste A puis dresse la liste des autres
c     trajets dont l'epicentre se trouve dans un rayon de 100 km du premier.
c     De cette liste il extrait la liste B des trajets enregistres ? la meme 
c     station.
c     Puis il fait la meme chose en continuant a partir du premier trajet suivant
c     qui n'appartient pas a la liste B.
c     
c      
      ncouch=32
      isor=12
      pi=4.*atan(1.)
      degtorad=pi/180.0
      radtodeg=180.0/pi
 
      write(*,*)'Lecture du fichier intomodesVs'
      write(*,*)'On suppose qu il ya 32 profondeurs dans intomodes'
      write(*,*)'en input'
      write(*,*)'Reinitialiser la variable ncouch si necessaire'
      write(*,*)''
      write(*,*)'Nom du fichier intomodes?'
      read(*,*)filenem
      open(10,status='old',file=filenem)
      open(12,file='intomodesVs.group')
      open(13,file='gcin.xy')
      write(13,'(1a)')'>'
      open(14,file='gcgroup.xy')
      write(14,'(1a)')'>'

      read(10,*) nsta,nbtra
      write(12,*) nsta,nbtra
      do 10 i=1,nsta
      read(10,1001)stat(i),lats(i),lons(i)
      write(12,1001)stat(i),lats(i),lons(i)
10    continue

      do 20 i=1,nbtra
      mask(i)=.TRUE.
      iclass(i)=0
      do 211 jp=1,ncouch
        sigmoy(i,jp)=0
211   continue
      read(10,1002) sta,nomtra(i)
      sta2(i)=sta(:4)
c     if(index(nomtra(i),sta2(i)).eq.0) then
c      write(*,*)sta2(i),nomtra(i)
c      stop'Erreur dans intomodesVs'
c     endif
      read(10,*)late(i),lone(i)
      write(13,*)late(i),lone(i)
      colate(i)=(90-late(i))*degtorad
      late(i)=late(i)*degtorad
      lone(i)=lone(i)*degtorad

      read(10,*) (vin(i,jp),jp=1,ncouch)
      read(10,*) (errin(i,jp),jp=1,ncouch)

      icomptsta=0
      do 21 ij=1,nsta
      if(index(stat(ij),sta2(i)).ne.0) then
         icomptsta=icomptsta+1
         write(13,*)lats(ij),lons(ij)
         write(13,'(1a)')'>'
         latsta(i)=lats(ij)
         lonsta(i)=lons(ij)
      endif
21    continue
      if(icomptsta.gt.1)stop'icomptsta> 1'
20    continue

      do 22 jp=1,ncouch
       vkksom(jp)=0.
22    continue

      kk=0
      latsom=0
      lonsom=0
      comptimin=999999
      comptimax=1
      do 25 i=1,nbtra
        compti=1
        latsom=late(i)*radtodeg
        lonsom=lone(i)*radtodeg
        tracer(compti)=i
        do 251 jp=1,ncouch
           vsom(jp)=vin(i,jp)
           errsom(jp)=errin(i,jp)
           errmax(jp)=errin(i,jp)
           sigma(jp)=0.
 251    continue
        if(mask(i).eqv..FALSE.) go to 25
        kk=kk+1
        mask(i)=.FALSE.

        do 26 j=i+1,nbtra
          if(mask(j).eqv..TRUE.)then
          sdmimj=cosdel(colate(i),lone(i),colate(j),lone(j))
            if(acos(sdmimj)*radtodeg.lt.2.)then
c             if(index(sta2(i),sta2(j)).ne.0) then
              if(sta2(i).eq.sta2(j)) then
                mask(j)=.FALSE.
                compti=compti+1
                tracer(compti)=j
                latsom=latsom+late(j)*radtodeg
                if(abs(lone(j)*radtodeg-lone(i)*radtodeg).gt.180)then
                    if(lone(j).lt.0)lone(j)=lone(j)+360
                    if(lone(j).gt.0)lone(j)=lone(j)-360
                endif
                lonsom=lonsom+lone(j)*radtodeg
                do 261 jp=1,ncouch
                  vsom(jp)=vsom(jp)+vin(j,jp)
                  errsom(jp)=errsom(jp)+errin(j,jp)
                  if(errin(j,jp).gt.errmax(jp))errmax(jp)=errin(j,jp)
 261          continue
              endif
            endif
          endif
26    continue

       iclass(compti)=iclass(compti)+1
       if(compti.le.comptimin)comptimin=compti
       if(compti.ge.comptimax)comptimax=compti
       do 252 jp=1,ncouch
          vmoy(jp)=vsom(jp)/compti
          do 253 it=1,compti
             sigma(jp)=sigma(jp) + 
     *(vin(tracer(it),jp)-vmoy(jp))**2
 253      continue
          sigma(jp)=sqrt(sigma(jp)/compti)
          sigmoy(compti,jp)=sigmoy(compti,jp)+ sigma(jp)
          sigma(jp)=sigma(jp)/sqrt(float(compti))
          sigma(jp)=max(errmax(jp),sigma(jp))
          vkksom(jp)=vkksom(jp)+ vmoy(jp)
 252   continue
       if (kk.lt.10) then
          write(12,1003) sta2(i),".z",kk," nbre de trajets:",compti
c    *,(tracer(ik),ik=1,compti)
        elseif (kk.lt.100) then
          write(12,1004) sta2(i),".z",kk," nbre de trajets:",compti
c    *,(tracer(ik),ik=1,compti)
        elseif (kk.lt.1000) then
          write(12,1005) sta2(i),".z",kk," nbre de trajets:",compti
c    *,(tracer(ik),ik=1,compti)
        elseif (kk.lt.10000) then
          write(12,1006) sta2(i),".z",kk," nbre de trajets:",compti
c    *,(tracer(ik),ik=1,compti)
        elseif (kk.lt.100000) then
          write(12,1007) sta2(i),".z",kk," nbre de trajets:",compti
c    *,(tracer(ik),ik=1,compti)
        elseif (kk.lt.1000000) then
          write(12,1008) sta2(i),".z",kk," nbre de trajets:",compti
c    *,(tracer(ik),ik=1,compti)
       endif
       latsom=latsom/compti
       lonsom=lonsom/compti
       
       write(14,*)latsom,lonsom
       write(14,*)latsta(i),lonsta(i)
       write(14,'(1a)')'>'
       write(12,*)latsom,lonsom
       write(12,'(33(f8.4,2x))')(vmoy(jp),jp=1,ncouch)
       write(12,'(33(f8.4,2x))')(sigma(jp),jp=1,ncouch)

25    continue

      write(*,*) kk ,' trajets retenus'
      write(*,*) '     Groupe minimum ',comptimin
      write(*,*) '     Groupe maximum ',comptimax

      write(*,*)'----------------------------------'
      prof=50.
      do 27 icouch=1,ncouch
         vkksom(icouch)=vkksom(icouch)/kk
         if(prof.lt.700.) then
           prof=prof+25
         else
           prof=prof+50
         endif
         write(*,1100)prof,vkksom(icouch)

         in=14+icouch
         write(toto,'(i4.4)') int(prof)
         nameout= 'sigmoycouch'//toto//'.xy'
         open(in,file=nameout)

c        write(in,*)comptimax-comptimin+1
         do 270 ic = comptimin,comptimax
         if(iclass(ic).gt.1) then
           sigmoy(ic,icouch)=sigmoy(ic,icouch)/iclass(ic)
           write(in,*) ic,sigmoy(ic,icouch) 
c        else
c          sigmoy(ic,icouch)=0
         endif
270   continue
          close(in)  
27    continue

      inn=14+ncouch+1
      open(inn,file='nclusterpergroup.xy')
c     write(inn,*)comptimax-comptimin+1
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
c1003  format(a,a,i1,5x,a,30i5)
1003  format(a,a,i1,5x,a,i5)
1004  format(a,a,i2,5x,a,i5)
1005  format(a,a,i3,5x,a,i5)
1006  format(a,a,i4,5x,a,i5)
1007  format(a,a,i5,5x,a,i5)
1008  format(a,a,i6,5x,a,i5)
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
