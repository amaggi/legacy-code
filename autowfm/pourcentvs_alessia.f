      real*8 latcoup(10),iCout(10)
      real*8 profs(34),c(360,720),er(360,720),prof(34),vref(34)
      real*8 profi(34),vs(34),ro(34)
      real*8 cp(360,720),ampG(360,720)
      real*8 da1(360,720),da2(360,720),alpha(360,720),phiao(360,720)
      real*8 sigda1(360,720),sigda2(360,720)
      real*8 alphmax,moy,moy2,moyampG,ampGmax, moyampGbox, moybox
      real*8 lat,lon
      integer sauti, sautj, deg_anis, nbox
      character*32 filein,fileanis,fileout
      character*32 fileaout,fileaout2, fileaoutG
      character*32 fileCout
      character*2  nofic,toto,nocoup            
      character*60 ligne

c---------------------------------------------------------
c Version adapted by Alessia 10/2003
c Removed gmt specific scaling: output is now in true G units
c This version assumes whole globe regiostack6 (now standard)
c---------------------------------------------------------
c VERSION adaptee a la regio de l'AUSTRALIE
c On recupere la densite (pour l'aniso crete-crete)dans PREM
c lat,lon adaptees.
c---------------------------------------------------------
c Ce pgm calcule la perturbation de Vs ou de Cr par rapport
c a la vitesse de reference choisie.
c Il ecrit les resultats dans des fichiers : 
c      -lon lat parametres(prof fixee)
c      -Prof lon parametres (lat fixee) pour faire des coupes en prof

c Traitement de l'anisotropie azimutale de la vitesse des ondes S:
c    En sortie de ce pgm, la vitesse isotrope et les parametres
c    Gc et Gs de Montagner et Nataf(1986)
c    On a da1, da2 qui correspondent au coefficient A1 et A2 de la 
c    formule de Smithi et Dahlen ( et non a un rapport) 
c    et l'amplitude crete-crete de G est calculee par :
c    amp(G)=2*ro*betaV*[(2*sqrt((da2**2)+(da1**2))/c)]
c    En fait on laisse tomber le facteur 2 car la barre sera dessinee 2 fois 
c    avec GMT afin de la centrer.
c    on sort 2 fichiers aniso:
c    un coord,angle par rapport a l'est, amplitude G 
c    un coord,angle+180, amplitude G (=demi amplitude crete-crete)
c    amplitude crete-crete= 2* CELUI QUI SORT DE CE PGM!!!
c
c    L'angle est sorti avec un atan2
c    Pour GMT on veut l'angle par rapport a l'EST
c    passage azimut->angle/est : angle=90-azimut
c    On a applique un facteur corr (=corr*2.54) pour que sur 
c    les cartes GMT un trait de 1cm =corr Pa d'amplitude
c    Cette version donne en plus la valeur max et la valeur
c    moyenne de l,anisotropie crete-crete
c---------------------------------------------------------
c   Cette version permet de traiter l'anisotropie azimutale
c   sur des fichiers inverses avec un pas variable. 
c   a l'ecriture des fichiers anisotropes seule une fleche
c   sur 3 (ou sur 6 suivant le pas) est dessinee.
c---------------------------------------------------------
             
      itestaniso2=0
      iout=10
      iaout=22
      iaout2=23
      iaoutG=24
      rad=0.017453293
      radi=0.5/rad

c     y0=90
      y0=0
      x0=-180
c     xpas=2
c     ili=90
c     icol=180
c     y=(y0/xpas)+1
c     x=(x0/xpas)-1

c     calculate values for anisotropy every deg_anis degrees
c     deg_anis=4

      print *,'Reference model: 1066a (1)'
      print *,'                 PREM  (2)'
      print *,'                 none  (3)'
      read(*,*)iref

      print *,'entrer le pas en degres du fichier initial le'
      print *,'nombre de lignes et de colonnes du fichier'
      read(*,*) xpas,ili,icol
      print *,'enter the multiple of the step for outputting anisotropy'
      read (*,*) deg_anis

      print *,'entrer le fichier vitesse '
      read(5,999)filein
      toto=filein(1:2)
      open(7,file=filein)
      print *,'entrer le fichier aniso '
      read(5,999)fileanis
      open(9,file=fileanis)
      print 1000,filein,fileanis
      print *,'entrez le nombre de prof des fichier'
      read(5,*)inprof

      if (inprof .gt. 34) then
       write(*,*)'Need to increase size of depth arrays to at least ',
     &            iprof
       stop
      endif

      print *,'entrez 1 pour faire des coupes en profondeur'
      print *,'(il faut des fichiers Interpolles)'
      print *,'(ne marche que dans le cas de Vs)'
      read(5,*)icoup

      if(icoup.eq.1) then
       print *,'entrez le nombre de coupes que vous voulez faire'
       print *,'entrez les latitudes correspondantes'
       read(5,*)nbcoup,(latcoup(jc),jc=1,nbcoup)

       do 10 jc=1,nbcoup
       iCout(jc)=24+jc
        write(nocoup,'(i2.2)')jc
        fileCout= toto//'couplat'//nocoup//'.xyz'
        open(iCout(jc),file=fileCout)
10     continue
      endif

c-----Lecture de la la densite et de Vs aux prof de 1066a-----
c-----ou a des profondeurs interpollees tout les 25 km--------
      print *,'Voulez-vous lire la densite et Vs aux prof de 1066a (1)'
      print *,'ou dans PREM smoothe et interpolle tout les  25 km (2)?'
      read(5,*)ilec 
      if (ilec.eq.1)then
      open(12,status='old',
     *file='des.1066a.rovs')
c    *file='/home1/eric/these/models/des.1066a.rovs')
      elseif (ilec.eq.2) then
      open(12,status='old',
     *file='/home/alessia/share/des.prem.rovs.int')
      endif
       read(12,*)
      read(12,*)
      read(12,*) npr
      if (npr .gt. 34) then
       write(*,*)'Need to increase size of depth arrays to at least ',
     &            iprof
       stop
      endif
      do 15 i=1,npr
       read(12,*) profi(i),ro(i),vs(i)
c      write(*,*) profi(i),ro(i),vs(i)
       ro(i)=ro(i)*1000
       vs(i)=vs(i)*1000
15    continue
      close(12)

c-------------------------------------------------------------

      do 20 iprof=1,inprof
       itestaniso=0
       write(nofic,'(i2.2)')iprof
       fileout= toto//'prof'//nofic//'.xyz'
       fileaout= 'an'//toto//'prof'//nofic//'.xyz'
       fileaout2= 'an2'//toto//'prof'//nofic//'.xyz'
       fileaoutG= 'Gprof'//nofic//'.xyz'
       open(iout,file=fileout)
       open(iaout,file=fileaout)
       open(iaout2,file=fileaout2)
       open(iaoutG,file=fileaoutG)
c      read tt in both files because need to advance reading 
c      point in both files
       read(7,*)tt
       read(9,*)tt 
       profs(iprof)=0.-tt
       print *,'--------------------------------------'
       print *,'prof',iprof,':',tt
       do 30 i=1,ili
   30  read(7,'(20(f10.6,1x))')(c(i,j),j=1,icol)
       do 40 i=1,ili
   40  read(7,'(20(f10.6,1x))')(er(i,j),j=1,icol)
       do 50 i=1,ili
   50  read(9,'(20(f10.6,1x))')(da1(i,j),j=1,icol)
       do 60 i=1,ili
   60  read(9,'(20(f10.6,1x))')(da2(i,j),j=1,icol)
       do 70 i=1,ili
   70  read(9,'(20(f10.6,1x))')(sigda1(i,j),j=1,icol)
       do 80 i=1,ili
   80  read(9,'(20(f10.6,1x))')(sigda2(i,j),j=1,icol)

       do 85 k=1,npr 
        if (abs(profi(k)-tt).lt.0.001)then
         densite=ro(k)
         beta=vs(k)
        endif
85     continue
    
c     calcul du pourcentage de perturbation par rapport a la valeur
c     moyenne de reference choisie.

       cref=0.
       if(iref.eq.1)then
        write(*,*)'1066a'
        open(8,file='mod1066a')
        read (8,*)ligne    
        write(*,*)ligne    
        read(8,*)ncouch  
        write(*,*)ncouch   
        do  90 icouch=1,ncouch 
         read(8,*)prof(icouch),vref(icouch) 
c        write(*,*)vref(icouch)
         if(prof(icouch).eq.tt) cref=vref(icouch)
 90     continue
        print*,'la valeur Vs de reference de 1066a est',cref
        close(8)
       else if (iref.eq.2) then
        write(*,*)'PREM'
        open(8,status='old',
     *  file='/home/alessia/share/des.prem.rovs.int')
        read(8,*)ligne
        read(8,*)ligne
        read(8,*) ncouch
        do 95 icouch=1,ncouch
         read(8,*) prof(icouch),bid,vref(icouch)
c        write(*,*) prof(icouch),vref(icouch)
         if(prof(icouch).eq.tt) cref=vref(icouch)
95      continue
        print*,'la valeur Vs de reference de PREM est',cref
        close(8)
       else
        do 100 i=1,ili
        do 100 j=1,icol
100     cref=cref+c(i,j)
        cref=cref/ili/icol 
        print*,'la valeur moyenne de Vs (reference) est ',cref
       endif
       print*,'la valeur de reference choisie pour Vs est',cref

       cmax=0
       cmin=0
       do 110 i=1,ili
       do 110 j=1,icol
        cp(i,j)=(c(i,j)-cref)/cref*100
        if(cp(i,j).gt.cmax)cmax=cp(i,j)
        if(cp(i,j).lt.cmin)cmin=cp(i,j)
110    continue
        print*,'perturbation Vs max et min',cmax,cmin

c     Conversion des indices colonnes en latitude et longitude
      moy=0.
      moy2=0.
      imoy2=0
      moyampG=0.
      alphmax=0.
      ampGmax=0.
      sauti=0.
      sautj=0.
      moyampGbox=0.
      moybox=0.
      nbox=0
      do 130 i=1,ili
       sauti=sauti+1
       if (sauti.eq.deg_anis)sauti=0
c      lat=(y-i)*xpas
       lat=y0+(i-1)*xpas+xpas/2.
       lat=90.-lat 
       do 130 j=1,icol
        sautj=sautj+1
        if (sautj.eq.deg_anis)sautj=0
c       lon=(x+j)*xpas   
        lon=x0+(j-1)*xpas+xpas/2.
        write(iout,'(f8.3,1x,f8.3,2x,f8.3)')lat,lon,cp(i,j)
        if(da1(i,j).ne.0..and.da2(i,j).ne.0.)then 
          itestaniso=itestaniso+1      
          alpha(i,j)=sqrt((da1(i,j)**2)+(da2(i,j)**2))
c         da1 et da2 sont en km/s, on multiplie par 1000
c         pour les passer en m/s
          ampG(i,j)=2*densite*beta*1000*alpha(i,j)
          alpha(i,j)= 100*(alpha(i,j)/c(i,j))
c         Extraction des valeurs max et moyennes de l'anisotropie
c         crete-crete
          moy=moy+2*alpha(i,j)
          moyampG=moyampG+2*ampG(i,j)
c         if we are within the box then add to box mean G amplitude
          if ((lat .le. 60 .and. lat .ge. -60) .and. 
     &       ((lon .ge. 120 .and. lon .le. 180) .or. 
     &        (lon .ge. -180 .and. lon .le. -70)) ) then
            moyampGbox=moyampGbox+2*ampG(i,j)
            moybox=moybox+2*alpha(i,j)
            nbox=nbox+1
          end if
          if(alpha(i,j).gt.alphmax)then
            alphmax=alpha(i,j)
            ampGmax=ampG(i,j)
            latmax=lat
            lonmax=lon
          endif
          phiao(i,j)= radi*atan2(da2(i,j),da1(i,j))
c         Phiao est l'angle par rapport au nord. Dans GMT
c         psxy plot les vecteurs par rapport a l'Est.
c         On convertit
          phiao(i,j)=90.-phiao(i,j)
          write(iaoutG,'(4(f8.3,1x))')lat,lon,2*ampG(i,j)/10**9
          if(sauti.eq.0.and.sautj.eq.0)then
            write(iaout,'(4(f8.3,1x))')lat,lon,phiao(i,j),alpha(i,j)
            write(iaout2,'(4(f8.3,1x))')lat,lon,(phiao(i,j)+180)
     *                                ,alpha(i,j)
           endif
        endif
130     continue
        moy=moy/(ili*icol)
        moyampG=moyampG/(ili*icol)
        if (nbox .gt. 0) then 
          write(*,*) 'Nbox: ', nbox
          moyampGbox=moyampGbox/nbox
          moybox=moybox/nbox
          sigampGbox=0.
          sigampbox=0.
          do 160 i=1,ili
c           lat=(y-i)*xpas
            lat=y0+(i-1)*xpas+xpas/2. 
            lat=90.-lat 
            do 161 j=1,icol
c             lon=(x+j)*xpas   
              lon=x0+(j-1)*xpas+xpas/2.
              if ((lat .le. 60 .and. lat .ge. -60) .and. 
     &           ((lon .ge. 120 .and. lon .le. 180) .or. 
     &            (lon .ge. -180 .and. lon .le. -70)) ) then
                sigampGbox=(2*ampG(i,j)-moyampGbox)**2
                sigampbox=(2*alpha(i,j)-moybox)**2
              endif
161         continue
160       continue
          sigampGbox=sqrt(sigampGbox/nbox)
          sigampbox=sqrt(sigampbox/nbox)
        end if
        write(*,*)'Pourcentage moyen anisotroopie crete-crete:'
        write(*,*)moy
        write(*,*)'Pourcentage moyen anisotroopie crete-crete (box):'
        write(*,*)moybox, '(', sigampbox,')'
        write(*,*)'Amplitude moyenne de Gc (crete-crete) en GPas:'
        write(*,*)moyampG/10**9
        write(*,*)'Amplitude moyenne de Gc (crete-crete) en GPas (box):'
        write(*,*)moyampGbox/10**9, '(', sigampGbox/10**9,')'
        write(*,*)'Pourcentage max anisotroopie crete-crete:'
        write(*,*)'plus latitude et longitude du point:'
        write(*,'(f8.3,2i5)')2*alphmax,latmax,lonmax
        write(*,*)'Amplitude max Gc crete-crete en GPas:'
        write(*,*)'pour ce meme point:'
        write(*,*)2*ampGmax/10**9
        close(iout)
        close(iaout)
        close(iaout2)
        close(iaoutG)

c     On ecrit dans les fichiers contenant les coupes.
c only if we are doing cross-sections!

        if (icoup .eq. 1 ) then
          do 140 jc=1,nbcoup
          do 140 i=1,ili
c           lat=(y-i)*xpas
            lat=y0+(i-1)*xpas+xpas/2. 
            lat=90.-lat 
            do 140 j=1,icol
c             lon=(x+j)*xpas
              lon=x0+(j-1)*xpas+xpas/2.
              if(icoup.eq.1.and.latcoup(jc).eq.lat)then
                write(iCout(jc),1200)profs(iprof),lon,cp(i,j)
              endif
140       continue
        endif
1200    format(f8.2,1x,f8.2,2x,f8.3)

        if(itestaniso.eq.(ili*icol))then
         itestaniso2=itestaniso2+1
        elseif(itestaniso.eq.0)then
         itestaniso2=itestaniso2-1
        endif

20    continue

      if(itestaniso2.eq.inprof)then
       write(*,*)'On traite anisotropie azimutale'
       write(*,*)'sur toute les profondeurs'
      elseif(itestaniso2.eq.(-1*inprof))then
       write(*,*)'On ne traite pas anisotropie azimutale'
       write(*,*)'pour toute les profondeurs'
      else
       write(*,*)'Certaines valeurs de A1 et A2 sont egale'
       write(*,*)'a zero. Ces points ne seront pas representes'
      endif 

c    On ferme les fichiers contenant les coupes

      if(icoup .eq. 1) then
      do 150 jc=1,nbcoup
       close(iCout(jc))
150   continue
      endif

      close(7)
      close(9)
999   format(a32)
1000  format(2a10)
1300  format(i5,f7.2,2f8.3)
      stop
      end
