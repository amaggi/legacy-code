      real*8 latcoup(10),iCout(10),iCerout(10)
      real*8 profs(30),c(360,720),er(360,720),prof(30),vref(30)
      real*8 profi(30),vs(30),ro(30)
      real*8 cp(360,720),ampG(360,720)
      real*8 da1(360,720),da2(360,720),alpha(360,720),phiao(360,720)
      real*8 sigda1(360,720),sigda2(360,720)
      real*8 alphmax,moy,moy2,moyampG,ampGmax 
      real*8 lat,lon
      integer sauti, sautj
      character*32 filein,fileanis,fileout
      character*32 fileaout,fileaout2,fileerout
      character*32 fileCout,fileCerout
      character*2  nofic,toto,nocoup            
      character*60 ligne

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
      ierout=11
      iaout=12
      iaout2=13
      rad=0.017453293
      radi=0.5/rad

      print *,'On prends 1066a comme modele de reference?(1=oui)'
      read(*,*)iref

      print *,'entrer latitude et longitude initiales'
      print *,'Pour l\'Australie, 20. 90. '
      print *,'Pour la corne de l\'Afrique, 55. 0. '
      read(*,*) y0,x0
 
      print *,'entrer le pas en degres du fichier initial le'
      print *,'nombre de lignes et de colonnes du fichier'
      read(*,*) xpas,ili,icol
      y=(y0/xpas)+1
      x=(x0/xpas)-1
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
      print *,'Entrez corr, facteur de correction pour le plot de'
      print *,'l\'anisotropie(corr=5->1cm = 5% d\'aniso sur les cartes)'
      read(*,*)corr     

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
       iCerout(jc)=34+jc
        write(nocoup,'(i2.2)')jc
        fileCout= toto//'couplat'//nocoup//'.xyz'
        fileCerout= 'Er'//toto//'couplat'//nocoup//'.xyz'
        open(iCout(jc),file=fileCout)
        open(iCerout(jc),file=fileCerout)
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
      do 15 i=1,npr
       read(12,*) profi(i),ro(i),vs(i)
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
       fileerout= 'Er'//toto//'prof'//nofic//'.xyz'
       open(iout,file=fileout)
       open(iaout,file=fileaout)
       open(iaout2,file=fileaout2)
       open(ierout,file=fileerout)
       read(9,*)tt
       read(7,*)tt
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
       if ((profi(k)-tt).lt.0.001)then
       densite=ro(k)
       beta=vs(k)
       endif
 85    continue
    
c     calcul du pourcentage de perturbation par rapport a la valeur
c     moyenne de reference choisie.

       cref=0.
       if(iref.eq.1)then
        write(*,*)'bonjour'
        open(8,file='mod1066a')
        read (8,*)ligne    
        write(*,*)ligne    
        read(8,*)ncouch  
        write(*,*)ncouch   
        do  90 icouch=1,ncouch 
         read(8,*)prof(icouch),vref(icouch) 
         write(*,*)vref(icouch)
         if(prof(icouch).eq.tt) cref=vref(icouch)
 90     continue
        print*,'la valeur Vs de reference de 1066a est',cref
        close(8)
       else
        do 100 i=1,ili
        do 100 j=1,icol
 100    cref=cref+c(i,j)
        cref=cref/ili/icol 
        print*,'la valeur moyenne  de Vs est',cref
c Valeurs de refence de PREM
c      if (tt.eq.75.)cref=4.409
c      if (tt.eq.100.)cref=4.413
c      if (tt.eq.125.)cref=4.420
c       if (tt.eq.150.)cref=4.425
c       if (tt.eq.175.)cref=4.432
c       if (tt.eq.200.)cref=4.436
c       if (tt.eq.225.)cref=4.619
c       if (tt.eq.250.)cref=4.663
c       if (tt.eq.275.)cref=4.683
c       if (tt.eq.300.)cref=4.700
c       if (tt.eq.325.)cref=4.718
c       if (tt.eq.350.)cref=4.730
c       if (tt.eq.375.)cref=4.750
c       if (tt.eq.400.)cref=4.849


c       if(tt.le.200)cref=4.5
c       if(tt.eq.60)cref=4.00
c       if(tt.eq.100)cref=4.10
       if (tt.eq.50.)cref=4.402
        print*,'la valeur de reference choisie pour Vs est',cref
       endif

       cmax=0
       cmin=0
       do 110 i=1,ili
       do 110 j=1,icol
        cp(i,j)=(c(i,j)-cref)/cref*100
        if(cp(i,j).gt.cmax)cmax=cp(i,j)
        if(cp(i,j).lt.cmin)cmin=cp(i,j)
c       if(cp(i,j).gt.10.5)cp(i,j)=10.5
c       if(cp(i,j).lt.-10.5)cp(i,j)=-10.5
 110	 continue
        print*,'perturbation Vs max et min',cmax,cmin

c     Conversion des indices colonnes en latitude et longitude
       moy=0.
       moy2=0.
       imoy2=0
       moyampG=0.
c      alphmax=0 car tjrs positif
       alphmax=0.
       ampGmax=0.
       sauti=0.
       sautj=0.
       do 130 i=1,ili
        sauti=sauti+1
        if (sauti.eq.4)sauti=0
        lat=(y-i)*xpas
       do 130 j=1,icol
         sautj=sautj+1
         if (sautj.eq.4)sautj=0
         lon=(x+j)*xpas   
         write(iout,'(f8.3,1x,f8.3,2x,f8.3)')lat,lon,cp(i,j)
         write(ierout,'(f8.3,1x,f8.3,2x,f8.3)')lat,lon,er(i,j) 
         if(da1(i,j).ne.0..and.da2(i,j).ne.0.)then 
          itestaniso=itestaniso+1      
          alpha(i,j)=sqrt((da1(i,j)**2)+(da2(i,j)**2))
c    da1 et da2 sont en km/s, on multiplie par 1000
c    pour les passer en m/s
          ampG(i,j)=2*densite*beta*1000*alpha(i,j)
          alpha(i,j)= 100*(alpha(i,j)/c(i,j))
c         bid=2*alpha(i,j)
c   Extraction des valeurs max et moyennes de l'anisotropie
c   crete-crete
          moy=moy+2*alpha(i,j)
c         if(lon.ge.30.and.lon.le.80.and.lat.ge.-10.
c    *and.lat.le.40) then
c           moy2=moy2+2*alpha(i,j)
c           imoy2=imoy2+1
c         endif
          moyampG=moyampG+2*ampG(i,j)
           if(alpha(i,j).gt.alphmax)then
             alphmax=alpha(i,j)
             ampGmax=ampG(i,j)
             latmax=lat
             lonmax=lon
           endif
c   Facteur d'echelle tel que sur les cartes corr%aniso=1cm
c   1 unite est representee par 2.54cm (1pouce) par gmt
          alpha(i,j)=alpha(i,j)/(2.54*corr)
          phiao(i,j)= radi*atan2(da2(i,j),da1(i,j))
c   Phiao est l'angle par rapport au nord. Dans GMT
c   psxy plot les vecteurs par rapport a l'Est.
c   On convertit
          phiao(i,j)=90.-phiao(i,j)
           if(sauti.eq.0.and.sautj.eq.0)then
c            write(iaout,'(2i5,2f8.3)')lat,lon,phiao(i,j),alpha(i,j)
c     write(iaout2,'(2i5,2f8.3)')lat,lon,(phiao(i,j)+180),alpha(i,j)
             write(iaout,'(4(f8.3,1x))')lat,lon,phiao(i,j),alpha(i,j)
             write(iaout2,'(4(f8.3,1x))')lat,lon,(phiao(i,j)+180)
     *,alpha(i,j)
c            write(*,*)lat,lon,da1(i,j),da2(i,j),phiao(i,j)
           endif
         endif
 130    continue
      moy=moy/(ili*icol)
c     moy2=moy2/imoy2
      moyampG=moyampG/(ili*icol)
      write(*,*)'Pourcentage moyen d\'anisotroopie crete-crete:'
      write(*,'(f8.3)')moy
c     write(*,*)'Pourcentage moyen d\'anisotroopie crete-crete'
c     write(*,*)'pour la zone (-10S<lat<40; 30E<lon<80E)'
c     write(*,'(f8.3)')moy2
      write(*,*)'Amplitude moyenne de Gc (crete-crete) en GigaPascals:'
      write(*,*)moyampG/10**9
      write(*,*)'Pourcentage max d\'anisotroopie crete-crete:'
      write(*,*)'plus latitude et longitude du point:'
      write(*,'(f8.3,2i5)')2*alphmax,latmax,lonmax
      write(*,*)'Amplitude max Gc crete-crete en GigaPascals:'
      write(*,*)'pour ce meme point:'
      write(*,*)2*ampGmax/10**9
      close(iout)
      close(iaout)
      close(iaout2)
      close(ierout)

c     On ecrit dans les fichiers contenant les coupes.

      do 140 jc=1,nbcoup
      do 140 i=1,ili
        lat=(y-i)*xpas
      do 140 j=1,icol
        lon=(x+j)*xpas
       if(icoup.eq.1.and.latcoup(jc).eq.lat)then
          write(iCout(jc),1200)profs(iprof),lon,cp(i,j)
          write(iCerout(jc),1200)profs(iprof),lon,er(i,j)
        endif
140   continue
1200  format(f8.2,1x,f8.2,2x,f8.3)
      if(itestaniso.eq.(ili*icol))then
      itestaniso2=itestaniso2+1
      elseif(itestaniso.eq.0)then
      itestaniso2=itestaniso2-1
      endif

20    continue

      if(itestaniso2.eq.inprof)then
       write(*,*)'On traite l\'anisotropie azimutale'
       write(*,*)'sur toute les profondeurs'
      elseif(itestaniso2.eq.(-1*inprof))then
       write(*,*)'On ne traite pas l\'anisotropie azimutale'
       write(*,*)'pour toute les profondeurs'
      else
       write(*,*)'Certaines valeurs de A1 et A2 sont egale'
       write(*,*)'a zero. Ces points ne seront pas representes'
      endif 

c    On ferme les fichiers contenant les coupes

      do 150 jc=1,nbcoup
       close(iCout(jc))
       close(iCerout(jc))
150   continue

      close(7)
      close(9)
999   format(a32)
1000  format(2a10)
1300  format(i5,f7.2,2f8.3)
      stop
      end
