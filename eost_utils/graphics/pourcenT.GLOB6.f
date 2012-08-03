      parameter(NCOUP=10,NPROF=34,NLAT=360,NLON=720)
      real*8 latcoup(NCOUP),iCout(NCOUP)
      real*8 prof(NPROF),vref(NPROF),vs(NPROF)
      real*8 profs(NPROF),profi(NPROF),ro(NPROF)
      real*8 c(NLAT,NLON),er(NLAT,NLON)
      real*8 cp(NLAT,NLON),ampG(NLAT,NLON),phiao(NLAT,NLON)
      real*8 da1(NLAT,NLON),da2(NLAT,NLON),alpha(NLAT,NLON)
      real*8 sigda1(NLAT,NLON),sigda2(NLAT,NLON)
      real*8 alphmax,moy,moy2,moyampG,ampGmax 
      real*8 lat,lon
      integer sauti, sautj
      character*32 filein,fileanis,fileout
      character*32 fileaout,fileaout2
      character*32 fileCout
      character*2  nofic,toto,nocoup            
      character*60 ligne

c---------------------------------------------------------
c 04/02/04 Code by Eric Debayle
c VERSION adaptee a regiostack6 (fichiers erreurs supprimes) 
c On recupere la densite (pour l'aniso crete-crete) dans PREM
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
c    formule de Smith et Dahlen ( et non a un rapport) 
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
      iaout=12
      iaout2=13
      rad=0.017453293
      radi=0.5/rad

      print *,'On prends 1066a comme modele de reference?(1=oui)'
      read(*,*)iref

      y0=0.
      x0=-180.
 
      print *,'entrer le pas en degres du fichier initial le'
      print *,'nombre de lignes et de colonnes du fichier'
      read(*,*) xpas,ili,icol
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

      nbcoup=0
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
      elseif (ilec.eq.2) then
      open(12,status='old',
     *file='des.prem.rovs.int')
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
       open(iout,file=fileout)
       open(iaout,file=fileaout)
       open(iaout2,file=fileaout2)
       read(9,*)tt
       read(7,*)tt
       profs(iprof)=0.-tt
       print *,'--------------------------------------'
       print *,'prof',iprof,':',tt
       do 30 i=1,ili
   30  read(7,'(20(f10.6,1x))')(c(i,j),j=1,icol)
       do 50 i=1,ili
   50  read(9,'(20(f10.6,1x))')(da1(i,j),j=1,icol)
       do 60 i=1,ili
   60  read(9,'(20(f10.6,1x))')(da2(i,j),j=1,icol)

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
c Valeurs de reference de PREM
c          if (tt.eq.40.)cref=4.400
c          if (tt.eq.50.)cref=4.402
c          if (tt.eq.75.)cref=4.409
c          if (tt.eq.100.)cref=4.413
c          if (tt.eq.125.)cref=4.420
c           if (tt.eq.150.)cref=4.425
c           if (tt.eq.175.)cref=4.432
c           if (tt.eq.200.)cref=4.436
c           if (tt.eq.225.)cref=4.619
c           if (tt.eq.250.)cref=4.663
c           if (tt.eq.275.)cref=4.683
c           if (tt.eq.300.)cref=4.700
c           if (tt.eq.325.)cref=4.718
c           if (tt.eq.350.)cref=4.730
c           if (tt.eq.375.)cref=4.750
c           if (tt.eq.400.)cref=4.849
c           if (tt.eq.450.)cref=5.078
c           if (tt.eq.500.)cref=5.224
c           if (tt.eq.550.)cref=5.370
c           if (tt.eq.600.)cref=5.516
c           if (tt.eq.650.)cref=5.552
c           if (tt.eq.700.)cref=6.019

        print*,'la valeur de reference choisie pour Vs est',cref
       endif

       cmax=0
       cmin=0
       do 110 i=1,ili
       do 110 j=1,icol
        cp(i,j)=(c(i,j)-cref)/cref*100
        if(cp(i,j).gt.cmax)cmax=cp(i,j)
        if(cp(i,j).lt.cmin)cmin=cp(i,j)
 110   continue
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
       fs=2
       do 130 i=1,ili
        sautj=0.
        sauti=sauti+1
        if (sauti.eq.2*fs)sauti=0
        lat=y0+(i-1)*xpas+xpas/2.
        lat=90.-lat 
       do 130 j=1,icol+5
         sautj=sautj+1
         if (sautj.eq.2*fs.and.abs(lat).ge.0.and.abs(lat).lt.50)sautj=0
         if (sautj.eq.3*fs.and.abs(lat).ge.50.and.abs(lat).lt.70)sautj=0
         if (sautj.eq.4*fs.and.abs(lat).ge.70.and.abs(lat).lt.80)sautj=0
         if (sautj.eq.6*fs.and.abs(lat).ge.80.and.abs(lat).le.90)sautj=0
         lon=x0+(j-1)*xpas+xpas/2.
         if (j.gt.icol)then
            jj=j-icol
         else
            jj=j
         endif
         write(iout,'(f8.3,1x,f8.3,2x,f8.3)')
     *lat,lon,cp(i,jj)
         if(da1(i,jj).ne.0..and.da2(i,jj).ne.0.)then 
          itestaniso=itestaniso+1      
          alpha(i,jj)=sqrt((da1(i,jj)**2)+(da2(i,jj)**2))
c    da1 et da2 sont en km/s, on multiplie par 1000
c    pour les passer en m/s
          ampG(i,jj)=2*densite*beta*1000*alpha(i,jj)
          alpha(i,jj)= 100*(alpha(i,jj)/c(i,jj))
c   Extraction des valeurs max et moyennes de l'anisotropie
c   crete-crete
          moy=moy+2*alpha(i,jj)
c         if(lon.ge.30.and.lon.le.80.and.lat.ge.-10.
c    *and.lat.le.40) then
c           moy2=moy2+2*alpha(i,jj)
c           imoy2=imoy2+1
c         endif
          moyampG=moyampG+2*ampG(i,jj)
           if(alpha(i,jj).gt.alphmax)then
             alphmax=alpha(i,jj)
             ampGmax=ampG(i,jj)
             latmax=lat
             lonmax=lon
           endif
c   Facteur d'echelle tel que sur les cartes corr%aniso=1cm
c   1 unite est representee par 2.54cm (1pouce) par gmt
          alpha(i,jj)=alpha(i,jj)/(2.54*corr)
          phiao(i,jj)= radi*atan2(da2(i,jj),da1(i,jj))
c   Phiao est l'angle par rapport au nord. Dans GMT
c   psxy plot les vecteurs par rapport a l'Est.
c   On convertit
          phiao(i,jj)=90.-phiao(i,jj)
           if(sauti.eq.0.and.sautj.eq.0)then
             write(iaout,'(4(f8.3,1x))')lat,lon,phiao(i,jj),alpha(i,jj)
             write(iaout2,'(4(f8.3,1x))')lat,lon,(phiao(i,jj)+180)
     *,alpha(i,jj)
           endif
         endif
 130    continue
      moy=moy/(ili*icol)
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
      write(*,*)'pour ce meme point :'
      write(*,*)2*ampGmax/10**9
      close(iout)
      close(iaout)
      close(iaout2)

c     On ecrit dans les fichiers contenant les coupes.

      do 140 jc=1,nbcoup
      do 140 i=1,ili
        lat=y0+(i-1)*xpas+xpas/2. 
        lat=90.-lat 
c       lat=(y-i)*xpas
      do 140 j=1,icol+5
c       lon=(x+j)*xpas
        lon=x0+(j-1)*xpas+xpas/2.
        if (j.gt.icol)then
            jj=j-icol
         else
            jj=j
         endif
       if(icoup.eq.1.and.latcoup(jc).eq.lat)then
          write(iCout(jc),1200)profs(iprof),lon,cp(i,jj)
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
150   continue

      close(7)
      close(9)
999   format(a32)
1000  format(2a10)
1300  format(i5,f7.2,2f8.3)
      stop
      end
