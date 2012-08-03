      parameter(NPROFMAX=50,NBTRAMAX=40000)
      parameter(NPTSMAX=18000,NCOUPMAX=120)
      real latc(NCOUPMAX) 
      real lat1(NCOUPMAX), lon1(NCOUPMAX)
      real lat2(NCOUPMAX), lon2(NCOUPMAX)
      real profs(NPROFMAX),prof(NPROFMAX),cref(NPROFMAX),depth(NPROFMAX)
      real vs(NPROFMAX),ro(NPROFMAX)
      real densite(NPROFMAX), beta(NPROFMAX)
      real*8 c(360,720), cp(NPROFMAX,360,720),ampG(NPROFMAX,360,720)
      real*8 da1(360,720),da2(360,720),alpha(NPROFMAX,360,720)
      real*8 phiao(NPROFMAX,360,720)
      real*8 alphmax,moy,moy2,moyampG,ampGmax, moyampGbox, moybox
      real*8 lat,lon
      integer sauti, sautj, deg_anis, nbox
      character*32 filein,fileanis,fileout
      character*32 fileaout,fileaout2, fileaoutG
      character*32 fileCout
      character*2  nofic,toto,nocoup            
      character*60 ligne
      logical want_depth

c---------------------------------------------------------
c Adapted by Alessia 2/2004
c New version of regiostack6 now no longer produces error maps, so
c removed reading of error maps from input.
c Added great circle path cuts for output.
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

      y0=0
      x0=-180

c     -----------------------------------------------------------------
c     basic input
c     -----------------------------------------------------------------
      print *,'Reference model: PREM  (1)'
      print *,'                 none  (2)'
      read(*,*)iref

      print *,'entrer le pas en degres du fichier initial le'
      print *,'nombre de lignes et de colonnes du fichier'
      read(*,*) xpas,ili,icol
      print *,'enter the multiple of the step for outputting anisotropy'
      read (*,*) deg_anis

      print *,'entrer le fichier vitesse '
      read(5,*)filein
      toto=filein(1:2)
      open(7,file=filein)
      print *,'entrer le fichier aniso '
      read(5,*)fileanis
      open(9,file=fileanis)
      print *,filein,fileanis
      print *,'entrez le nombre de prof des fichier'
      read(5,*)inprof

      if (inprof .gt. NPROFMAX) then
       write(*,*)'Need to increase size of depth arrays to at least ',
     &            inprof
       stop
      endif

c     -----------------------------------------------------------------
c     set up output stuff
c     -----------------------------------------------------------------
      print *,'Number of depth layers (max ', NPROFMAX,')'
      print *,'and depthlist '
      read(*,*) nbdepth,(depth(ij),ij=1,nbdepth)

      print *,'Number of latitude cross-sections (max ', NCOUPMAX,')'
      print *,'and lat of the cross-sections '
      write(*,*)'(e.g. 2,15,17)'
      read(*,*)nbcouplat,(latc(ij),ij=1,nbcouplat)

      print *,'Number of great circle path cross-sections ',
     *        '(max ', NCOUPMAX,')'
      read(*,*)nbcoup
      do 4 i=1,nbcoup
        print *,'Coord of 2 points to define the minor arc along which'
        write(*,*)'the cross-section ',i,' is plotted '
        write(*,*)'(lat1,lon1,lat2,lon2)'
        write(*,*)'(e.g. 25,-15,52,25)'
        read(*,*)lat1(i),lon1(i),lat2(i),lon2(i)
4     continue



c     -----------------------------------------------------------------
c     read reference model (if nothing else, for later use in Gamp)
c     -----------------------------------------------------------------
      open(8,status='old',
     *       file='/home/alessia/share/des.prem.rovs.int')
      read(8,*)ligne
      read(8,*)ligne
      read(8,*) ncouch
      do 95 i=1,ncouch
        read(8,*) prof(i),ro(i),vs(i)
c       ro and vs will later be used to calculate amplitude of G
        ro(i) = ro(i)*1000
        vs(i) = vs(i)*1000
95    continue
      close(8)
       
c     -----------------------------------------------------------------
c     START LOOP OVER DEPTHS FOR READING FILES
c     -----------------------------------------------------------------

      do 20 iprof=1,inprof
c       read tt in both files because need to advance reading 
c       point in both files
        read(7,*)tt
        read(9,*)tt 
c       get the depth from the file
        profs(iprof)=tt 
        print *,'--------------------------------------'
        do 30 i=1,ili
30        read(7,'(20(f10.6,1x))')(c(i,j),j=1,icol)
        do 50 i=1,ili
50        read(9,'(20(f10.6,1x))')(da1(i,j),j=1,icol)
        do 60 i=1,ili
60        read(9,'(20(f10.6,1x))')(da2(i,j),j=1,icol)

c       find the correct density and beta for this depth
        do 85 k=1,ncouch 
          if (abs(prof(k)-tt).lt.0.001)then
           densite(iprof)=ro(k)
           beta(iprof)=vs(k)
          endif
85      continue

c       find reference velocity
        cref(iprof)=0
        if(iref.eq.1) then 
          cref(iprof) = beta(iprof)/1000.0
c         print*,'la valeur Vs de reference de PREM est',cref(iprof)
        else
          do 100 i=1,ili
            do 100 j=1,icol
100           cref(iprof)=cref(iprof)+c(i,j)
          cref(iprof)=cref(iprof)/ili/icol 
c         print*,'la valeur moyenne de Vs (reference) est ',cref(iprof)
        endif
c       print*,'la valeur de reference choisie pour Vs est',cref(iprof)
        print *,'prof',iprof,':',tt, '  ref_vel : ', cref(iprof)


c       calcul du pourcentage de perturbation par rapport a la valeur
c       moyenne de reference choisie.
        cmax=0
        cmin=0
        do 110 i=1,ili
        do 110 j=1,icol
          cp(iprof,i,j)=(c(i,j)-cref(iprof))/cref(iprof)*100
          if(cp(iprof,i,j).gt.cmax)cmax=cp(iprof,i,j)
          if(cp(iprof,i,j).lt.cmin)cmin=cp(iprof,i,j)
110     continue
c       print*,'perturbation Vs max et min',cmax,cmin

c       calculation of anisotropic parameters
        do 130 i=1,ili
          do 130 j=1,icol
          if(da1(i,j).ne.0..and.da2(i,j).ne.0.)then 
            alpha(iprof,i,j)=sqrt((da1(i,j)**2)+(da2(i,j)**2))
c           da1 et da2 sont en km/s, on multiplie par 1000
c           pour les passer en m/s
            ampG(iprof,i,j)=2*densite(iprof)*beta(iprof)*1000*
     *                      alpha(iprof,i,j)
c           turn alpha into a percentage
            alpha(iprof,i,j)= 100*(alpha(iprof,i,j)/c(i,j))
            phiao(iprof,i,j)= radi*atan2(da2(i,j),da1(i,j))
c           Phiao est l'angle par rapport au nord. Dans GMT
c           psxy plot les vecteurs par rapport a l'Est.
c           On convertit
            phiao(iprof,i,j)=90.-phiao(iprof,i,j)
          endif
130     continue

20    continue

c     -----------------------------------------------------------------
c     END LOOP OVER DEPTHS ON INPUT
c     -----------------------------------------------------------------

      close(7)
      close(9)

c     -----------------------------------------------------------------
c     START OUTPUT STUFF
c     -----------------------------------------------------------------

c     -----------------------------------------------------------------
c     do all the horizontal cuts first
c     -----------------------------------------------------------------
      do 200 iprof=1,inprof
        want_depth=.false.
c       do we want this depth?
        do 201, k=1,nbdepth
          if (abs(depth(k)-profs(iprof)).lt.0.001) want_depth=.true.
201     continue
        if(want_depth) then
c         set up output files 
          write(nofic,'(i2.2)')iprof
          fileout= toto//'prof'//nofic//'.xyz'
          fileaout= 'an'//toto//'prof'//nofic//'.xyz'
          fileaout2= 'an2'//toto//'prof'//nofic//'.xyz'
          fileaoutG= 'Gprof'//nofic//'.xyz'
          open(iout,file=fileout)
          open(iaout,file=fileaout)
          open(iaout2,file=fileaout2)
          open(iaoutG,file=fileaoutG)
          do 230 i=1,ili
            sauti=sauti+1
            if (sauti.eq.deg_anis)sauti=0
            lat=90-(y0+(i-1)*xpas)
            do 230 j=1,icol
              sautj=sautj+1
              if (sautj.eq.deg_anis)sautj=0
              lon=x0+(j-1)*xpas
              write(iout,'(f8.3,1x,f8.3,2x,f8.3)')lat,lon,cp(iprof,i,j)
              write(iaoutG,'(4(f8.3,1x))')lat,lon,
     *                                    2*ampG(iprof,i,j)/10**9
              if(sauti.eq.0.and.sautj.eq.0)then
                write(iaout,'(4(f8.3,1x))')lat,lon,phiao(iprof,i,j),
     *                                     alpha(iprof,i,j)
                write(iaout2,'(4(f8.3,1x))')lat,lon,phiao(iprof,i,j)+180
     *                                     ,alpha(iprof,i,j)
              endif
230       continue
          close(iout)
          close(iaout)
          close(iaout2)
          close(iaoutG)
        endif
200   continue

c     -----------------------------------------------------------------
c     the vertical cross-sections (consant latitude)
c     -----------------------------------------------------------------
      do 300 ic=1,nbcouplat

        write(nofic,'(i2.2)')ic
        fileout= toto//'lat'//nofic//'.xyz'
        open(iout,file=fileout)

        i=int((90-latc(ic)-y0)/xpas)+1
        do 311 iprof=1,inprof
          do 311, j=1,icol
            lon=x0+(j-1)*xpas
            write(iout,1200) -profs(iprof),lon,cp(iprof,i,j)
311     continue

        close(iout)
300   continue


c     -----------------------------------------------------------------
c     the vertical cross-sections (great circle path)
c     -----------------------------------------------------------------
      do 400 ic=1,nbcoup
        write(nofic,'(i2.2)')ic
        fileout= toto//'gcp'//nofic//'.xyz'
        open(iout,file=fileout)

c       initialise path geometry
        alat1=lat1(ic)
        alon1=lon1(ic)
        alat2=lat2(ic)
        alon2=lon2(ic)
        call mygrt(alat1,alon1,alat2,alon2,disd,az12,az21,gc)
        delta=disd
        azes=az12
        do while (azes.lt.0.) 
          azes=azes+360.
        enddo
        do while (azes.gt.360.) 
          azes=azes-360.
        enddo
        dpas= xpas/2.
        npas= int(delta/dpas)
        npat= npas+1

c       Boucle sur chaque point de l'azimuth epicentre-station
        do 410 jcle=1,npat
          dstcle=dpas*float(jcle-1)*111.195
          alat1=lat1(ic)
          alon1=lon1(ic)
          call mygds(alat1,alon1,azes,dstcle,blat,blon)
          i=int((90-(blat-y0))/xpas)+1
          j=int((blon-x0)/xpas)+1
          do 411 iprof=1,inprof
            write(iout,*)-profs(iprof),jcle, cp(iprof,i,j)
411       continue
410     continue

        close(iout)
400   continue

c     -----------------------------------------------------------------
c     The reference model
c     -----------------------------------------------------------------
      open(iout,file='ref-model.dat')
      do 500 iprof=1,inprof
        write(iout,'(f8.2,1x,f8.3)') -profs(iprof), cref(iprof)
500   continue
      close(iout)

c     -----------------------------------------------------------------
c     END OUTPUT STUFF
c     -----------------------------------------------------------------

1200  format(f8.2,1x,f8.2,2x,f8.3)
      end
