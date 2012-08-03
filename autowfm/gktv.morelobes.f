      parameter(NPERMAX=30,NPARMAX=7,NCOUCHMAX=35,NPMAX=1024,NPMAX2=513)
c     parameter(NFILTRMAX=5,NDATAMAX=5,NMODMAX=10)

c     modif eric Oct 97 inversion plusieurs lobes : NDATAMAX=10
c     parameter(NFILTRMAX=5,NDATAMAX=10,NMODMAX=10)
c     modif eric 16 Oct 2000 inversion plus de 5 periodes : NFILTRMAX=10
      parameter(NFILTRMAX=10,NDATAMAX=10,NMODMAX=10)

      dimension pper(NFILTRMAX),ppar(NFILTRMAX),ndata(NFILTRMAX),
     *          idata(NFILTRMAX,NDATAMAX),nf1(NFILTRMAX),
     *          nband(NFILTRMAX), refr(NPMAX),refi(NPMAX),
     *          syntr(NPMAX,NMODMAX), synti(NPMAX,NMODMAX)
      dimension str(NFILTRMAX,NDATAMAX),
     *          sti(NFILTRMAX,NDATAMAX), dreelr(NDATAMAX),
     *          dreeli(NDATAMAX), dsyntr(NDATAMAX),dsynti(NDATAMAX)
      character entr*30, titre*80
c
c Programme extrayant les donnees reelles et les erreurs-donnees
c du fichier-resultat de cross8 en vue de leur integration dans
c le fichier-donnee de l'inversion (tv.0)
c
      write(0,*)' Programme extrayant les donnees reelles et'
      write(0,*)' les erreurs-donnees du fichier-resultat'
      write(0,*)' de cross8 en vue de leur integration dans'
      write(0,*)' le fichier-donnee de l''inversion (tv.0)'
      igk=  11
c
    2 write(0,'("fichier contenant les gk^ ? ",$)')
      read(*,'(a)')entr
      open(igk,status='old',err=2,file=entr)
      call lecgk(igk,km,nfiltr,np,ip,nf1min,nf2max,dnu,
     *                                  dt,dist,ienv,ierr)
      write(*,*) km,nfiltr,ienv,'    nb modes, nb periodes, ienveloppe'
c     ienv = 1 si on a travaille sur l'enveloppe du cross-correl.
c              sauf si on a picke une seule donnee
c          = 0 si on a travaille sur le cross-correlogramme lui-meme.
      if(ierr.eq.1)then
                close(igk)
                stop 'Erreur dans le sous-programme "lecgk"'
                endif
c
      do 200 ks=1,km
c
c     syntr et synti = parties reelle et imaginaire de la TF du
c                    signal synthetique pour les differents modes
  200 read(igk,*)(syntr(i,ks),synti(i,ks),i=nf1min,nf2max)
c
c
c
      do 10 k=1,km
c
c
      read(igk,*)(refr(i),refi(i),i=nf1min,nf2max)
c
c     refr et refi = parties reelle et imaginaire de la TF du
c                    signal de reference pour le mode k
c
c
      do 101 l=1,nfiltr
c
                    read(igk,*) npp,nf1(l),nband(l),pper(l),ppar(l)
                    read(igk,*)ndata(l),(idata(l,i),i=1,ndata(l))
                    read(igk,*)ndata(l),(str(l,i),sti(l,i),i=1,ndata(l))
  101 continue
c
   10 continue
c
      do 40 k=1,km
      do 40 l=1,nfiltr
      read(igk,*)ndat,(dreelr(i),dreeli(i),i=1,ndat)
      if(ndat.eq.1.and.dreelr(1).eq.0..and.dreeli(1).eq.0.) then
      write(*,*)'  1     0.   0.    pas de donnee'
      else
      write(*,'(i3,(10e14.6))')ndat,(dreelr(i),dreeli(i),i=1,ndat)
      endif
   40 continue
      do 42 k=1,km
      do 42 l=1,nfiltr
      read(igk,*)ndat,(dsyntr(i),dsynti(i),i=1,ndat)
   42 continue
      do 43 i=1,ndat
      if(dsyntr(i).ne.str(nfiltr,i).or.dsynti(i).ne.sti(nfiltr,i))
     *    stop 'incoherence des donnees synthetiques'
   43 continue
      read(igk,'(a)') titre
      do 44 k=1,km
      do 44 l=1,nfiltr
      read(igk,*)ndat,(dsyntr(i),i=1,ndat)
      write(*,'(i3,(10e14.6))')ndat,(dsyntr(i),i=1,ndat)
   44 continue
      close(igk)
      end
c ----------------------------------------------------------------------
      subroutine lecgk(igk,km,nfiltr,np,ip,nf1min,nf2max,dnu,dt,
     *                 dist,ienv,ierr)
      character*80 bcd(4)
      character*30 entr
      character*30 entref
      read(igk,'(a)')(bcd(i),i=1,4)
      read(igk,'(a)')entr
      read(igk,'(a)')entref
      ierr=0
      write(0,'(a)') ' sismogramme reel : ', (bcd(i),i=1,4)
      write(0,'(2a)') '      synthetique = ', entr
      write(0,'(2a)') '      reference   = ', entref
1     read(igk,*)km,nfiltr,np,ip,nf1min,nf2max,dnu,dt,dist,ienv
      return
      end
