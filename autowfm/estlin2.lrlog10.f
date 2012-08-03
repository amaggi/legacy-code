      parameter(NPERMAX=30,NPARMAX=7,NCOUCHMAX=35,NMODMAX=10,mp=200)
      dimension dcdp(NPERMAX,NCOUCHMAX,NPARMAX,NMODMAX)
      dimension nper(NMODMAX),t(NPERMAX,NMODMAX), c(NPERMAX,NMODMAX)
      dimension alpha(NPERMAX,NMODMAX)
      dimension P0(mp),P(mp)
      integer seirei,seilov
      character sortr1*30,sortr2*30,sortl1*30,sortl2*30
      character sort1*30,sort2*30,entr*30
      character entrr*30,entrl*30,titre*80
c
c     Programme d'estimation lineaire des vitesses de phase
c     et de l'attenuation (coefficient alpha) pour un modele perturbe.
c     Les derivees partielles de alpha par rapport a Q**-1 doivent
c     correspondre au dernier parametre de la matrice dcdp.
c     Il n'est pas prevu d'inverser pour la seule attenuation: on doit
c     obligatoirement avoir au moins un paquet de der. part. de la
c     vitesse de phase avant les der. part. de l'attenuation.
c     On suppose (cf ss-programme de lecture) que le starting model
c     de l'inversion est identique a l'a priori model.
c---------------------------------------------------------------
c       Ce programme a ete adapte a l'inversion de plusieurs
c       trajets sur les deux composantes love et Rayleigh(2/06/93).
c       Derniere modif le 22/09/94
c       Il a aussi ete modifie le 6 mai 1997
c---------------------------------------------------------------
c
      pi=3.14159
      idcdpr=10
      idcdpl=14
      imodr=11
      imodl=15
      iattr=12
      iattl=16
      iexray=0
      iexlov=0
c
      write(*,*) 'entrez le nombre de seismes:'
      read(*,*) nbs
      write(*,*) 'nombre de seismes inverses pour rayleigh et love:'
      read(*,*) nlr
      nbt=nbs+nlr 
      write(*,*) ' le nombre de trajets total est:',nbt
      write(*,*) 'nombre de seismes inverses pour rayleigh:'
      read(*,*) seirei
      seilov=nbt-seirei
      write(*,*) ' Le nombre de seismes inverses pour love est:',seilov
 
c   Fichiers en entree
      write(0,'(a,$)') 'fichier contenant les dc/dp Rayleigh? :'
      read(*,*)entrr
      write(0,'(a,$)') 'fichier contenant les dc/dp love? :'
      read(*,*)entrl
      write(0,'("fichier de sortie dispersion Rayleigh
     *(disp.ray.itn) ? :",$)')
      read(*,'(a)')sortr1 
      write(0,'("fichier de sortie dispersion love
     *(disp.lov.itn) ?:",$)')
      read(*,'(a)')sortl1 
      write(0,'("fichier de sortie attenuation Rayleigh
     *(att.ray.itn) ? : ",$)')
      read(*,'(a)')sortr2
      write(0,'("fichier de sortie attenuation love
     *(att.lov.itn) ? : ",$)')
      read(*,'(a)')sortl2

      if (seirei.ne.0) iexray=1
      if (seilov.ne.0) iexlov=1
     
      call lecmod(P0,P,mp,titre,nvar,ncouch,nmode,nbs)
      write(*,*) nmode,' mode(s)'

      do 20 ilovray=1,iexray+iexlov


c   Le fichier dc/dp sera t-il rayleigh ou love?
      if (ilovray.eq.1.and.iexray.eq.1) then
      idcdp=idcdpr
      imod=imodr  
      iatt=iattr
      entr=entrr
      sort1=sortr1
      sort2=sortr2
      else
      idcdp=idcdpl
      imod=imodl  
      iatt=iattl
      entr=entrl
      sort1=sortl1
      sort2=sortl2
      endif
      
c  Le choix est fait, il ne reste plus qu'a ouvrir les fichiers.
 
      open(idcdp,status='old',file=entr)
      read(idcdp,'(a)') titre
      read(idcdp,*) nmode
      write(*,*) titre
      open(imod,status='new',file=sort1)
      open(iatt,status='new',file=sort2)


      read(idcdp,*)  npar, nprof, iq
      if(npar.ne.nvar. or .nprof.ne.ncouch) stop 'Incompatibilite'

      call lecdp(dcdp,t,c,alpha,npar,nprof,nper,nmode,idcdp,ierr)
      if(ierr.eq.1)then
                close(idcdp)
                stop 'Erreur dans le sous-programme "lecdp"'
                endif
c
      write(imod,'(a)') titre
      write(iatt,'(a)') titre
      write(imod,*) nmode
      m=nvar*ncouch
c --on ecrit le moment sismique et non le log10 du Moment sismique
      write(iatt,*) nmode,(10**P0(isei), 10**P(isei),isei=m+1,m+nbs)
     *,'    Moment'
      do 10 k=1,nmode
      write(imod,*) nper(k)
      write(iatt,*) nper(k)
      do 10 n=1,nper(k)
      cest=0.
      attest=0.
      l=0
      if (iq.eq.1) then

      do 11 i=1,ncouch
      do 12 j=1,nvar-1
      l=l+1
      cest=cest+(P(l)-P0(l))*dcdp(n,i,j,k)
12    continue
      l=l+1
c  Attention: le dernier paquet de derivees partielles
c      doit etre relatif a d(attenuation)/dQ**-1
      attest=attest+dcdp(n,i,nvar,k)*(P(l)-P0(l))
11    continue

      else

      do 111 i=1,ncouch
      do 111 j=1,nvar
      l=l+1
      cest=cest+(P(l)-P0(l))*dcdp(n,i,j,k)
111   continue

      endif

      cest=c(n,k)+cest
      attest=alpha(n,k)+ attest*2*pi/t(n,k)/(c(n,k)*c(n,k))
c     on a lu dans dcdp des pseudo der.part.,
c     ie 0.5*C*(1.33*beta2/alpha2*alpha/C*dC/dalpha+ beta/C*dC/dbeta)
c        au lieu de pi/C/T*(1.33 etc..) qui est dalpha/dQbeta**-1
      write(imod,*) t(n,k), cest,  '0.00'
      write(iatt,*) t(n,k), attest
10    continue

      close(imod)
      close(iatt)
      close(idcdp)
20    continue
      end
c ----------------------------------------------------------------------
      subroutine lecdp(dp,t,c,alpha,npar,ncouch,nper,kmode,lu,ierr)
      parameter(NPERMAX=30,NPARMAX=7,NCOUCHMAX=35,NPMAX=1024,NPMAX2=513)
      parameter(NMODMAX=10)
*****     LECTURE DES DERIVEES PARTIELLES (SAITO PAR EX.)    *****
*     NPERMAX=      nombre maxi de periodes auxquelles existent
*                     les derivees partielles pour un mode isole
*     NPARMAX=     nombre maxi de parametres du modele de
*                     Terre (5 elastiques + densite = 6)
*     NCOUCHMAX=   nombre de couches du modele de Terre
      dimension nper(NMODMAX),t(NPERMAX,NMODMAX), c(NPERMAX,NMODMAX)
      dimension alpha(NPERMAX,NMODMAX)
      dimension dp(NPERMAX, NCOUCHMAX, NPARMAX,NMODMAX)
c     dimension ep(NCOUCHMAX),amodel(NPARMAX,NCOUCHMAX)
      do 10 k=1,kmode
      write(0,*) ' Lecture der. part. du mode',k,'  sur',kmode
      read(lu,*) nper(k)
*  LECTURE DES DERIVEES PARTIELLES
      do 1 i=1,nper(k)
      read(lu,*,end=9,err=99) t(i,k), c(i,k), bid, alpha(i,k)
      do 1 ipar=1,npar
c changed eric for linux (2/4/2002)
      read(lu,'(1x,5e16.8)',end=9,err=99)(dp(i,icouch,ipar,k),
     *                                 icouch=1,ncouch)
c     read(lu,'(5e16.8)',end=9,err=99)(dp(i,icouch,ipar,k),
c    *                                 icouch=1,ncouch)
    1 continue
   10 continue
*  PASSAGE AUX DERIVEES PARTIELLES ABSOLUES A PARTIR DES RELATIVES
c     read(lu,*)(ep(icouch),i=1,ncouch)
c     read(lu,*)((amodel(ipar,icouch),icouch=1,ncouch),ipar=1,npar)
c     do 2 k=1,kmode
c     do 2 i=1,nper(k)
c     do 2 ipar=1,npar
c     do 2 icouch=1,ncouch
c     dp(i,icouch,ipar,k)=c(i,k)/amodel(ipar,icouch)*dp(i,icouch,ipar,k)
c   2 continue
      return
    9 write(*,*)'Fin de fichier detectee dans lecture der. part., mode '
     *           ,kmode
      ierr=1
      return
   99 write(*,*)'Erreur detectee dans lecture der. part., mode ', kmode
      ierr=1
      return
      end
c----------------------------------------------------------------------
      Subroutine lecmod(P0,Pinit,mp,titre,nvar,ncouch,nmode,nbs)
c     modifiee le 13/05/93 afin de lire des fichiers tv.n contenant
c     plusieurs moments sismiques(cas de l'inversion de plusieurs trajets).

      Parameter (mmp=200)
      Dimension P0(mp), Pinit(mp), prof(mmp),DP0(mmp)
      character name*40, titre*80

      if(mmp.ne.mp) stop 'erreur de dimension dans lecmod'

c   Attention: ne pas confondre nombre de seismes et nbre de trajets.
c   voir commentaires procedures init.lr.run
 
      il1=9
      write(*,*) 'entrez le nom du fichier-modele (tv.n)'
      read(*,*) name
      open(il1, file=name, status='old', err=999)

      read(il1,'(a)') titre
      read(il1,*) nvar, ncouch
      m=nvar*ncouch

c   lecture du modele a priori
      do 1 i=1,ncouch
1     read(il1,*) prof(i),(P0(k),k=(i-1)*nvar+1,i*nvar)
      do 11 k=m+1,m+nbs
11    read(il1,*) P0(k)

c   lecture des ecarts-type a priori pour ce modele
      do 2 i=1,ncouch
2     read(il1,*) prof(i),(DP0(k),k=(i-1)*nvar+1,i*nvar)
      do 22 k=m+1,m+nbs
22    read(il1,*) DP0(k)

c   lecture du modele resultant de l'iteration precedente
      do 3 i=1,ncouch
3     read(il1,*) prof(i),(Pinit(k),k=(i-1)*nvar+1,i*nvar)
      do 33 k=m+1,m+nbs
33    read(il1,*) Pinit(k)

      read(il1,*) distcor
      do 4 i=1,nvar
    4 read(il1,*) coupl
      read(il1,*) nmod2
      if(nmod2.lt.nmode) nmode=nmod2

      close(il1)

      return
  999 stop 'erreur de lecture fichier'
      end
