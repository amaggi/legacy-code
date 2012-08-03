c  ce pgm est a compiler avec les instructions suivantes en se logant
c  sur queyras (28 avril 1993)
c  f77 -misalign -g -c convdp2.lr.f convdp2.lr.o
c  f77 -misalign -g convdp2.lr.o des.cross.o suncore.ps.o -o convdp2.lr 
c  -lcore77 -lcore -lsunwindow -lpixrect -lm


c     Version permettant l'inversion du moment sismique
c------------------------------------------------------------------------
c     modifiee le 28/04/93 pour permettre l'inversion simultanee
c     de plusieurs trajets (meme seisme enregistre sur sa composante love
c     et Rayleigh ou plusieurs seismes ayant des trajets proches)
c     les modifs sont encadrees par des lignes de pointilles.
c     modif le 7 May 1997 pour permettre l'inversion de log10Q
c     on ajoute un facteur (Q**-1)*ln(10) aux pseudo derivees partielles 
c     de l'attenuation
c ----------------------------------------------------------------------- 

      parameter(NPERMAX=30,NPARMAX=7,NCOUCHMAX=100,NMODMAX=10)
      implicit real*8 (a-h,o-z)
      real*4 t,c,u,alpha
      dimension am0(10),dam0(10)
      dimension dcdp(NPERMAX,NCOUCHMAX,NPARMAX,NMODMAX)
      dimension nper(NMODMAX),t(NPERMAX,NMODMAX), c(NPERMAX,NMODMAX)
      dimension u(NPERMAX,NMODMAX), alpha(NPERMAX,NMODMAX)
      dimension a(NCOUCHMAX), indic(NPARMAX)
      dimension ncouch1(NCOUCHMAX), ncouch2(NCOUCHMAX)
      dimension ep(NCOUCHMAX), xx(NCOUCHMAX)
      dimension amod(NPARMAX,NCOUCHMAX), amod2(NPARMAX,NCOUCHMAX)
      dimension ectyp(NPARMAX),coupl(NPARMAX,NPARMAX)
      character entr*80, sort1*80, titre*80, carte*80, form*80
c
      idcdp=10
      itv=11
      io=   12
c
    1 write(0,'(a,$)') 'fichier contenant les dc/dp ? :'
      read(*,'(a)')entr
	write(*,*)entr
      open(idcdp,status='old',err=1,file=entr)
      read(idcdp,'(a)') titre
      read(idcdp,'(a)') carte
      read(carte,*) nmode
      write(*,*) titre(:lnblnk(titre))
      write(*,*) nmode,' mode(s)'
    2 write(0,'(a,$)') 'fichier de sortie du modele ? :'
      read(*,*)entr
      write(*,*)entr
c---ce fichier ne peut etre declare 'new'
      open(itv,err=2,file=entr)
      write(itv,'(a)') titre(:lnblnk(titre))
    3 write(0,'("fichier de sortie des der. part. abs. ? : ",$)')
      read(*,'(a)')sort1
      write(*,'(a)')sort1
      open(io,status='new',err=3,file=sort1)
c
      read(idcdp,*)  npar, nprof
      write(0,*) 'Modele a',npar,' parametres, ',nprof,' couches'
      if(nprof.gt.NCOUCHMAX) stop 'redimensionner le tableau NCOUCHMAX
     *dans convdp2.lr.f'
c
      write(0,*)'Entrer 0 ou 1 (pour exclure ou conserver le parametre)'
      write(0,*)' ncouch1,  ncouch2 (indices des couches extremes a cons
     *erver)'
      nvar=0
      ncmin=nprof
      ncmax=1
      do 10 i=1,npar
      write(0,'("Parametre  ",i3," : ",$)') i
      read(*,*) indic(i),ncouch1(i),ncouch2(i)
      nvar=nvar+indic(i)
      if(indic(i).eq.0) go to 10
      if(ncouch1(i).lt.ncmin) ncmin=ncouch1(i)
      if(ncouch2(i).gt.ncmax) ncmax=ncouch2(i)
10    continue
      ncouch=ncmax-ncmin+1
      write(0,*)'Entrer 1 si le dernier parametre conserve est le',
     *          ' facteur de qualite'
      read(*,*) iq
      
      write(0,*) 'Debut de lecture des der. part.'
c
      call lecdp(dcdp,t,c,u,alpha,npar,nprof,nper,nmode,idcdp,ierr)

      if(ierr.eq.1)then
                close(idcdp)
                stop 'Erreur dans le sous-programme "lecdp"'
                endif


      read(idcdp,*) xx(1), xx(2), (ep(j),j=2,nprof)
      call prof(ep,xx,ncmin,ncmax,ncmax)
      do 101 i=1,npar
      read(idcdp,*) (amod(i,j),j=1,nprof)
  101 continue
      l=0
      do 102 i=1,npar
      if(indic(i).eq.0) goto 102
      l=l+1
      do 103 j=ncmin,ncmax
      amod2(l,j)=amod(i,j)
  103 continue
  102 continue


c modif eric--------------------------------------------------
c On ecrit les derivees partielles absolues dans le fichier dpa.
c Pour Q on multiplie les pseudo dp par un facteur ln(10)*1/Q 
c (inversion de log10(Q))

      write(io,'(a)') titre(:lnblnk(titre))
      write(io,'(a)') carte(:lnblnk(carte))
      write(io,*)  nvar, ncouch, iq
      do 11 k=1,nmode
      write(io,*)  nper(k)
      do 12 n=1,nper(k)
c       write(io,*) t(n,k), c(n,k), u(n,k), alpha(n,k)
c changed by eric for linux  (2/04/2002)
        write(io,1002) t(n,k), c(n,k), u(n,k), alpha(n,k)
1002  format(f5.0,2(1x,f11.8),1x,e15.8)
        do 13 i=1,npar
          if(indic(i).eq.0) goto 13
          l=0
          do 14 j=ncmin,ncmax
             if (i.eq.npar) then 
                dcdp(n,j,i,k)=-(dcdp(n,j,i,k)*amod(i,j))*alog(10.)
             endif
             l=l+1
             a(l)=dcdp(n,j,i,k)
             if(j.lt.ncouch1(i).or.j.gt.ncouch2(i)) a(l)=0.
14        continue
c changed by eric for linux  (2/04/2002)
c         write(io,'(5d16.8)') (a(ii),ii=1,l)
          write(io,1003)(a(ii),ii=1,l)
1003  format('D',5d16.8)
13      continue
12    continue
11    continue
c--------------------------------------------------------------
c--pour l'ecriture du modele a priori dans tv.init on ecrit
c--maintenant log10 Q au lieu de 1/Q.

      write(itv,*) nvar, ncmax-ncmin+1
      do 104 j=ncmin,ncmax
      amod2(nvar,j)=log10(1./amod2(nvar,j))
      write(itv,'(f10.4,7f10.6)') xx(j),(amod2(i,j),i=1,nvar)
  104 continue
      write(*,*) ' Entrer  1 pour des ecarts-types absolus,'
      write(*,*) ' Entrer -1 pour des ecarts-types relatifs'
      read(*,*) iec
      if(iec.eq.1) then
      write(*,*) ' Entrer les ecarts-type (indep. de la prof.):'
      write(*,*) ' (On inverse log10Q donc l\'ecart-type '
      write(*,*) ' a entrer est un sigma sur log10Q ) '
      read(*,*) (ectyp(i),i=1,nvar)
      else
      write(*,*) ' Entrer les facteurs pour calculer les ec-type:'
      read(*,*) (ectyp(i),i=1,nvar)
      endif
c---------------------------------------------------------------------
      write(*,*)' Entrer le nombre de seismes total que vous'
      write(*,*)' voulez inverser:'
      read (*,*) nbsei 

      write(*,*)' Entrer le moment sismique et son incertitude'
      write(*,*)'                        (relative ou absolue)'
      read(*,*) (am0(isei), dam0(isei),isei=1,nbsei)
      do 120 isei=1,nbsei
      if(dam0(isei).lt.0.) dam0(isei)=am0(isei)*abs(dam0(isei))
      write(itv,*) am0(isei), '            Moment sismique'
120   continue

c---------------------------------------------------------------------
	if(iec.eq.1) then
      do 105 j=ncmin,ncmax
      write(itv,'(f10.4,7f10.6)') xx(j),(ectyp(i),i=1,nvar)
  105 continue
      	else
      do 106 j=ncmin,ncmax
      write(itv,'(f10.4,7f10.6)') xx(j),(amod2(i,j)*ectyp(i),i=1,nvar)
  106 continue
	endif
c---------------------------------------------------------------------
      do 130 isei=1,nbsei
      write(itv,*) dam0(isei), '            erreur sur moment sismique'
130   continue
c---------------------------------------------------------------------
c     Ecriture du starting model

      do 107 j=ncmin,ncmax
      write(itv,'(f10.4,7f10.6)') xx(j),(amod2(i,j),i=1,nvar)
  107 continue
c---------------------------------------------------------------------
      do 140 isei=1,nbsei
      write(itv,*) am0(isei), '            moment sismique'
140   continue
c--------------------------------------------------------------------- 
      write(*,*)' Entrer la distance de correlation'
      read(*,*) distcorr
      write(itv,*) distcorr, '            distance de correlation'
      write(*,*)' Entrer les couplages (i,j) entre parametres'
	write(form,'(a,i3,a)') '(',nvar,'f8.4,a,i3)'
      do 108 i=1,nvar
      read(*,*) (coupl(i,j),j=1,nvar)
      write(itv,form) (coupl(i,j),j=1,nvar),'         couplages',i
  108 continue
      
      end
c ----------------------------------------------------------------------
      subroutine prof(ep,xx,ncouc1,ncouc2,m)
      implicit real*8 (a-h,o-z)
*     passage des epaisseurs aux profondeurs
      parameter (NCOUCHMAX=100)

      dimension ep(NCOUCHMAX),xx(NCOUCHMAX)
***      xx(1)=ep(1)
***      xx(2)=ep(1)+ep(2)/2.
***      pro=ep(1)+ep(2)
***	pro=xx(2)+ep(2)/2.
***      do 1 i=3,ncouc2
***      xx(i)=2.*pro-xx(i-1)
***    1 pro=pro+ep(i)
***      xx(1)=2.*xx(1)-xx(2)

c	xx(1)=     lue dans le fichier dp
c	xx(2)=hmin lue dans le fichier dp
	do 1 i=3,ncouc2
	xx(i)=xx(i-2)+2*ep(i-1)
    1	continue

      do 2 i=ncouc2+1,m+1
      j=i-(ncouc2-ncouc1+1)
    2 xx(i)=xx(j)
 
      return

      end
c ----------------------------------------------------------------------
      subroutine lecdp(dp,t,c,u,alpha,npar,ncouch,nper,kmode,lu,ierr)
      parameter(NPERMAX=30,NPARMAX=7,NCOUCHMAX=100,NPMAX=1024)
      parameter(NPMAX2=513,NMODMAX=10)
      implicit real*8 (a-h,o-z)
      real*4 t,c,u,alpha
*****     LECTURE DES DERIVEES PARTIELLES (SAITO PAR EX.)    *****
*     NPERMAX=      nombre maxi de periodes auxquelles existent
*                     les derivees partielles pour un mode isole
*     NPARMAX=     nombre maxi de parametres du modele de
*                     Terre (5 elastiques + densite = 6)
*                                         + Qbeta   = 7)
*     NCOUCHMAX=   nombre de couches du modele de Terre
      dimension nper(NMODMAX),t(NPERMAX,NMODMAX), c(NPERMAX,NMODMAX)
      dimension u(NPERMAX,NMODMAX), alpha(NPERMAX,NMODMAX)
      dimension dp(NPERMAX, NCOUCHMAX, NPARMAX,NMODMAX)
c     dimension ep(NCOUCHMAX),amodel(NPARMAX,NCOUCHMAX)
      do 10 k=1,kmode
      write(0,*) ' Lecture der. part. du mode',k,'  sur',kmode
      read(lu,*) nper(k)
*  LECTURE DES DERIVEES PARTIELLES
      do 1 i=1,nper(k)
      read(lu,*,end=9,err=99) t(i,k), c(i,k), u(i,k), alpha(i,k)
      do 1 ipar=1,npar
      read(lu,'(5d16.8)',end=9,err=99)(dp(i,icouch,ipar,k),
     *                                 icouch=1,ncouch)
c
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
     *           ,k
      write(*,*)'Periode', t(i,k), ' Parametre', ipar,' Couche',icouch
      ierr=1
      return
c
   99 write(*,*)'Erreur detectee dans lecture der. part., mode ', k
      write(*,*)'Periode', t(i,k), ' Parametre', ipar,' Couche',icouch
      ierr=1
      return
      end
