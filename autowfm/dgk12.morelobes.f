      parameter(NPERMAX=30,NPARMAX=7,NCOUCHMAX=35,NPMAX=1024,NPMAX2=513)
c     parameter(NFILTRMAX=5,NDATAMAX=5,NMODMAX=10)
c     modif eric Oct 97 inversion +s lobes: NDATAMAX=10
c     parameter(NFILTRMAX=5,NDATAMAX=10,NMODMAX=10)

c     modif eric 10 Oct 2000 : NFILTRMAX=10
      parameter(NFILTRMAX=10,NDATAMAX=10,NMODMAX=10)

      dimension dcdp(NPERMAX,NCOUCHMAX,NPARMAX,NMODMAX),
     *        work(NPMAX),w(NPMAX2),fgauss(NFILTRMAX,NPMAX2)
      dimension nper(NMODMAX),t(NPERMAX,NMODMAX),c(NPERMAX,NMODMAX)
      dimension pper(NFILTRMAX),ppar(NFILTRMAX),ndata(NFILTRMAX),
     *          idata(NFILTRMAX,NDATAMAX),nf1(NFILTRMAX),nf2(NFILTRMAX),
     *          nband(NFILTRMAX), refr(NPMAX),refi(NPMAX),
     *          syntr(NPMAX,NMODMAX), synti(NPMAX,NMODMAX)
      dimension gr(NPMAX),gi(NPMAX),
     *          dgr(NPMAX),dgi(NPMAX), str(NFILTRMAX,NDATAMAX),
     *          sti(NFILTRMAX,NDATAMAX), dreelr(NDATAMAX),
     *          dreeli(NDATAMAX), dsyntr(NDATAMAX),dsynti(NDATAMAX)
      character entr*30, sort1*30, titre*80
      real*8 pi, pi2
c
c Programme de calcul des derivees partielles des intercorrelogrammes
c entre signaux multimodes et modes purs
c
      idcdp=10
      igk=  11
      io=   12
c
      pi=3.1415926536
      pi2=2*pi
    1 write(0,'(a,$)') 'fichier contenant les dc/dp ? :'
      read(*,*)entr
      open(idcdp,status='old',err=1,file=entr)
      read(idcdp,'(a)') titre
      read(idcdp,*) nmode
      write(*,*) titre
      write(*,*) nmode,' mode(s)'
    2 write(0,'("fichier contenant les gk^ ? ",$)')
      read(*,'(a)')entr
      open(igk,status='old',err=2,file=entr)
    3 write(0,'("fichier de sortie des dgk/dp ? : ",$)')
      read(*,'(a)')sort1
      open(io,status='new',err=3,file=sort1)
      call lecgk(igk,km,nfiltr,np,ip,nf1min,nf2max,dnu,
     *                                  dt,dist,ienv,ierr)
c     ienv = 1 si on a travaille sur l'enveloppe du cross-correl.
c              sauf si on a picke une seule donnee
c          = 0 si on a travaille sur le cross-correlogramme lui-meme.
      if(ierr.eq.1)then
                close(igk)
                stop 'Erreur dans le sous-programme "lecgk"'
                endif
c
      if(nmode.lt.km) stop 'der. part. insuffisantes '
      read(idcdp,*)  npar, nprof, iq
      write(io,*) entr,km,nfiltr,npar,nprof,ienv
      write(*,*)
     *     'distance,enveloppe,nb-modes,nb-per-filtr.,npar,nprof:'
      write(*,*)dist,ienv,km,nfiltr,npar,nprof
c
C     write(0,*)'entrer le nb max de parametres (Vs,ro,Vp,XI,PHI,ETA)'
C     read(*,*) nparlim
C     write(0,*)'entrer le nb max de modes de reference pour les d.p.'
C     read(*,*) klim
c
      call lecdp(dcdp,t,c,npar,nprof,nper,km,idcdp,ierr)
      if(ierr.eq.1)then
                close(idcdp)
                stop 'Erreur dans le sous-programme "lecdp"'
                endif
c
      do 200 ks=1,km
c                                                                       ks
c     syntr et synti = parties reelle et imaginaire de la TF du         ks
c                    signal synthetique pour les differents modes       ks
  200 read(igk,*)(syntr(i,ks),synti(i,ks),i=nf1min,nf2max)              ks
c
c
      do 86 ks=1,km
                freq1=(nf1min-1)*dnu                                    ks
                tnu1=1/freq1                                            ks
                freq2=(nf2max-1)*dnu                                    ks
                tnu2=1/freq2                                            ks
                tmin=min(t(1,ks),t(nper(ks),ks))                        ks
                tmax=max(t(1,ks),t(nper(ks),ks))                        ks
                write(*,*)'Mode synt.',ks                               ks
                write(*,*)'Tdpmin= ',tmin,' Tdpmax= ',tmax              ks
                write(*,*)'Tharmin=',tnu2,' Tharmax=',tnu1              ks
	if(nper(ks).lt.5)write(*,*)'Attention, il n''y a que',nper(ks),
     *	'periodes pour le mode',ks,'<spline> en exige 5.'
ccc   80    if(tmin.le.1/freq2.or.freq2.lt.freq1) goto 81                  ks
ccc         write(*,*) 'Mode',ks,' : interpolation impossible a T=',1/freq2  ks
ccc         freq2=freq2-dnu                                                ks
ccc         goto 80                                                        ks
ccc   81    if(tmax.ge.1/freq1.or.freq1.eq.freq2) goto 86                  ks
ccc         write(*,*) 'Mode',ks,' : interpolation impossible a T=',1/freq1  ks
ccc         freq1=freq1+dnu                                                ks
ccc         goto 81                                                        ks
   86 continue                                                          ks
c
      do 10 k=1,km
      inew=1                                                            k
c                                                                       k
c                                                                       k
      read(igk,*)(refr(i),refi(i),i=nf1min,nf2max)                      k
c                                                                       k
c     refr et refi = parties reelle et imaginaire de la TF du           k
c                    signal de reference pour le mode k                 k
c                                                                       k
      write(0,*)' CALCUL POUR LE MODE DE REFERENCE',k,' sur',km         k
c                                                                       k
      do 10 ipar=1,npar                                                 k
      do 10 iprof=1,nprof                                               k ipar
c                                                                       k ipar iprof
      do 87 i=1,npp                                                     k ipar iprof
      gr(i)=0.                                                          k ipar iprof i
      gi(i)=0.                                                          k ipar iprof i
   87 continue                                                          k ipar iprof i
c                                                                       k ipar iprof
      do 88 ks=1,km                                                     k ipar iprof
c     Boucle sur les modes synthetiques                                 k ipar iprof ks
c                                                                       k ipar iprof ks
      call spline(ks,nper,ipar,iprof,t,c,dcdp,nf1min,nf2max,dnu,w)      k ipar iprof ks
      do 8 i=nf1min,nf2max                                              k ipar iprof ks
      f= pi2*(i-1)*dnu*w(i)*dist                                        k ipar iprof ks i
      gr(i)=gr(i) - f*(syntr(i,ks)*refi(i)-synti(i,ks)*refr(i))         k ipar iprof ks i
      gi(i)=gi(i) - f*(syntr(i,ks)*refr(i)+synti(i,ks)*refi(i))         k ipar iprof ks i
    8 continue                                                          k ipar iprof ks i
   88 continue                                                          k ipar iprof ks
      if (iq.eq.1.and.ipar.eq.npar) then                                k ipar iprof
c       on traite le cas  dgk/dQbeta**-1                                k ipar iprof
c          on tient maintenant compte du i (imaginaire pur)             k ipar iprof
c          qui manque dans la pseudo derivees partielle lue dans lecdp  k ipar iprof
                do 9 i=nf1min,nf2max                                    k ipar iprof
                gsave=gr(i)                                             k ipar iprof i
                gr(i)=-gi(i)                                            k ipar iprof i
                gi(i)=gsave                                             k ipar iprof i
    9           continue                                                k ipar iprof i
                endif                                                   k ipar iprof
c                                                                       k ipar iprof
      do 101 l=1,nfiltr                                                 k ipar iprof
c                                                                       k ipar iprof l
c     pper et ppar = periode centrale et largeur du filtre              k ipar iprof l
c     idata = indices temporels des points auxquels on                  k ipar iprof l
c             inverse (et donc auxquels on calcule les d.p.)            k ipar iprof l
c     str = partie relle du cross-corr. temporel                        k ipar iprof l
c     sti = partie imaginaire du cross-correlogramme temporel           k ipar iprof l
c                                                                       k ipar iprof l
      if(inew.eq.1) then                                                k ipar iprof l
                    read(igk,*) npp,nf1(l),nband(l),pper(l),ppar(l)     k ipar iprof l
                    read(igk,*)ndata(l),(idata(l,i),i=1,ndata(l))       k ipar iprof l
                    read(igk,*)ndata(l),(str(l,i),sti(l,i),i=1,ndata(l))  k ipar iprof l
                    nf2(l)=nf1(l)+nband(l)-1                            k ipar iprof l
                    call szgauss(l,npp,nf1(l),nband(l),fgauss)          k ipar iprof l

c		der. part. du moment sismique
		do 30 i=1,ndata(l)
		dgr(idata(l,i))=str(l,i)
		dgi(idata(l,i))=sti(l,i)
   30		continue
          if(idata(l,1).eq.0) then                                      k ipar iprof l
c cas de la non-validation du pickage automatique sur le cross-corr.    k ipar iprof l
c      la derivee partielle est forcee a 0.                             k ipar iprof l
                    dgr(1)=0.                                           k ipar iprof l
                    dgi(1)=0.                                           k ipar iprof l
      write(io,*) k,pper(l),ppar(l),ipar,iprof,ndata(l),                k ipar iprof l
     *           (idata(l,i),i=1,ndata(l)), dgr(1),dgi(1)               k ipar iprof l
           else                                                         k ipar iprof l
      if(ienv.eq.1.and.ndata(l).ne.1)                                   k ipar iprof l
     *      call dpenv(l,dgr,dgi,ndata(l),idata,str,sti)                k ipar iprof l
      if(ienv.eq.1.and.ndata(l).eq.1)                                   k ipar iprof l
     *      call dpphase(l,dgr,dgi,ndata(l),idata,str,sti)              k ipar iprof l
      write(io,*) k,pper(l),ppar(l),ipar,iprof,                         k ipar iprof l
     *            ndata(l),(idata(l,i),i=1,ndata(l))                    k ipar iprof l
      write(io,'(20e14.6," M0")')                                       k ipar iprof l
     *            (dgr(idata(l,i)),dgi(idata(l,i)),i=1,ndata(l))        k ipar iprof l

                               endif                                    k ipar iprof l

                    endif                                               k ipar iprof l
      do 89 i=1,npp                                                     k ipar iprof l
      dgr(i)=gr(i)*fgauss(l,i)                                          k ipar iprof l i
      dgi(i)=gi(i)*fgauss(l,i)                                          k ipar iprof l i
   89 continue                                                          k ipar iprof l i
      do 90 i=npp+1,np                                                  k ipar iprof l
      dgr(i)=0.                                                         k ipar iprof l i
      dgi(i)=0.                                                         k ipar iprof l i
   90 continue                                                          k ipar iprof l i
      call nlogn(ip,dgr,dgi,+1.)                                        k ipar iprof l
c                                                                       k ipar iprof l
c  dgr et dgi sont maintenant des fonctions temporelles                 k ipar iprof l
c                                                                       k ipar iprof l
      call sztrans(work,dgr,np)                                         k ipar iprof l
      call sztrans(work,dgi,np)                                         k ipar iprof l
          if(idata(l,1).eq.0) then                                      k ipar iprof l
c cas de la non-validation du pickage automatique sur le cross-corr.    k ipar iprof l
c      la derivee partielle est forcee a 0.                             k ipar iprof l
                    dgr(1)=0.                                           k ipar iprof l
                    dgi(1)=0.                                           k ipar iprof l
      write(io,*) k,pper(l),ppar(l),ipar,iprof,ndata(l),                k ipar iprof l
     *           (idata(l,i),i=1,ndata(l)), dgr(1),dgi(1)               k ipar iprof l
           else                                                         k ipar iprof l
      if(ienv.eq.1.and.ndata(l).ne.1)                                   k ipar iprof l
     *      call dpenv(l,dgr,dgi,ndata(l),idata,str,sti)                k ipar iprof l
      if(ienv.eq.1.and.ndata(l).eq.1)                                   k ipar iprof l
     *      call dpphase(l,dgr,dgi,ndata(l),idata,str,sti)              k ipar iprof l
      write(io,*) k,pper(l),ppar(l),ipar,iprof,                         k ipar iprof l
     *            ndata(l),(idata(l,i),i=1,ndata(l))                    k ipar iprof l
      write(io,'(10e14.6)')                                             k ipar iprof l
     *            (dgr(idata(l,i)),dgi(idata(l,i)),i=1,ndata(l))        k ipar iprof l
           endif                                                        k ipar iprof l
  101 continue                                                          k ipar iprof l
c                                                                       k ipar iprof
      inew=0                                                            k ipar iprof
c                                                                       k ipar iprof
   10 continue                                                          k ipar iprof


      do 40 k=1,km
      do 40 l=1,nfiltr                                                  k
      read(igk,*)ndat,(dreelr(i),dreeli(i),i=1,ndat)                    k l
      write(io,'(i3,(10e14.6))')ndat,(dreelr(i),dreeli(i),i=1,ndat)     k l
   40 continue                                                          k l
      do 42 k=1,km
      do 42 l=1,nfiltr                                                  k
      read(igk,*)ndat,(dsyntr(i),dsynti(i),i=1,ndat)                    k l
      write(io,'(i3,(10e14.6))')ndat,(dsyntr(i),dsynti(i),i=1,ndat)     k l
   42 continue                                                          k l
      do 43 i=1,ndat
      if(dsyntr(i).ne.str(nfiltr,i).or.dsynti(i).ne.sti(nfiltr,i))      i
     *    stop 'incoherence des donnees synthetiques'                   i
   43 continue                                                          i
      close(igk)
      close(idcdp)
      close(io)
      end
c ----------------------------------------------------------------------
      subroutine lecdp(dp,t,c,npar,ncouch,nper,kmode,lu,ierr)
      parameter(NPERMAX=30,NPARMAX=7,NCOUCHMAX=35,NPMAX=1024,NPMAX2=513)
      parameter(NMODMAX=10)
*****     LECTURE DES DERIVEES PARTIELLES (SAITO PAR EX.)    *****
*     NPERMAX=      nombre maxi de periodes auxquelles existent
*                     les derivees partielles pour un mode isole
*     NPARMAX=     nombre maxi de parametres du modele de
*                     Terre (5 elastiques + densite = 6)
*     NCOUCHMAX=   nombre de couches du modele de Terre
      dimension nper(NMODMAX),t(NPERMAX,NMODMAX), c(NPERMAX,NMODMAX)
      dimension dp(NPERMAX, NCOUCHMAX, NPARMAX,NMODMAX)
      dimension u(NPERMAX,NMODMAX),att(NPERMAX,NMODMAX)
c     dimension ep(NCOUCHMAX),amodel(NPARMAX,NCOUCHMAX)
      do 10 k=1,kmode
      write(0,*) ' Lecture der. part. du mode',k,'  sur',kmode
      read(lu,*) nper(k)
*  LECTURE DES DERIVEES PARTIELLES

c Modif le 05/06/2002 (anne) j'ai ajoute la lecture de u et att pour etre sur
c en plus j'ai change l'enonciation du format de lecture (mis en place
c 1003 format('D',5d16.8) ) Cela est directement tire du format avec
c lequel est ecrit dpa* par convdp2.log10

      do 1 i=1,nper(k)
      read(lu,*,end=9,err=99) t(i,k), c(i,k),u(i,k),att(i,k)
      do 1 ipar=1,npar
      read(lu,1003,end=9,err=99)(dp(i,icouch,ipar,k),
     *                                 icouch=1,ncouch)
1003  format('D',5d16.8)
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
    9 write(0,*)'Fin de fichier detectee dans lecture der. part., mode '
     *           ,kmode
      ierr=1
      return
   99 write(0,*)'Erreur detectee dans lecture der. part., mode ', k
      ierr=1
      return
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
      write(*,'(a)') ' sismogramme reel : ', (bcd(i),i=1,4)
      write(*,'(2a)') '      synthetique = ', entr
      write(*,'(2a)') '      reference   = ', entref
      write(0,'("les fichiers sont-ils corrects?")')
      if(nonoui().ne.1)then
                ierr=1
                return
                endif
1     read(igk,*)km,nfiltr,np,ip,nf1min,nf2max,dnu,dt,dist,ienv
      return
      end
c ----------------------------------------------------------------------
      function nonoui()
      character*10 yo(14),yn(8)
      character c*10
      data yo/'oui','ouai','da','ya','yes','o','y','si','yo','ye','yeeh'
     1,'ok','positif','affirmatif'/
      data yn/'non','n','no','niet','nicht','neni','nein','pas'/
C
C Cette fonction retourne 1 si reponse positive
C                         0 si reponse negative
C                        -1 si reponse incorrecte apres 3 essais
C
      nonoui=-1
      k=0
3     read(*,'(a)')c
      do 1 i=1,14
      if(c.eq.yo(i))nonoui=1
1     continue
      do 2 i=1,8
      if(c.eq.yn(i))nonoui=0
2     continue
      if(nonoui.eq.-1)then
                if(k.ge.3)then
                          write(0,'("advienne ce qu il adviendra")')
                          return
                else
                          k=k+1
                write(0,'("repondez a la question, rappel no ",i2,$)')k
                          go to 3
                          endif
      else 
                return
                endif
      end  
c ----------------------------------------------------------------------
      subroutine sztrans(w,y,n)
      parameter(NPMAX=1024)
      dimension w(NPMAX),y(NPMAX)
c
c Ce sous-programme translate le tableau y,
c l'indice 1 se retrouve en nc=n/2+1
c
      nc=n/2+1
      do 1 i=1,n
1     w(i)=y(i)
      ns=nc-1
      do 2 i=1,ns
      ii=i+ns
2     y(ii)=w(i)
      do 3 i=nc,n
      ii=i-ns
3     y(ii)=w(i)
      return
      end
c ----------------------------------------------------------------
      SUBROUTINE NLOGN(N,xr,xi,SIGN)
      parameter(NPMAX=1024)
C ALGORITHME COOLEY TUKEY.
C      SIGN=-1.0  EXP (-2*I*PI*...)
C      SIGN =1.0 (1/Q) * EXP(2*I*PI*...)
C      NMAX=PLUS GRANDE VALEUR DE N AVEC
C      DIMENSION X(2**N) ,M(NMAX)
      DIMENSION XR(NPMAX),XI(NPMAX),M(20)
      LX=2**N
      DO 1 I=1,N
    1 M(I)=2**(N-I)
      DO 4 L=1,N
      NBLOC=2**(L-1)
      LBLOC=LX/NBLOC
      LBHAF=LBLOC/2
      K=0
      DO 4 IBLOC=1,NBLOC
      FK=K
      FLX=LX
      V=SIGN*6.2831853*FK/FLX
      WKR=COS(V)
      WKI=SIN(V)
      ISTAT=LBLOC*(IBLOC-1)
      DO 2 I=1,LBHAF
      J=ISTAT+I
      JH=J+LBHAF
      QR=XR(JH)*WKR-XI(JH)*WKI
      QI=XI(JH)*WKR+XR(JH)*WKI
      XR(JH)=XR(J)-QR
      XI(JH)=XI(J)-QI
      XR(J)=XR(J)+QR
      XI(J)=XI(J)+QI
    2 CONTINUE
      DO 3 I=2,N
      II=I
      IF (K-M(I)) 4,3,3
    3 K=K-M(I)
    4 K=K+M(II)
      K=0
      DO 8 J=1,LX
      IF (K-J) 5,6,6
    6 HOLDR=XR(J)
      HOLDI=XI(J)
      XR(J)=XR(K+1)
      XI(J)=XI(K+1)
      XR(K+1)=HOLDR
      XI(K+1)=HOLDI
    5 DO 7 I=1,N
      II=I
      IF (K-M(I)) 8,7,7
    7 K=K-M(I)
    8 K=K+M(II)
      IF (SIGN) 11,9,9
    9 DO 10 I=1,LX
      XR(I)=XR(I)/FLX
   10 XI(I)=XI(I)/FLX
   11 RETURN
      END
C-----------------------------------------------------------------------
      subroutine spline(k,nper,ipar,icouch,tt,c,dp,n1,n2,dnu,dpint)
      parameter(NPERMAX=30,NPARMAX=7,NCOUCHMAX=35,NPMAX=1024,NPMAX2=513)
      parameter(NMODMAX=10)
*****    INTERPOLATION DES DERIVEES PARTIELLES         *****
*****    AUX FREQUENCES HARMONIQUES DE LA TF.          *****
*****    On lit des der. part. de la vitesse de phase  *****
*****    et on calcule des d.p. du nbre d'onde / omega *****
*****    1/omega * dk/dp = -1/c**2 * dc/dp             *****
*     NPERMAX=     nombre maxi de periodes auxquelles existent
*                    les derivees partielles pour un mode isole
*     NPARMAX=     nombre maxi de parametres du modele de
*                    Terre (en general: 5 elastiques + densite =6)
*     NCOUCHMAX=   nombre de couches du modele de Terre
*     NPMAX    =   nombre de points maxi des TF
      dimension nper(NMODMAX),tt(NPERMAX,NMODMAX), c(NPERMAX,NMODMAX)
      dimension dp(NPERMAX,NCOUCHMAX,NPARMAX,NMODMAX)
      dimension as(NPERMAX), bs(NPERMAX), cs(NPERMAX), ds(NPERMAX)
      dimension dpint(NPMAX), ibon(NPMAX), t(NPERMAX)
*  LECTURE DES CARACTERISTIQUES DE LA TF
      if(nper(k).lt.5) then
			call flush(6)
                write(0,*) 'Interpolation impossible pour le mode ', k
                write(0,*) '(moins de 5 periodes disponibles)'
                stop ' Fin anormale de dgk12, subroutine <spline>'
                endif
      do 33 i=1,nper(k)
      t(i)=tt(i,k)
      as(i)=(-1./c(i,k)/c(i,k))*dp(i,icouch,ipar,k)
   33 continue
      call inter3(nper(k),t,as,bs,cs,ds)
*  RECHERCHE DES 2 POINTS ENCADRANT LE POINT A INTERPOLER
      if(t(1).gt.t(nper(k))) then
                          i1=1
                          i2=nper(k)
                          inc=1
                else 
                          i1=nper(k)
                          i2=1
                          inc=-1
                endif
      do 6 j=n1,n2
      frq=(j-1)*dnu
      tnu=1./frq
      do 3 i=i1,i2-inc,inc
      a=(t(i)-tnu)*(t(i+inc)-tnu)
      if(a.gt.0.) goto 3
      ibon(j)=i
      if(inc.eq.-1) ibon(j)=i+inc
      goto 4
    3 continue
      ibon(j)=0
    4 continue
*  INTERPOLATION PROPREMENT DITE
      if(ibon(j).eq.0) then
                   dpint(j)=0.
                else
                   dt=tnu-t(ibon(j))
                   ibo=ibon(j)
             dpint(j)=as(ibo)+bs(ibo)*dt+cs(ibo)*dt*dt+ds(ibo)*dt*dt*dt
                endif
    6 continue
C     if(icouch.eq.10) then
C             write(0,*)' Spline: T, C, d.p. Saito pour la couche 10'
C             write(0,*)(t(i),i=1,nper(k))
C             write(0,*)(c(i,k),i=1,nper(k))
C             write(0,*)(dp(i,icouch,ipar,k),i=1,nper(k))
C             nav=(n1+n2)/2
C             f=1./((nav-1)*dnu)
C             write(0,*)'dk/dm interpolee a',f,dpint(nav)
C             endif
      return
      end
c -------------------------------------------------------------
C
      SUBROUTINE INTER3(N,Z,A,B,C,D)
      parameter(NPERMAX=30,NPARMAX=7,NCOUCHMAX=35,NPMAX=1024,NPMAX2=513)
      DIMENSION Z(NPERMAX),A(NPERMAX),B(NPERMAX),C(NPERMAX),D(NPERMAX)
C
      INTEGER P,Q
C     N : NOMBRE DE POINTS.
C     Z : COORDONNEES DES POINTS.
C     A : VALEURS DE LA FONCTION EN CES POINTS (COEFFS CONSTANTS
C         DES POLYNOMES).
C     B,C,D : COEFFS DES POLYNOMES RESPECTIVEMENT DE DEGRE 1,2 ET 3
C     LES COEFFS B,C,D SONT DEFINIS POUR (N-1) POINTS (PAS LE DERNIER)
C
C     SI ZZ EST DANS L'INTERVALLE (Z(I),Z(I+1))
C     F(ZZ) = A(I)+B(I)*(ZZ-Z(I))+C(I)*(ZZ-Z(I))**2+D(I)*(ZZ-Z(I))**3
C       EN PARTICULIER : F(Z(I))= A(I)
C     F,F' ET F'' SONT CONTINUES DE Z(1) JUSQU'A Z(N).
C
C     ALGORITHME CALCUL DES C (ORDRE 2 ) POUR CUBIC SPLINE .
C      PREMIERE PARTIE : A GAUCHE.
      N1=N-1
      N2=N-2
      DO 500 I=2,N1
      I2=I+2
      R1=Z(I)-Z(I-1)
      S1=2.*(Z(I+1)-Z(I-1))
      T1=Z(I+1)-Z(I)
      U1=3.*((A(I+1)-A(I))/(Z(I+1)-Z(I))-(A(I)-A(I-1))/(Z(I)-Z(I-1)))
      IF (I.EQ.2) GO TO 510
C
      R2=Z(I-1)-Z(I-2)
      S2=2.*(Z(I)-Z(I-2))
      T2=(Z(I)-Z(I-1))
      U2=3.*((A(I)-A(I-1))/(Z(I)-Z(I-1))-(A(I-1)-A(I-2))/
     *(Z(I-1)-Z(I-2)))
      R3=-R1*R2/S2
      S3=S1-R1*T2/S2
      U3=U1-R1*U2/S2
      T3=T1
C
      R1=R3
      S1=S3
      T1=T3
      U1=U3
C
      IF (I.EQ.3) GO TO 510
C
      I3=I-3
C     DO 600 NPP=1,I3
      P=I-2-NPP+1
C
      T2=-(Z(P+1)-Z(P))*T2/S2
      U2=3.*((A(P+1)-A(P))/(Z(P+1)-Z(P))-(A(P)-A(P-1))/
     *(Z(P)-Z(P-1)))-(Z(P+1)-Z(P))*U2/S2
      S2=2.*(Z(P+1)-Z(P-1))-(Z(P+1)-Z(P))*R2/S2
      R2=Z(P)-Z(P-1)
C
      R3=-R1*R2/S2
      S3=S1-R1*T2/S2
      U3=U1-R1*U2/S2
      T3=T1
C
      R1=R3
      S1=S3
      T1=T3
      U1=U3
C
600   CONTINUE
  510 continue
C
C    DEUXIEME PARTIE DE L'ALGORITHME : A DROITE.
C
      IF (I.EQ.N1) GO TO 511
      R2=(Z(I+1)-Z(I))
      S2=2.*(Z(I+2)-Z(I))
      T2=Z(I+2)-Z(I+1)
      U2=3.*((A(I+2)-A(I+1))/(Z(I+2)-Z(I+1))-(A(I+1)-A(I))/
     *(Z(I+1)-Z(I)))
C
      R3=R1
      S3=S1-T1*R2/S2
      T3=-T1*T2/S2
      U3=U1-T1*U2/S2
C
      R1=R3
      S1=S3
      T1=T3
      U1=U3
C
      IF (I.EQ.N2) GO TO 511
C
      DO 700 Q=I2,N1
C
C
      R2=-(Z(Q)-Z(Q-1))*R2/S2
      U2=3.*((A(Q+1)-A(Q))/(Z(Q+1)-Z(Q))-(A(Q)-A(Q-1))/
     *(Z(Q)-Z(Q-1)))-(Z(Q)-Z(Q-1))*U2/S2
      S2=2.*(Z(Q+1)-Z(Q-1))-(Z(Q)-Z(Q-1))*T2/S2
      T2=(Z(Q+1)-Z(Q))
C
      R3=R1
      S3=S1-T1*R2/S2
      T3=-T1*T2/S2
      U3=U1-T1*U2/S2
C
      R1=R3
      S1=S3
      T1=T3
      U1=U3
C
700   CONTINUE
C
511   C(I)=U1/S1
C
500   CONTINUE
C
C      FIN DE L'ALGORITHME DE CALCUL DES COEF C ......
C
100   CONTINUE
C
      C(1)=0.
      C(N)=0.
      B(1)=(A(2)-A(1))/(Z(2)-Z(1))-C(2)*(Z(2)-Z(1))/3.
      D(1)=C(2)/(3.*(Z(2)-Z(1)))
C
C    CAL CUL DES COEFFS B ET A.
C
      DO 200 K=3,N1
C
      B(K-1)=(A(K)-A(K-1))/(Z(K)-Z(K-1))-C(K-1)*(Z(K)-Z(K-1))
     *-(C(K)-C(K-1))*(Z(K)-Z(K-1))/3.
C
      D(K-1)=(C(K)-C(K-1))/(3.*(Z(K)-Z(K-1)))
200   CONTINUE
C
      B(N-1)=(A(N)-A(N-1))/(Z(N)-Z(N-1))-2.*C(N-1)*(Z(K)-Z(K-1))/3.
      D(N-1)=-C(N-1)/(3.*(Z(N)-Z(N-1)))
C
      RETURN
      END
c ----------------------------------------------------------------------
      subroutine dpenv(l,gr,gi,ndata,idata,str,sti)
c     parameter(NFILTRMAX=5,NPMAX=1024,NDATAMAX=5)
c     modif eric Octobre 97 inversion plusieurs lobes:NDATAMAX=10
c     parameter(NFILTRMAX=5,NPMAX=1024,NDATAMAX=10)

c     modif eric 10 Oct 2000 :NFILTRMAX=10
      parameter(NFILTRMAX=10,NPMAX=1024,NDATAMAX=10)
      dimension idata(NFILTRMAX,NDATAMAX), gr(NPMAX), gi(NPMAX),
     *          str(NFILTRMAX,NDATAMAX),sti(NFILTRMAX,NDATAMAX)
c
c     Calcul des derivees partielles de l'enveloppe
c     a partir des derivees partielles de gk^, des valeurs de
c     gk^, et des valeurs du module de gk^.
c
      do 1 i=1,ndata
      j=idata(l,i)
      gr(j)=(str(l,i)*gr(j) + sti(l,i)*gi(j)) /
     *        sqrt(str(l,i)*str(l,i)+sti(l,i)*sti(l,i))
    1 continue
      return
      end
c ----------------------------------------------------------------------
      subroutine dpphase(l,gr,gi,ndata,idata,str,sti)
	parameter(echpha=1.e4)
c     parameter(NFILTRMAX=5,NPMAX=1024,NDATAMAX=5)
c     modif eric Octobre 97 inversion plusieurs lobes:NDATAMAX=10
c     parameter(NFILTRMAX=5,NPMAX=1024,NDATAMAX=10)

c     modif eric 10 Oct 2000 : NFILTRMAX=10
      parameter(NFILTRMAX=10,NPMAX=1024,NDATAMAX=10)
      dimension idata(NFILTRMAX,NDATAMAX), gr(NPMAX), gi(NPMAX),
     *          str(NFILTRMAX,NDATAMAX),sti(NFILTRMAX,NDATAMAX)
c
c     Calcul des derivees partielles de la phase
c     a partir des derivees partielles de gk^, des valeurs de
c     gk^, et des valeurs du module de gk^.
c
      do 1 i=1,ndata
      j=idata(l,i)
      gr(j)=(str(l,i)*gi(j) - sti(l,i)*gr(j)) /
     *        (str(l,i)*str(l,i)+sti(l,i)*sti(l,i))
c     Mise a l'echelle numerique
	gr(j)=echpha*gr(j)
    1 continue
      return
      end
c ----------------------------------------------------------------
      subroutine szgauss(l,npp,nf1,nband,ff)
c     parameter(NPMAX2=513,NFILTRMAX=5)

c     Modif eric 10 Oct 2000 NFILTRMAX=10
      parameter(NPMAX2=513,NFILTRMAX=10)
      dimension f(NPMAX2), ff(NFILTRMAX,NPMAX2)
c
c Calcul une fonction de Gauss tronquee a -30 db et prolongee par des
c zeros a droite et a gauche.
c
c*** npp est le nombre total de points
c*** nf1 est l'indice du premier point non nul
c*** nband est la largeur totale de la gaussienne
c*** f contient le resultat
c
c Cara mars 1986
c
      n=nf1-1
      do 1 i=1,n
1     f(i)=0.
      n=nf1+nband
      if(n.gt.npp.or.nf1.gt.npp) write(0,'("erreur dans les arguments
     1 de szGauss")')
      do 2 i=n,npp
2     f(i)=0.
      n=nband/2.0001
      beta=3.*alog(10.)/n/n/2
      n=nf1+n
      nn=n-1
      do 3 i=nf1,nn
      d=n-i
3     f(i)=exp(-beta*d*d)
      f(n)=1.
      n=n+1
      nn=nf1+nband-1
      i0=(n-1)*2
      do 4 i=n,nn
      ii=i0-i
4     f(i)=f(ii)
      do 5 i=1,npp
    5 ff(l,i)=f(i)
      return
      end
