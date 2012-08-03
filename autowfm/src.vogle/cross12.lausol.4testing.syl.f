c Version de cross12.lr.auto.F faite pour tourner avec la bibliotheque graphique Vogle
c  Version pour LINUX ( 04/06/2002)
 	parameter(echpha=1.e4)
      parameter (ndim=4096,nmode=10,nfdim=10,nfr=2049,invd=10)
c     parameter (ndim=20000,nmode=10,nfdim=10,nfr=2049,invd=10)
      parameter (nlob=10)

C  nfr = ndim/2 +1
      dimension xr(ndim),xi(ndim),yr(ndim),yi(ndim),xrs(nfr),
     1xis(nfr),f(nfr),yrm(nfr),yim(nfr),n1(nfdim),nb(nfdim),
     1ppar(nfdim),pper(nfdim),xsave(ndim)
C     dimension work(ndim),workk(ndim),idata(invd)
      dimension work(ndim),workk(ndim),idata(invd)
      dimension qualit(nmode,nfdim,nlob),ndatam(nmode,nfdim),
     *          idatam(nmode,nfdim,invd),qual(nlob),qualo(nmode,nfdim)
      dimension ey(nmode,nfdim)
      dimension ratiosynt(nmode,nfdim)
      character sort*30,entr*30,entref*30,bcd(4)*80,titre*80,tit1*80
      real*8 pi2
      character stockage*80
      data pi2 /6.28318530718/
      integer echy,npts
      common/c1/ratiosynt
      common/c2/imodmax,nperpi                
c
c   Programme de comparaison graphique entre synthetique multimode
c   et sismogramme reel:
c                        soit au niveau du signal
c                        soit au niveau du l'intercorrelogramme 
c                        soit au niveau de l'intercorrelogramme filtre
c
c Cara avril 1986
c
c	cross9.f : version utilisant et dessinant la partie reelle de l'int.
c	cross10.f: version utilisant et dessinant la phase de l'intercorr.
c	cross11.f: version inversant de plus le moment sismique
c--Modif eric 1992----------------------------------------------
c	cross12.lr.f: version de cross12.f modifiee pour l'application
c                   a plusieurs trajets.
c--Modif eric 1997----------------------------------------------
c     un systeme de piquage automatique des intercorrelogrames a ete
c     developpe est et introduit dans cette version de l'algo.
c     cross12.lr.anew.F
c     Dans cross12.lr.amorelobes.F on selectione automatiquement
c     de 2 lobes sur les intercorrelogrames. On peut aussi introduire 
c     une erreur differente sur chacun des lobes, dependant du rapport
c     signal/bruit (changement du tabelau qualit, introduction
c     du tableau qual). 
c     Dans cross12.lr.amltest.F on a change par rapport a
c     cross12.lr.amorelobes.F :
c     la selection des lobes : un lobe n'est pas selectione
c     s'il est trop proche de la bordure du cross-correlograme
c     (et trop loin du temps de reference).
c     Derniere Version : cross12.lr.auto.F
c     On utilise le fichier energiemode pour reperer
c     a chacune des periodes ou on va filtrer les intercorrelogrammes 
c     le mode le plus energetique ainsi que la contribution
c     de chaque mode a l'energie totale du sismogramme. 
c      
c--------------------------------------------------------------
      inr=7
      ist=8
      ioo=9
      io=10
      is=11
      iref=12
      i1=0
      icross=0
      isegment=1
c     write(0,*) 'L\'option dessin est necessaire pour :'
c     write(0,*) '    - calcul de l\'enveloppe'
c     write(0,*) '    - pickage et stockage en vue du calcul des d.p.'
      write(*,'("voulez-vous dessiner? ",$)')
c
c Initialisation du dessin eventuel
c
      if(nonoui(0).eq.1)then
      ides=1
      xg=0.
      xd=100.
      dx=xd-xg
      yt=100.
      y0=0.
      iview=7
      input=0
      idim=0
      ibase=0
      icoul=8
      call initdes(xg,xd,y0,yt,iview,input,idim,icoul,isegment)
	call CreateRetainSeg(1)
      yt=yt-10.
      endif
c
c  Preparation d'une inversion des intercorrelogrammes
c
98    ider=0
      ipic=0
      write(*,'("voulez-vous piquer les valeurs de l\'intercorrelogramme
     1? ",$)')
      if(nonoui(0).eq.1)then
      write(*,'("entrer coef (1 est une valeur standard): ",$)')
      read(*,*)coef
      ipic=1
      write(*,'("voulez-vous creer un fichier pour calcul des derivees p
     1artielles? ",$)')
         if(nonoui(0).eq.1)then
         ider=1
         if(ider.eq.1.and.ides.eq.0)stop 'il fallait dessiner!...'
         open(io,form='unformatted',status='scratch')
         write(*,'("nom du fichier? ",$)')
         read(*,'(a)')sort
         open(ioo,status='new',err=98,file=sort)
         endif
      endif
c
c Lecture des trois fichiers d'entree
c
1     call szopen(is,'("fichier contenant le sismogramme reel? ",$)')
      read(is,'(a)')(bcd(i),i=1,4)
c     write(*,'(a)')(bcd(i),i=1,4)
      write(*,'("est-ce le bon sismogramme? ",$)')
      if(nonoui(0).eq.0)go to 1
      read(is,'(2i3,f6.2)')iho,imo,so
      read(is,'(2i3,f6.2,a)')ihd,imd,sd,tit1
      read(is,*) n,pas,distr
      if(n.gt.ndim)go to 99
        titre=bcd(1)
        titre=titre(1:28)//'     '//tit1(4:16)
        write(tit1,'(i5.2,a,i2.2)') iho,'h',imo
        titre=titre(1:48)//tit1(1:14)
      td=(ihd-iho)*3600.+(imd-imo)*60.+sd-so
      read(is,*)(xr(i),i=1,n)
      close(is)
c
22    write(*,'("fichier contenant le signal synthetique? ",$)')
      read(*,'(a)')entr
      open(is,status='old',form='unformatted',err=22,file=entr)
      read(is) npp,dnu,t1,dt,dlt,km,ip,dists
      np=2*(npp-1)
      nppp=npp+1
c
23    write(*,'("fichier contenant le signal de reference? ",$)')
      read(*,'(a)')entref
      open(iref,status='old',form='unformatted',err=23,file=entref)
      read(iref) nppref,dnuref,t1ref,dtref,dltref,kmref,ipref,distref
      if(npp.ne.nppref.or.dnu.ne.dnuref.or.t1.ne.t1ref.or.
     *   dt.ne.dtref.or.dlt.ne.dltref.or.km.ne.kmref.or
     *   .ip.ne.ipref.or.dists.ne.distref)then
      write(*,'(" incoherence des fichiers synt. et ref.")')
      close(iref)
      close(is)
      go to 22
      endif
      if(abs(distr-distref).gt.0.5)stop" incoherence dans les distances 
     *epicentrales"
      if(npp.gt.ndim/2+1)go to 99
c
c Ecriture des titres du sismogramme reel et des noms des fichiers
c synthetique et reference dans les fichiers contenant les
c intercorrelogrammes temporels et frequentiels
c
      if(ider.eq.1)then
      write(ioo,'(a)')(bcd(i),i=1,4)
      write(ioo,'(a)')entr
      write(ioo,'(a)')entref
      endif
c
      if(ides.eq.1)then
      dhc=(yt-y0)/4./km
      dhs=(yt-y0)/4
      endif
c
c Normalisation du sismogramme enregistre
      call sznorm(n,pas,td,npp,dt,t1,xr,xi,ip,dnu)

c--------------------------------------------------------
c 8 Janvier 1998---version automatique du piquage
c Les periodes de filtrage des intecorrelogrammes ont ete
c choisies lors de la selection des donnees. Elles sont
c supposees etre ecrites dans 2 fichiers appelles
c normalement "pick.ray.1" pour rayleigh et "pick.lov.1"
c pour love.
c Les periodes de filtrage, largeur de bande sont lues 
c dans le fichier pick.
c les donnees de l'inversion et les parametres relatifs a leur
c stockage sont selectionnees automatiquement dans la subroutine
c des et ce a chaque run du pgm.
c ici le parametre lecpic sert a reperer l'iteration 0.
c qd lecpic=0 (iteration 0 du pgm) on va recuperer
c le "facteur d'echelle" qui est utilise comme critere
c de rejet d'un sismogramme dans la procedure automatique
c On ecrit les parametres relatifs au stockage des 
c intercorrelogrames a la suite du fichier pick.
c Ceci n'a aucune utilite du point de vue informatique
c mais permet de tester visuellement quelle erreur est donnee
c a telle ou telle donnee de l'inversion.
c--------------------------------------------------------
97    write(*,*)'entrez le nom du fichier pick contenant les periodes'
      write(*,*)'de filtrage des intercorrelogrammes.'
      read (*,*) stockage
      write(*,*)'Entrer lecpic :'
      write(*,*)'-lecpic=0 : 1er run'
      write(*,*)'-lecpic=1 : 2eme run'
      read(*,*)lecpic
      open(ist,status='old',err=97,file=stockage)
      read(ist,*)nmodpi,nperpi
c--------------------------------------------------------
c si echy=1 les cross-corr d'une colonne seront a l'echelle du 
c fondamental (comme dans le wfm classique), sinon, pour chaque 
c mode il seront a l'echelle de la premiere frequence pour
c laquelle ce mode est filtre.
c--------------------------------------------------------
      write(*,*)'voulez-vous les cross-corr d\'une colonne a l\'echelle 
     *de celle du fondamental(1=oui)?'
      read (*,*)echy
c choix du menu
c
c NB: nonoui(1) va lire la reponse dans le fichier pick
c     nonoui(0) prend la reponse au clavier
c
      if(ipic.eq.0) then
      write(*,'("voulez-vous calculer les intercorrelogrammes? ",$)')
      if(nonoui(0).ne.1) go to 7
      endif
      icross=1
      if(ides.eq.1) then
      ienv=0
      write(*,'("voulez-vous dessiner l\'enveloppe? ",$)')
      if(nonoui(1).eq.1) then
      ienv=1
      endif
c Si on ne dessine pas l'enveloppe et que l'on lit la reponse a l'ecran
c on ne l'ecrit pas dans le fichier de stockage(ce cas n'etant jamais
c rencontre en waveform modelling).
      endif
      pper(1)=0.
      ppar(1)=0.
      write(*,'("voulez-vous filtrer les intercorrelogrammes? ",$)')
      if(nonoui(1).ne.1) go to 7
c cas ou on filtre les intercorrelogrammes
c
c Calcul des parametres des filtres passe-bande gaussien
c pour les intercorrelogrammes
c

      eps=1.e-6
      nf1min=npp
      nf2max=0
13    write(*,'("entrer la periode centrale en sec: ",$)')
      read(ist,*) per
      write(*,*)'PERIODE', per
3     write(*,'("entrer la largeur de bande relative a 30db: ",$)')
      read(ist,*) par
      write(*,*)'PAR', par
      fr=1./per/dnu
      ir=(fr*par+1.5)/2
      nband=2*ir+1
      nf1=fr+1.5-float(ir)
      if(nf1.lt.1.or.nf1+nband.gt.npp)go to 3
      i1=i1+1
      ppar(i1)=par
      pper(i1)=per
      n1(i1)=nf1
      nb(i1)=nband
      write(*,*)'NB',nb(i1),' periode',i1
      nf2=nf1+nband-1
      if(nf1.lt.nf1min) nf1min=nf1
      if(nf2.gt.nf2max) nf2max=nf2
      write(*,'("encore",i2," periodes possibles, on continue? ",$)')nf
     1dim-i1
      nonououi=nonoui(1)
      write(*,*)'ENCORE PERIODE',nonououi
      if (nonououi.eq.1.and.i1.lt.nfdim) then
      go to 13
      endif
      dx=(xd-xg)/i1
7     continue
c modif eric 4 Juilllet 1997
c     write(*,*)'derniere periode entree=phase? (1=oui, 0=non)'
      write(*,*)'Combien de periodes=phases ? '
      read(ist,*) ipha
      write(*,*) 'phase=',ipha

      nper=max(i1,1)
      if(ider.eq.1)write(io)km,nper
      if(ider.eq.1)
     *        write(ioo,*)km,nper,np,ip,nf1min,nf2max,dnu,dt,dists,ienv
      do 4 i=1,npp
      xrs(i)=0.                                                         i
4     xis(i)=0.                                                         i
c
c Calcul de la TF de la trace synthetique somme
c
      do 111 k=1,km
      read(is)(yrm(i),yim(i),i=1,npp)                                   k
      if(ider.eq.1) write(ioo,'(10e14.6)')                              k
     *              (yrm(i),yim(i),i=nf1min,nf2max)                     k
      do 9 i=1,npp                                                      k
      xrs(i)=xrs(i)+yrm(i)                                              k i
9     xis(i)=xis(i)+yim(i)                                              k i
111   continue                                                          k
      close(is)

c
      if(icross.eq.0) go to 17
c-modif eric 8 Janvier 98-----------------------
c ratiosyntmax sert a evaluer le mode le plus
c energetique a longue periode(ie 60 ou 80s).
c c'est ce mode qui sera inverse lorsque les 
c periodes lp seront selectionees.
c lecture du fichier energiemode
      ratiosyntmax=0
      open(inr,status='old',file='energiemode')
      do 71 imodpi=1,nmodpi
      read(inr,*)ibid,(bid,ratiosynt(imodpi,iper),iper=1,nperpi)
      write(*,*)'RATIOSYNT'
      write(*,*)(ratiosynt(imodpi,iper),iper=1,nperpi)
      if (ratiosynt(imodpi,1).gt.ratiosyntmax) then
       ratiosyntmax=ratiosynt(imodpi,1)
       imodmax=imodpi
      endif
71     continue
      write(*,*)'IMODMAX',imodmax
      close(inr)
c--------------------------------------------
c
c Calcul des intercorrelogrammes avec les branches de reference
c
      do 10 k=1,km
      read(iref)(yrm(i),yim(i),i=1,npp)                                 k
c                                                                       k
c TF des intercorrelogrammes reels gk(per,t) ou gk(t)                   k
c                                                                       k
      do 10 l=1,nper                                                    k
      do 8 i=1,npp                                                      k l
      yr(i)=xr(i)*yrm(i)+xi(i)*yim(i)                                   k l i
8     yi(i)=xi(i)*yrm(i)-xr(i)*yim(i)                                   k l i
      if(i1.ne.0) then                                                  k l
      nf1=n1(l)                                                         k l
      nband=nb(l)                                                       k l
      call szgauss(npp,nf1,nband,f)
      do 6 i=1,npp                                                      k l
      yi(i)=yi(i)*f(i)                                                  k l i
6     yr(i)=yr(i)*f(i)                                                  k l i
      endif                                                             k l
c                                                                       k l
c Dessin et stockage des intercorrelogrammes reels                      k l
c                                                                       k l
      do 5 i=nppp,np                                                    k l
      yi(i)=0.                                                          k l i
5     yr(i)=0.                                                          k l i
      call nlogn(ip,yr,yi,+1.)                                          k l
      call sztrans(work,yr,np)                                          k l
      call sztrans(work,yi,np)                                          k l
      if(ides.eq.1) then                                                k l
      yb=(yt-y0)/2+(k-1)*2*dhc                                          k l
      yh=yb+dhc*0.8                                                     k l
      ey(k,l)=0.
                                                            
****    Si echy =1, on trace tous les cross-corr. d'une colonne
****    a la meme echelle que celle du fondamental.
      if (echy.eq.1) ey(k,l)=ey(1,l)

****      JJL, 28 Septembre 1988                                        k l
      xgg=xg+(l-1)*dx+dx*0.2                                            k l
      xdd=xgg+dx*0.8                                                    k l
c        ecriture des cross-corr. reels sur le fichier scratch          k l
      if(ider.eq.1) write(io)np,dt,pper(l),ppar(l),                     k l
     1(yr(i),yi(i),i=1,np)                                              k l
c                                                                       k l
      ibase=ienv
      ndata=0
c modif eric inversion automatique
      if((ipha.eq.1.and.l.eq.nper) 
     *.or.(ipha.eq.2.and.l.eq.(nper-1))
     *.or.(ipha.eq.2.and.l.eq.(nper)))then
      npts=1
      else
      npts=3
      endif
      call des(yi,yr,np,xgg,xdd,yb,yh,ey(k,l),dt,ibase,idata,ndata,pper, 
     *ppar,k,l,ipic,coef,qual,ienv,lecpic,npts,isegment) 
      ndatam(k,l)=ndata
      qualo(k,l)=qual(1 )                                               k l
      ilob=1
      do 24 i=1,ndata                                                   k l
      qualit(k,l,i)=qual(ilob)                                          k l
      if(mod(i,npts).eq.0)ilob=ilob+1
c     write(*,*)'QUALIT,mode,periode,indice',k,l,i,qualit(k,l,i)
24    idatam(k,l,i)=idata(i)                                            k l i
      endif                                                             k l
10    continue
      read(iref) amoment
      close(iref)
c
c Calcul,dessin et stockage des intercorrelogrammes synthetiques
c gk^(per,t) ou gk^(t)
c
      open(iref,status='old',form='unformatted',file=entref)
      read(iref) npp,dnu,t1,dt,dlt,km,ip
      do 20 k=1,km
c     lecture et ecriture des modes de reference
      read(iref)(yrm(i),yim(i),i=1,npp)                                 k
      if(ider.eq.1) write(ioo,'(10e14.6)')                              k
     *              (yrm(i),yim(i),i=nf1min,nf2max)                     k
      do 20 l=1,nper                                                    k
      do 18 i=1,npp                                                     k l
      yr(i)=xrs(i)*yrm(i)+xis(i)*yim(i)                                 k l i
18    yi(i)=xis(i)*yrm(i)-xrs(i)*yim(i)                                 k l i
      if(i1.ne.0) then                                                  k l
      nf1=n1(l)                                                         k l
      nband=nb(l)                                                       k l
      nf2=nf1+nband-1                                                   k l
      call szgauss(npp,nf1,nband,f)                                     k l
      do 16 i=1,npp                                                     k l
      yi(i)=yi(i)*f(i)                                                  k l i
16    yr(i)=yr(i)*f(i)                                                  k l i
      endif                                                             k l
      if(ider.eq.1)  write(ioo,*)npp,nf1,nband,pper(l),ppar(l)          k l
      do 15 i=nppp,np                                                   k l
      yi(i)=0.                                                          k l i
15    yr(i)=0.                                                          k l i
      call nlogn(ip,yr,yi,+1.)                                          k l
      call sztrans(work,yr,np)                                          k l
      call sztrans(work,yi,np)                                          k l
c                                                                       k l
      if(ides.eq.1) then                                                k l
      yb=(yt-y0)/2+(k-1)*2*dhc+dhc*0.8                                  k l
      yh=yb+dhc*0.8                                                     k l
****  On trace le cross-correl. synth. a la meme echelle que le reel    k l
****  ey=0.              (a dater du 20 Septembre 1988). JJL.           k l
      xgg=xg+(l-1)*dx+dx*0.2                                            k l
      xdd=xgg+dx*0.8                                                    k l
      do 199 i=1,np                                                     k l
      workk(i)=yi(i)
199   work(i)=yr(i)                                                     k l i
c        ecriture des cross-corr. synt. sur le fichier scratch          k l
      if(ider.eq.1) write(io)np,dt,pper(l),ppar(l),                     k l
     1(yr(i),yi(i),i=1,np)                                              k l
C     if(ider.eq.1)then                                                 k l
C     write(ioo,'(10e14.6)')'provisoire',np,(work(i),i=1,np)            k l
C     endif                                                             k l
*     call SetLinestyle(DOTTED)                                         k l
      ienv0=ienv
      ndata=ndatam(k,l)
      eyy=ey(k,l)
      if(ienv.eq.1.and.ndatam(k,l).eq.1.and.idatam(k,l,1).ne.0) then
      	eyy=ey(k,l)/2.
      	endif
      ibase=ienv0
      ipic0=0                                                           k l
	idata(1)=idatam(k,l,1)
      call des(workk,work,np,xgg,xdd,yb,yh,eyy,dt,ibase,idata,ndata,pper   
     *,ppar,k,l,ipic0,coef,qual,ienv0,0,npts,isegment) 
      call SetLinestyle(SOLID)                                          k l
      if(ider.eq.1) then                                                k l
c        Ecriture des cross-corr. synth. aux memes indices              k l
c        que ceux pickes sur les cross-corr. reels                      k l
        ndata=ndatam(k,l)                                               k l
        do 26 i=1,ndata                                                 k l
26      idata(i)=idatam(k,l,i)                                          k l i
        write(ioo,*)ndata,(idata(i),i=1,ndata)                          k l
        if(idata(1).ne.0) then                                          k l
                write(ioo,'(i3,(10e14.6))')                             k l
     *              ndata,(yr(idata(i)),yi(idata(i)),i=1,ndata)         k l
                else                                                    k l
*               cas ou on n'a rien picke                                k l
                write(ioo,*)ndata,'   0.   0.    pas de donnee'         k l
                endif                                                   k l
        endif                                                           k l
      endif                                                             k l
20    continue                                                          k l
      close(iref)
      if(ider.ne.1) go to 25
      rewind(io)
      read(io)km,nper
c   Lecture des cross-correlogrammes reels sur le fichier scratch
      do 40 k=1,km
      do 40 l=1,nper                                                    k
      read(io)np,dt,pper(l),ppar(l),(yr(i),yi(i),i=1,np)                k l
      ndata=ndatam(k,l)                                                 k l
      do 41 i=1,ndata                                                   k l
      workk(i)=yi(idatam(k,l,i))                                        k l i
41    work(i)=yr(idatam(k,l,i))                                         k l i
      if(idatam(k,l,1).ne.0) then                                       k l
                write(ioo,'(i3,(10e14.6))')                             k l
     *          ndata,(work(i),workk(i),i=1,ndata)                      k l
                else                                                    k l
*               cas ou on n'a rien picke                                k l
                write(ioo,*)ndata,'   0.   0.    pas de donnee'         k l
                endif                                                   k l
40    continue                                                          k l
c   Lecture des cross-correlogrammes synt. sur le fichier scratch
      do 42 k=1,km
      do 42 l=1,nper                                                    k
      read(io)np,dt,pper(l),ppar(l),(yr(i),yi(i),i=1,np)                k l
      ndata=ndatam(k,l)                                                 k l
      do 43 i=1,ndata                                                   k l
      workk(i)=yi(idatam(k,l,i))                                        k l i
43    work(i)=yr(idatam(k,l,i))                                         k l i
      if(idatam(k,l,1).ne.0) then                                       k l
                write(ioo,'(i3,(10e14.6))')                             k l
     *          ndata,(work(i),workk(i),i=1,ndata)                      k l
                else                                                    k l
*               cas ou on n'a rien picke                                k l
                write(ioo,*)ndata,'   0.   0.    pas de donnee'         k l
                endif                                                   k l
42    continue                                                          k l
      rewind(io)
c   Lecture des cross-correlogrammes reels sur le fichier scratch
c   en vue de calculer l'erreur sur les donnees a inverser
      read(io)km,nper
      write(ioo,*)'Erreurs-donnees estimees visuellement:'
      do 440 k=1,km
      do 440 l=1,nper                                                   k
      read(io)np,dt,pper(l),ppar(l),(yr(i),yi(i),i=1,np)                k l
      ndata=ndatam(k,l)                                                 k l
      do 441 i=1,ndata                                                  k l
      if(ienv.eq.1.and.ndata.ne.1) then                                 k l i
        ybidr=yr(idatam(k,l,i))                                         k l i
        ybidi=yi(idatam(k,l,i))                                         k l i
        work(i)=sqrt(ybidr*ybidr+ybidi*ybidi)*qualit(k,l,i)               k l i
        else                                                            k l i
c	pour une donnee de phase, on choisit l'erreur
c	en centiemes de periode
        work(i)=echpha*pi2*qualit(k,l,i)                                    k l i
        endif                                                           k l i
441   continue                                                          k l i
      if(idatam(k,l,1).eq.0) then                                       k l
          work(1)=1.e5                                                  k l
          endif                                                         k l
      write(ioo,'(i3,(10e14.6))')ndata,(work(i),i=1,ndata)              k l
440   continue                                                          k l
25    continue
17    continue
      if(ides.ne.1)go to 33
c
c Dessin des traces temporelles reelles et synthetiques
c
      call fildes(xr,xi,npp,dnu)
      call nlogn(ip,xr,xi,+1.)
      yb=y0
      dhs=dhs*2./3.
      yh=yb+dhs

c      calcul et sauvegarde du signal reel, qui va etre modifie dans des(..)
      ienv=0
      do 30 i=1,np
      xr(i)=2*xr(i)                                                     i
30    xsave(i)=xr(i)                                                    i
	eyloc=0.
      ibase=0
      call des(xi,xr,np,xg,xd,yb,yh,eyloc,dt,ibase,idata,ndata,pper,
     *ppar,1,1,0,coef,qual,ienv,0,npts,isegment)
      eyreel=eyloc
      do 11 i=nppp,np
      xr(i)=0.                                                          i
11    xi(i)=0.                                                          i
      do 12 i=1,npp
      xr(i)=xrs(i)                                                      i
12    xi(i)=xis(i)                                                      i
      call fildes(xr,xi,npp,dnu)
      call nlogn(ip,xr,xi,+1.)
      yb=yh
      yh=yb+dhs
      do 31 i=1,np
  31  xr(i)=2*xr(i)                                                     i
c modif eric 8 Aout 1997: on calcule l'energie entre
c 3.5 et 6 km/s
      smreel=0.
      smresi=0.
      do 32 i=1,np
      vtest=(dists/(t1+(i-1)*dt))
c     write(*,*)'VTEST',vtest 
      if(vtest.gt.3.5.and.vtest.lt.6)    
     *smreel=smreel+xsave(i)*xsave(i)                                   i
      xsave(i)=xsave(i)-xr(i)                                           i
      if(vtest.gt.3.5.and.vtest.lt.6)    
     *smresi=smresi+xsave(i)*xsave(i)                                   i
32    continue
      eyloc=0.
      ibase=0
      call des(xi,xr,np,xg,xd,yb,yh,eyloc,dt,ibase,idata,ndata,pper,
     *ppar,1,1,0,coef,qual,ienv,0,npts,isegment)
      write(*,*)' Echelle de la trace reelle:',eyreel
      write(*,*)' Echelle de la trace synth.:',eyloc
      amult=eyloc/eyreel
      yb=yh
      yh=yb+dhs
      enrgrs=smresi/smreel
      eyloc=eyreel
      ibase=0
      call des(xi,xsave,np,xg,xd,yb,yh,eyloc,dt,ibase,idata,ndata,
     *pper,ppar,1,1,0,coef,qual,ienv,0,npts,isegment)
      write(*,*)' Echelle de la trace residu:',eyloc
      write(*,*)' Energie relative de 3.5 a 6 km/s(residu/reel):',enrgrs
      if(ider.eq.1)
     *write(ioo,*)' Energie relative de 3.5 a 6 km/s(res./reel):',enrgrs
c-----------------------------------------------------------------
c modif eric juillet 97: si lepic=0 on ecrit amult a la fin du fich.
c de sortie amult sera utilise lors de la selection des donnees.
      if(lecpic.eq.0) write(ioo,*)'facteur d\'echelle',amult
      write(*,*)
      write(*,*)' Moment sismique initial:',amoment
      write(*,*)
      if(lecpic.eq.0.and.(amult.gt.5.or.amult.lt.0.2)) then
      open(13,file="ECH")
      write(13,*)amult
      close(13)
      endif
      write(*,*) "finished if statement"
c modif eric pour le dessin on n'indiquera seulement l'erreur
c Cd qualo du premier point selectione.
      call desu(np,xg,xd,y0,yt,dists,t1,dt,dhs,entr,
     *                  titre,nper,pper,km,qualo,ey,amult)
      close(ist)

      write(*,*) "ended subr. desu"
c
c Fin
c
33    call flush(6)
         call CloseRetainSeg(isegment)
      if(ider.eq.1) then
      close(ioo)
      close(io)
      endif
      write(*,*)
      write(*,'(" voulez-vous un autre sismogramme? ",$)')
      call flush(6)
      if(mouinon().eq.1)then
*dess       if(nonoui().eq.1)then
      i1=0
      icross=0
      do 34 i=1,nfdim
      n1(i)=0                                                           i
      nb(i)=0                                                           i
      ppar(i)=0.                                                        i
34    pper(i)=0.                                                        i
      if(ides.eq.1)then
      ibase=0
      dx=xd-xg
      call DelAllRetainSegs()
      call NewFrame
      call CreateRetainSeg(1)
      endif
      go to 98
      endif
      if(ides.eq.1) then
	call CloseRetainSeg(0)
	call findes
	endif
	if(ider.eq.1)write(0,*)' Mettre a jour tv.xx (par gktv par ex.)'
      stop
99    write(*,'("depassement de dimension")')
      if(ides.eq.1)call findes
      stop
      end
c ----------------------------------------------------------------
      subroutine szgauss(npp,nf1,nband,f)
      dimension f(1)
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
      if(n.gt.npp.or.nf1.gt.npp) write(*,'("erreur dans les arguments
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
      return
      end
c --------------------------------------------------------------
      subroutine sznorm(n,pas,td,npp,dt,t1,xr,xi,ip,dnu)
      dimension xr(1),xi(1)
c
c   Interpolation du sismogramme reel aux instants d'echantillonnage
c   du signal synthetique, transformee de Fourier et restitution
c   de la TF du signal analytique associe, calculee sur le meme 
c   nombre de points et aux memes frequences que le signal synthetique
c
c                 n    nombre de point du signal reel
c                 pas  pas en temps du signal reeel
c                 td   instant du debut du signal reel
c                 npp  indice de la frequence de Nyquist du signal
c                      synthetique
c                 t1   instant du debut du signal synthetique
c                 xr(k)tableau contenant en entree le signal reel
c                      en sortie la partie reele de la TF du 
c                      signal analytique associe
c                 xi(k)tableau contenant en sortie la partie
c                      imaginaire de la TF du signal analytique
c                      k=1,2**ip
c                 ip   puissance de 2 superieure au nombre de points
c                      du signal synthetique
c                 dnu  pas en frequence du synthetique
c
      np=2*(npp-1)
      eps=0.000001
      if(pas.gt.dt) go to 10
c
c Filtrage passe bas
c
c     write(*,'("filtrage anti-repliement du sismogramme reel:")')
      call szpu2(n,nnp,iip)
      dnn=1./pas/nnp
      if(n.lt.nnp)then
      nn=n+1
      do 9 i=nn,nnp
9     xr(i)=0.
      endif
      do 8 i=1,nnp
8     xi(i)=0.
      call nlogn(iip,xr,xi,-1.)
      per=1./dnu/(npp-1)
c     write(*,'("la periode de coupure du filtre est ",f10.4)')per
      if2=1./per/dnn+1.5
      do 17 i=if2+1,nnp
      xr(i)=0.
17    xi(i)=0.
      xr(1)=xr(1)/2
      call nlogn(iip,xr,xi,+1.)
      do 6 i=1,n
      xr(i)=2*xr(i)
6     xi(i)=0.
10    continue
c
c Interpolation parabolique
c
      dt1=t1-td
      dt2=(n-1)*pas+td-(np-1)*dt-t1
      a1=abs(dt1)
      a2=abs(dt2)
      L1=a1/dt-eps
      L2=a2/dt-eps
      aps=eps*pas
      if(a1.le.aps) then
      m1=2
      else
      if(dt1.gt.aps) then
      m1=1
      else
      m1=L1+2
      endif      
      endif      
      if(a2.le.aps)then
      m2=np-1
      else
      if(dt2.gt.aps) then
      m2=np
      else
      m2=np-L2-1
      endif
      endif
      print*,m1,m2,np,ip,npp
      do 1 j=m1,m2
      t=t1+(j-1)*dt
      i2=(t-td)/pas+1.5
      i1=i2-1
      i3=i2+1
      call szint(i1,i2,i3,pas,xr,t,td,xint)
1     xi(j)=xint
      do 2 j=m1,m2
2     xr(j)=xi(j)
      if(m1.gt.1)then
      m1=m1-1
      do 3 j=1,m1
3     xr(j)=0.
      endif
      if(m2.lt.np)then
      m2=m2+1
      do 4 j=m2,np
4     xr(j)=0.
      endif
      do 14 j=1,np
14    xi(j)=0.
      call nlogn(ip,xr,xi,-1.)
      nppp=npp+1
      do 5 i=nppp,np
      xr(i)=0.
5     xi(i)=0.
      xr(npp)=xr(npp)/2
      xi(npp)=xi(npp)/2
      xr(1)=xr(1)/2
      return
      end
c ------------------------------------------------------------
      subroutine szint (i1,i2,i3,pas,a,x,x0,aint)
      dimension a(1)
      x1=(i1-1)*pas+x0
      x2=x1+pas
      x3=x2+pas
      p1=(x-x3)*(x-x2)/((x1-x3)*(x1-x2))
      p2=(x-x1)*(x-x3)/((x2-x1)*(x2-x3))
      p3=(x-x1)*(x-x2)/((x3-x1)*(x3-x2))
      aint=a(i1)*p1+a(i2)*p2+a(i3)*p3
      return
      end
c -------------------------------------------------------------
      function nonoui(lecpic)
      character*10 yo(14),yn(8)
      character c*10
      data yo/'oui','ouai','da','ya','yes','o','y','si','yo','ye','yeeh'
     1,'ok','positif','affirmatif'/
      data yn/'non','n','no','niet','nicht','neni','nein','pas'/
C
C Cette fonction retourne 1 si reponse positive
C                         0 si reponse negative
C                        -1 si reponse incorrecte apres 3 essais
C modifiee le 25/05/93 pour lecture de la reponse dans un fichier
c si lecpic=1.

      ist=8
      nonoui=-1
      k=0
3     if (lecpic.eq.0)then
      read(*,'(a)')c
      else
      read(ist,'(a)')c
      write(*,*)'CARACTERE LU ', c
      endif
      do 1 i=1,14
      if(c.eq.yo(i))nonoui=1
1     continue
      do 2 i=1,8
      if(c.eq.yn(i))nonoui=0
2     continue
      if(nonoui.eq.-1)then
      if(k.ge.3)then
      write(*,'("advienne ce qu il adviendra")')
      return
      else
      k=k+1
      write(*,'("repondez a la question, rappel no ",i2,$)')k
      go to 3
      endif
      else
      return
      endif
      end
c -------------------------------------------------------------
      SUBROUTINE NLOGN(N,xr,xi,SIGN)
C ALGORITHME COOLEY TUKEY.
C      SIGN=-1.0  EXP (-2*I*PI*...)
C      SIGN =1.0 (1/Q) * EXP(2*I*PI*...)
C      NMAX=PLUS GRANDE VALEUR DE N AVEC
C      DIMENSION X(2**N) ,M(NMAX)
      DIMENSION XR(2),XI(2),M(20)
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
c ---------------------------------------------------------------
      subroutine des(x,y,n,xg,xd,y1,y2,ey,pas,ibase,idata,ndata,pper,
     1ppar,kmoderic,l,ipic,coef,qual,ienv,lecpic,npts,isegment)
      parameter (ndim=4096,nlob=10)
c     parameter (ndim=20000,nlob=10)
      parameter (nmode=10,nfdim=10)
      dimension xw(ndim), yw(ndim)
      dimension qual(nlob)
      dimension x(1),y(1),pper(1),ppar(1),idata(1)
c     dimension jmax(100)
      dimension ymaxi(10),imaxi(10),jjs(10),npics(10),ratiomax(10)
      dimension ratiosynt(nmode,nfdim)
      integer npts
      real ytop,ycourrant
      common/c1/ratiosynt
      common/c2/imodmax,nperpi

c       Subroutine de dessin
c       ey    = 0       ==> on calcule l'echelle a partir du min et du max
c       ibase = 1       ==> la ligne de base est en bas du dessin
c             = 0       ==>                      au y moyen
c       ipic  = 1       ==> on picke le maximum et eventuellement
c                            les points definissant le lobe
 
c ist estle numero d'unite logique ou l'on stockera les parametres des
c intercorellogrammes piques.    
      ist=8
      iredes=0
      call CloseRetainSeg(isegment)
      isegment=isegment+1
      call CreateRetainSeg(isegment)
        nsave=n
        idsave=idata(1)

      if(ienv.eq.0) then
c       dessin du cross-corr. reel
        ibase=0
        do 15 i=1,n
15      yw(i)=y(i)
      else
	if(ndata.ne.1.or.idata(1).eq.0) then
c		dessin enveloppe
        	do 14 i=1,n
14      	yw(i)=sqrt(y(i)*y(i)+x(i)*x(i))
	     else
c		dessin 1 arche pour phase
		ibase=0
      		ldt=coef*pper(l)/ppar(l)/pas+0.5
		ii=0
		idd=max(1,idata(1)-ldt)
		iff=min(nsave,idata(1)+ldt)
		do 16 i=idd,iff
        	ii=ii+1
        	yw(ii)=y(i)
16		continue
c		on change le nb de pts a dessiner
c		et la position du tiret.
        	n=ii
        	idata(1)=idata(1)-idd+1
	     endif
	endif



  100 continue

      yreduc=0.9
      dh=(y2-y1)*yreduc
      yb=y1+dh*(1.-yreduc)*0.5
      yh=y2-dh*(1.-yreduc)*0.5
      imax=1
      ymin=yw(imax)
      ymax=ymin
      do 1 i=1,n
      if(yw(i).lt.ymin) ymin=yw(i)
      if(yw(i).gt.ymax)then
      ymax=yw(i)
      imax=i
      endif
1     continue
      if(ymin.eq.0..and.ymax.eq.0.) then
      npic=0
      ndata=1
      idata(ndata)=0
      write(0,'(a,f6.1)')' Pas de cross-corr. pour ce mode a la periode'
     *   ,pper(l)
      return
      endif
      if(ienv.eq.1) ymin=0.
      if(ymax.ne.ymin) then
      	if(ey.eq.0.) ey=dh/(ymax-ymin)
      else
	write(0,*) ' ymin = ymax, echelle en Y non recalculee'
      	endif
4     ex=(xd-xg)/(n-1)
      yc=(yb+yh)/2
      ys=0.
      do 3 i=1,n
3     ys=ys+yw(i)
      ys=ys/n
      y0=yc
      if(ibase.eq.1)then
      ys=0.
      if(ey.eq.0.) ey=dh/ymax
      y0=yb
      endif
      do 2 i=1,n
      yw(i)=(yw(i)-ys)*ey+y0
2     xw(i)=(i-1)*ex+xg
      if(ibase.eq.0)then
c     	Dessin de la ligne de base
CCCC  	call MoveAbs2(xd,yc-dh/2)
CCCC  	call LineAbs2(xd,yc+dh/2)
      else
      	ybase=yb+(1-ibase)*(yc-yb)
      	call MoveAbs2(xg,ybase)
      	call LineAbs2(xd,ybase)
      	xc=xg+(n-1)*ex/2.
      	call MoveAbs2(xc,yc-dh/2)
      	call LineAbs2(xc,yc-dt)
c     	Dessin de l'echelle-temps (60 s)
      	ycc=yc+dh/4
      	top=dh/20
      	xmn=60*ex/pas+xg
      	call MoveAbs2(xg+top,ycc+top)
      	call LineAbs2(xg+top,ycc)
      	call LineAbs2(xmn+top,ycc)
      	call LineAbs2(xmn+top,ycc+top)
      endif
c     Dessin du signal
      call MoveAbs2(xw(1),yw(1))
      call PolyLineAbs2(xw,yw,n)
c
c Piquage des points pour l'intercorrelogramme a inverser
c
      if(ipic.eq.0) goto 999
      ldt=coef*pper(l)/ppar(l)/pas+0.5
      if(ibase.eq.0) ldt=pper(l)/pas/4.+0.5
      if(iredes.eq.1) goto 101

      ndata=1
      idata(ndata)=0
c
c Reperage de tous les maximums du cross-correl.
c
      call CloseRetainSeg(isegment)
c Le segment ci-dessous ne doit jamais etre sur la
c sortie papier (PostScript) ==> call minuscules
c------------ modif eric 11 Aout 1997---------------------
c modif eric 11 Aout 97:le rapport signal/bruit de chaque enveloppe
c correspond au rapport entre max de l'enveloppe et moyenne de tous
c les minimums.    
      yminiinf=9999999999.  
      yminisup=9999999999.  
      yminimoy=0.
      ysom=0.
      icompteur=0
      if(yw(1).lt.yw(2)) then
        ycourrant=((yw(1)-y0)/ey)+ys
        yminiinf=ycourrant
      endif
      do 12 i=2,n-1
      if((yw(i).lt.yw(i-1).and.yw(i).le.yw(i+1)).or.
     *   (yw(i).lt.yw(i-1).and.yw(i).le.yw(i+1))) then
          ycourrant=((yw(i)-y0)/ey)+ys
          ysom=ysom+ycourrant
          icompteur=icompteur+1
      endif
12    continue
      if(yw(n).lt.yw(n-1)) then
        ycourrant=((yw(n)-y0)/ey)+ys
        yminisup=ycourrant
      endif
      
      if(icompteur.eq.0) then
         yminimoy=(yminiinf+yminisup)/2.
      elseif (icompteur.ne.0)then
         yminimoy=ysom/icompteur
      endif

c modif eric 11 Aout 97: le maximum de l'enveloppe est
c selectionne parmis les max/(moyenne des min)>3. On pique
c celui pour lequel Amp/(dist) est maximum.

      call createretainseg(999)
      do 27 ilob=1,10
      jjs(ilob)=0
      npics(ilob)=0
      ymaxi(ilob)=0.    
      imaxi(ilob)=0
      ratiomax(ilob)=0
27    continue

      ii=0
      ytop=0
      indicsimaxselected=0
      do 11 i=2,n-1
      ratio=0
c     recherche du milieu du plateau
      if(yw(i).gt.yw(i-1).and.yw(i).ge.yw(i+1)) then
         icherch=i-1
111      icherch=icherch+1
         if(icherch.eq.n) goto 121                          
         if(yw(icherch).eq.yw(i))goto 111
c        Si ratio=max/(moyenne des min)>3 pour le fondamental
c        et 3 sur les hamrmoniques le max est retenu.
c        (On passe a 6 au deuxieme run)
c        Parmis les max retenus on selectionne les deux qui presentent
c        le meilleur rapport signal/bruit.
         ycourrant=((yw(i)-y0)/ey)+ys
         if (yminimoy.eq.0.) then
            ratio=ycourrant
         else 
            ratio=ycourrant/yminimoy
         endif
121      continue
C        Calcul du rapport max/dist
c Modif eric 8 Janvier 1998---------------------
c    Inversion des enveloppes LP: on ne selectionne les lobes que
c    pour le modes le plus excite, si le ratio est>3.
c    Inversion des enveloppes LP+CP avec ou sans phase
c    On va selectioner  un lobe pour une periode et un mode donne si:
c    -l'energie du mode considere filtre a la  periode consideree
c     contribue pour plus de 1% a l'energie du sismo symthe total
c     filtre a la meme periode
c    -le ratio est>3 (on peut clairement identifier un lobe principal)
c Modif eric 3 Avril 98
c    Dans un deuxieme temps il est prevu de refaire tourner le wfm auto
c    sur les sismogrammes rejetes au premiers run. On est alors plus 
c    severe sur la selection des donnees et on demande ratio>6 pour
c    selectioner un lobe avec 30% d'erreur et ratio>9 pour 10%
c -----------------------------------------------
      lper=l
      if(l.gt.nperpi) lper=3-(l-nperpi)
      if(npts.eq.1) then
      write(*,*)'RATIOSY,KM,Lper',ratiosynt(kmoderic,lper),kmoderic,lper
      endif
c -------------------------------------
c   A utiliser pour le premier run
        if ((pper(1).gt.50.and.ratio.ge.3.
     *       .and.ratiosynt(kmoderic,lper).gt.1..and.npts.eq.3)
     *   .or.(pper(1).le.50.and.ratio.ge.3.and.
     *        ratiosynt(kmoderic,lper).gt.1))then
c -------------------------------------
c   A utiliser pour le deuxieme run
c       if ((pper(1).gt.50.and.ratio.ge.6.
c    *       .and.ratiosynt(kmoderic,lper).gt.1..and.npts.eq.3)
c    *   .or.(pper(1).le.50.and.ratio.ge.6.and.
c    *        ratiosynt(kmoderic,lper).gt.1))then

c            write(*,*)'MODE',kmoderic,'PERIODE',pper(l)
c    *      ,'YCOUR',ycourrant,'YMINIMOY',yminimoy,'RATIO',ratio
       	ii=ii+1
            icourrant=i+int((icherch-i)/2)
            call moveabs2(xw(icourrant),yw(icourrant))
            call lineabs2(xw(icourrant),yw(icourrant)+dh/5)
            call moveabs2(xw(icourrant),ybase)
            call lineabs2(xw(icourrant),ybase-dh/5)
c           jmax(ii)=icourrant
c           write (*,*)abs(icourrant-(n/2))
            if(abs(icourrant-(n/2)).eq.0) then
c             ytop=ycourrant/(pas**0.5)
              ytop=ycourrant/(pas)
            else 
c             ytop=ycourrant/((pas*abs(icourrant-(n/2)))**0.5)
              ytop=ycourrant/(pas*abs(icourrant-(n/2)))
            endif

c       selection de 1 ou 2 lobes presentant un ytop max
c       s'ils ne sont pas en bordure du cross-correlograme
c       Marche si on selectionne 1 ou 3 points
        if(npts.gt.3)
     *   stop 'Modifier cross12 pour selectioner
     *         plus de 3 points'
c         write(*,*)'ICOURR',icourrant,ldt
          if (((icourrant-ldt).gt.0).and.((icourrant+ldt).lt.nsave))then
            if(npts.eq.3.and.ytop.ge.ymaxi(1))then
               if(indicsimaxselected.eq.1)indicsimaxselected=2
               if(indicsimaxselected.eq.0)indicsimaxselected=1
               ratiomax(2)=ratiomax(1)
               ymaxi(2)=ymaxi(1)
               imaxi(2)=imaxi(1) 
               jjs(2)=jjs(1)
               npics(2)=npics(1)
               ratiomax(1)=ratio
               ymaxi(1)=ytop
               imaxi(1)=icourrant
               jjs(1)=ii
               npics(1)=npts
            else if (npts.eq.3.and.ytop.gt.ymaxi(2)) then
               indicsimaxselected=2
               ratiomax(2)=ratio
               ymaxi(2)=ytop
               imaxi(2)=icourrant
               jjs(2)=ii
               npics(2)=npts
c           else if(npts.eq.1.and.ytop.gt.ymaxi(1)) then
c      ratio>6 (ou 9 suivant run) car on veut etre sur de ne fitter la phase
c      que lorsqu'on a mis 10% au max d'erreur sur l'enveloppe.
c      A utiliser pour premier run
            else if(npts.eq.1.and.ratio.gt.6.and.ytop.gt.ymaxi(1)) then
c      A utiliser pour deuxieme run
c           else if(npts.eq.1.and.ratio.gt.9.and.ytop.gt.ymaxi(1)) then
c             write(*,*)'BONJNPTS=1,indice',ii
              indicsimaxselected=1
              ratiomax(1)=ratio 
              ymaxi(1)=ytop      
              imaxi(1)=icourrant
              jjs(1)=ii
              npics(1)=npts
            endif
          endif
         endif
      endif
11    continue

c     Modif Vogle test eric
      call CloseRetainSeg(999)

c On classe les 2 lobes par indice croissant
      if((ymaxi(1).ne.0.and.ymaxi(2).ne.0).and.
     *  (jjs(2).lt.jjs(1))) then
         bidratio=ratiomax(2) 
         bidymax=ymaxi(2) 
         bidimax=imaxi(2) 
         bidjj=jjs(2)    
         bidnp=npics(2)
         ratiomax(2)=ratiomax(1)
         ymaxi(2)=ymaxi(1)
         imaxi(2)=imaxi(1)
         jjs(2)=jjs(1)    
         npics(2)=npics(1)
         ratiomax(1)=bidratio    
         ymaxi(1)=bidymax 
         imaxi(1)=bidimax  
         jjs(1)=bidjj     
         npics(1)=bidnp   
      endif 

  5    continue
c 5    write(0,'(a)') " piquage ,periode",l
c      write(*,*)'MODE,PERIODE,RATIOS',kmoderic,pper(l)
c    *,(ratiomax(ij),ij=1,indicsimaxselected)

       if(indicsimaxselected.eq.0) then
        jj=0
        npic=0
c NB eric 4/06/2002 necessiterai une declaration qual(0:nlob)
c qual(0) n'est pas utilise ulterieurement
c        qual(0)=0.
        if(lecpic.ne.1)write(ist,*)0,0,0.
       endif

      do 28 ilob=1,indicsimaxselected
       if(npts.eq.3) then
        jj=jjs(ilob)
        npic=npics(ilob)
c   A utiliser pour premier run
        if (ratiomax(ilob).ge.6.) then
c   A utiliser pour deuxieme run
c       if (ratiomax(ilob).ge.9.) then
          qual(ilob)= 10.
        else
          qual(ilob)=30.
        endif
       elseif(npts.eq.1) then
        jj=jjs(1)
        npic=npics(1)
        qual(ilob)=5.
       endif
       if(jj.gt.ii.or.
     * (npic.ne.0.and.npic.ne.1.and.npic.ne.3.and.npic.ne.5)) then
 	write(0,*) " Valeur interdite "
 	go to 5
       endif
  28    continue

      if(indicsimaxselected.ne.0) then
      if(lecpic.ne.1)write(ist,*)
     *(jjs(ilob),ilob=1,indicsimaxselected),npts,
     *(qual(ilob),ilob=1,indicsimaxselected)
      endif

c
c       Effacement des reperes des maximas
c
        call delretainsegment(999)
        isegment=isegment+1
        call CreateRetainSeg(isegment)
c
      do 29 ilob=1,indicsimaxselected
         if(npts.eq.3) then
          jj=jjs(ilob)
          npic=npics(ilob)
c   A utiliser pour premier run
        if (ratiomax(ilob).ge.6.) then
c   A utiliser pour deuxieme run
c       if (ratiomax(ilob).ge.9.) then
            qual(ilob)=10.
          else
            qual(ilob)=30.
          endif
         elseif(npts.eq.1) then
          jj=jjs(1)
          npic=npics(1)
          qual(ilob)=5.
         endif
c
        if(jj.le.0.or.npic.le.0) return
c       if(npic.eq.1)write(0,*)'==> pickage phase'
c       imax=jmax(jj)
        imax=imaxi(ilob)
        qual(ilob)=qual(ilob)*0.01
        npic=npic/2+1
        if(npic.ne.1)ldt=ldt/(npic-1)
c       Reperage des indices pour tous les points selectiones
c       Marche si 3pts selectiones par lobe. Pas teste pour
c       5 pts par lobe.
        if(ilob.eq.1) then
          do 6 i=1,npic
          ii=i+npic-1
          idata(ii)=imax+ldt*(i-1)
          iii=npic-i+1
6         idata(iii)=imax-ldt*(i-1)
        endif
        if(ilob.eq.2) then
          do 66 i=1,npic
          ii=i+4
          idata(ii)=imax+ldt*(i-1)
          iii=6-i
66        idata(iii)=imax-ldt*(i-1)
        endif
        if((indicsimaxselected.gt.2).or.(npts.gt.3))
     *   stop 'Modifier cross12 pour piquer plus de 
     *        lobes ou selectioner plus de 3 points'
29    continue

        if(npts.eq.1) then
           ndata=1
        elseif (npts.eq.3.and.indicsimaxselected.ne.0) then
           ndata=npts*indicsimaxselected
        endif
        if(indicsimaxselected.eq.0) then
           ndata=1
           idata(1)=0
        endif

       write(*,*)'IMAXI,N',(imaxi(i),i=1,indicsimaxselected)
       write(*,*)'IDATA,N',(idata(i),i=1,ndata),n
      do 8 i=1,ndata
88      continue
        if(idata(i).lt.0)then
             idata(i)=idata(i)+nsave
             stop' PAS NORMAL:les donnees en bodures 
     *             des cross-corr doivent normalement etre rejetes'
        endif
        if(idata(i).gt.nsave) then
             idata(i)=idata(i)-nsave
             stop' PAS NORMAL:les donnees en bodures  
     *             des cross-corr doivent normalement etre rejetes'
        endif
        if(idata(i).lt.0.or.idata(i).gt.nsave)go to 88
    8 continue
      idsave=idata(1)

  101 continue

c	dessin des tirets de pickage
      do 7 i=1,ndata
      call MoveAbs2(xw(idata(i)),yw(idata(i)))
      call LineAbs2(xw(idata(i)),yw(idata(i))-dh/3)
   7  continue
c
      if(ienv.eq.1.and.npic.eq.1) then
c
c	Effacement du dessin de l'enveloppe, des tirets de pickage,
c	 et remplacement par le dessin du cross-correl. reel
c
	iredes=1
	call CloseRetainSeg(isegment)
      	call DelRetainSegment(isegment-1)
      	call DelRetainSegment(isegment)
      	isegment=isegment+1
      	call CreateRetainSeg(isegment)
      	npic=0
	ii=0
	idd=max(1,imax-ldt)
	iff=min(nsave,imax+ldt)
	do 312 i=idd,iff
	ii=ii+1
	yw(ii)=y(i)
  312	continue
	n=ii
	idsave=idata(1)
	idata(1)=idata(1)-idd+1
      	ey=ey/2.
	ibase=0
      	goto 100
      	endif
999	continue

	n=nsave
	idata(1)=idsave
      if(iredes.eq.1) ey=ey*2.

      return
      end
c --------------------------------------------------------------
      subroutine szopen(l,c)
      character*160c
c
c Ce sous-programme ouvre l unite logique "l" en lecture
c apres reponse a la question posee dans la chaine de caracteres "c"
c 
c     character*30f
      character*100f
1     write(0,c)
      read(*,'(a)')f
      open(l,status='old',err=2,file=f)
      return
2     write(0,'("fichier inexistant")')
      go to 1
      end
c ---------------------------------------------------------------
      SUBROUTINE SZPU2 (N,NP,I)
      I=1
    1 NP=2**I
      IF(NP-N) 2,3,3
    2 I=I+1
      GO TO 1
    3 RETURN
      END
c------------------------------------------------------------------
      subroutine sztrans(w,y,n)
      dimension w(1),y(1)
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
c-------------------------------------------------------------------
      subroutine fildes(xr,xi,npp,dnu)
      dimension xr(1),xi(1)
c
c          Filtrage passe-bande d'un sismogramme sous forme TF
c
c                 npp  indice de la frequence de Nyquist du signal
c                 xr   partie reelle de la TF
c                 xi   partie imaginaire de la TF
c                 dnu  pas en frequence de la TF
c
c
c Filtrage passe bande
c
      per=1./dnu/(npp-1)
15    write(*,'("la periode de Nyquist doit etre superieure a",
     *      f10.4)')per
16    write(*,'("entrer les periodes extremes du plateau du filtre: ",$)
     *')
      read(*,*)p2,p3
      if(p3.gt.p2) then
      pbid=p2
      p2=p3
      p3=pbid
      endif
      write(*,'("entrer les periodes ou le filtre revient a zero: ",$)
     *')
      read(*,*)p1,p4
      if(p4.gt.p1) then
      pbid=p1
      p1=p4
      p4=pbid
      endif
      if (p4.lt.per) goto 15
      if(p2.gt.p1.or.p4.gt.p3) then
      write(*,*) 'Les periodes definissant le filtre sont incoherentes'
      goto 16
      endif
      if1=1./p1/dnu+1.5
      if2=1./p2/dnu+1.5
      if3=1./p3/dnu+1.5
      if4=1./p4/dnu+1.5
        if(if1.le.1) if1=2
        if(if2.le.1) if2=2
        if(if3.gt.npp) if3=npp
        if(if4.gt.npp) if4=npp
      f1=dnu*(if1-1)
      f2=dnu*(if2-1)
      f3=dnu*(if3-1)
      f4=dnu*(if4-1)
      write(*,*)'frequences reelles du filtre:',f1,f2,f3,f4
      write(*,*)'periodes reelles du filtre:', 1./f4,1./f3,1./f2,1./f1
      id2=if2-if1
      id3=if4-if3
      if(id2.lt.1) goto 88
      do 7 i=if1,if2
      u=float(i-if2)/float(id2)
      u=u*u
      u=(1.-u)*(1.-u)
      xr(i)=xr(i)*u
7     xi(i)=xi(i)*u
88    if(id3.lt.1) goto 19
      do 71 i=if3,if4
      u=float(i-if3)/float(id3)
      u=u*u
      u=(1.-u)*(1.-u)
      xr(i)=xr(i)*u
71    xi(i)=xi(i)*u
19    continue
c--- Modif eric 4/06/2002 compilation Linux
c     if(if1.le.1) goto 17
      if(if1.le.1) goto 1711
c------------------------------------------
      do 17 i=1,if1-1
      xr(i)=0.
      xi(i)=0.
17    continue
1711  continue
c--- Modif eric 4/06/2002 compilation Linux
c     if(if4.ge.npp) goto 171
      if(if4.ge.npp) goto 1710
c------------------------------------------
      do 171 i=if4+1,npp
      xr(i)=0.
      xi(i)=0.
171   continue
1710   continue

      return
      end
c ---------------------------------------------------------------
      subroutine desu(n,xg,xd,y0,yt,dist,t1,dt,dh,entr,
     *                  titre,nper,pper,km,qualit,ey,amult)
      parameter (nmode=10,nfdim=10)


c     Dessin d'un axe gradue en vitesse de groupe, avec
c     un pas de 0.5 km/s.

      dimension pper(1),qualit(nmode,nfdim)
      dimension ey(nmode,nfdim)
      character texte*20, entr*(*), titre*80, forma*80

      call SetCharPrecision(CHARACTER)
      write(*,*) "called SetCharPrecision"

      call SetFont(STICK)
      call SetCharSize((xd-xg)*.015,dh*.15)

      call MoveAbs2(xg,yt+8.5)
      call Text(titre)

      call SetCharSize((xd-xg)*.01,dh*.10)

      call MoveAbs2(xd-0.5*dh,y0+2.6*dh)
      call Text('resid.')
c     Dessin du nom du fichier synthetique pour l'iteration en cours
      call MoveAbs2(xd-0.5*dh,y0+1.6*dh)
      call Text(entr)
      write(titre,'(a,f5.2)') 'x',amult
      call MoveAbs2(xd-0.5*dh,y0+1.3*dh)
      call Text(titre)
      call MoveAbs2(xd-0.5*dh,y0+0.6*dh)
      call Text('data')

c     Dessin des graduations en vitesse de groupe
      kmax=int(2*dist/t1)
      kmin=int(2*dist/(t1+(n-1)*dt)-0.001)+1

      ex=(xd-xg)/(n-1)

      write(texte,'(f3.1)') kmax/2.
      call MoveAbs2(xg+(xd-xg)*.005,y0+0.10*dh)
      call Text(texte)

      call MoveAbs2(xg,y0)
      write(*,*) "called moveAbs"

      do 1 k=kmax,kmin,-1
      u=k/2.
      i=(dist/u-t1)/dt+1
      y=y0+0.07*dh*(2-(k-2*(k/2)))
      xx=(i-1)*ex+xg
	xx=xg+(dist/u-t1)*(xd-xg)/((n-1)*dt)
      call LineAbs2(xx,y0)
      call LineAbs2(xx,y)
      call LineAbs2(xx,y0)
    1 continue

      call LineAbs2(xd,y0)

      if((xd-xx).lt.(xd-xg)*0.075) xx=xd-(xd-xg)*0.080
      write(texte,'(f3.1,a)') kmin/2., ' km/s'
      call MoveAbs2(xx+(xd-xg)*.005,y0+0.10*dh)
      call Text(texte)

c     Dessin de l'echelle-temps (60 s)
       ycc=y0+dh*0.75
       top=dh/20
       xmn=xg+(xd-xg)*60/((n-1)*dt)
       call MoveAbs2(xg+top,ycc+top)
       call LineAbs2(xg+top,ycc)
       call LineAbs2(xmn+top,ycc)
       call LineAbs2(xmn+top,ycc+top)
      call MoveAbs2(xg+(xd-xg)*.02,y0+0.83*dh)
      call Text('60 s')

c     Trace des periodes des cross-corr.
      if(pper(1).eq.0) return
      xcal=xg+(xd-xg)/(3*nper)
      xdec=(xd-xg)/nper
      do 2 i=1,nper
      call SetCharSize((xd-xg)*.01,dh*.10)
      call MoveAbs2(xcal+(i-1)*xdec,yt+2.)
      write(titre,'(a,f5.1,a)')'T =',pper(i),' s'
      call Text(titre)
      iq=0
      iey=0
      do 22 j=1,km
c  	iq=0 ==> qualites identiques
      if(abs(qualit(j,i)-qualit(1,i)).gt.1.e-4) iq=1
      if(ey(j,i).ne.ey(1,i)) iey=1
   22 continue

      call SetCharSize((xd-xg)*.006,dh*.06)
      write(*,*) "called SetCharSize"

      write(*,*) iq, abs(qualit(1,i))
      if(iq.eq.0.and.abs(qualit(1,i)).lt.0.0001) goto 299
c       On n'ecrit ce titre que si on a picke des donnees sur les
c       cross-correlogrammes.
	ncar=3
       write(*,*) "ncar is", ncar
	do 123 j=1,km
	if(qualit(j,i).ne.0.) nc=log10(int(100*qualit(j,i))+0.000001)+1
	if(nc.gt.ncar) ncar=nc
  123	continue

c        write(*,*) "hello"

c PROBLEM---------------------- \/
      ncar=3
      write(forma,'(a,i5,a)')'("  Cd = ",10(f',ncar+2,'.1,"% "))'
       write(*,*) forma, ncar
      if(iq.eq.0) then
c	Cas ou les incertitudes relatives sont identiques pour toute
c	la colonne.
      	write(titre,forma) 100.*qualit(1,i)
      else
c        write(*,*) "hello1"
c	cas ou les incertitudes sont differentes d'un mode a l'autre.
      	write(titre,forma) (100.*qualit(j,i),j=1,km)
c        write(*,*) "hello2"
      endif
     
c PROBLEM----------------------  /\
      write(*,*) "finished if statements"

      call MoveAbs2(xg+0.2*(xd-xg)/nper+(i-1)*xdec,yt-0.5)
      call Text(titre)
  299 continue

c      write(*,*) "hello", iey
      if(iey.eq.0) then
c       On n'ecrit ce titre que si l'echelle est la meme pour toute
c       la colonne.
      call MoveAbs2(xg+0.4*(xd-xg)/nper+(i-1)*xdec,y0+(yt-y0)*0.485)
      write(titre,'(a,f7.4)') 'Y-SCALE :',ey(1,1)/ey(1,i)
      call Text(titre)
      endif
c      write(*,*) "end if "

c       On n'ecrit ce titre que si on a picke des donnees sur les
    2 continue

c     Trace des numeros de mode a gauche des cross-corr.
      ycal=y0+(yt-y0)/2*(1+0.2/km)
      ydec=(yt-y0)/(2*km)
      call SetCharSize((xd-xg)*.008,dh*.08)
      call MoveAbs2(xg,ycal+(km-0.5)*ydec)
      call Text('mode')
      call SetCharSize((xd-xg)*.012,dh*.12)
      do 3 i=1,km
      call MoveAbs2(xg,ycal+(i-1)*ydec)
      write(titre,'(i2.1)') i-1
      call Text(titre)
    3 continue

c     Trace des indicateur de cross-corr.
      xcal=(xd-xg)/nper*1.03
      xdec=(xd-xg)/nper
      ycal=y0+(yt-y0)/2
      ydec=(yt-y0)/(2*km)
      ycal=ycal+0.1*ydec
      yyd=ydec*0.5*0.8
      call SetCharSize((xd-xg)*.008,dh*.08)
      do 4 i=1,nper-1
      do 4 j=1,km
      call MoveAbs2(xcal+(i-1)*xdec,ycal+(j-1)*ydec)
      call Text('data')
      call MoveAbs2(xcal+(i-1)*xdec,ycal+(j-1)*ydec+yyd)
      call Text('synt')
    4 continue

      return
      write(*,*) "finished subroutine"
      end
c ------------------------------------------------------------
      function mouinon()
c     fonction ne lisant pas l'entree standard mais le clavier
      character*30 ttynam, name
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
      mouinon=-1
      k=0
      name=ttynam(6)
      if(name.eq.'')then
                   mouinon=0
                   return
                   endif
      open(1,file=name)
3     read(1,'(a)')c
      do 1 i=1,14
      if(c.eq.yo(i))mouinon=1
1     continue
      do 2 i=1,8
      if(c.eq.yn(i))mouinon=0
2     continue
      if(mouinon.eq.-1)then
         if(k.ge.3)then
            write(*,'("advienne ce qu il adviendra")')
            close(1)
            return
            else
            k=k+1
            write(*,'("repondez a la question, rappel no ",i2,$)')k
            go to 3
            endif
      else
          close(1)
            return
      endif
      end

