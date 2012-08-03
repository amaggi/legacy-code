c---------------------------------------------------------------
C       5/06/93       @(#)lecdgk12.phase.lr.f	1.4
* sous-programme de lecture des fichiers modele et donnees adapte
* a l'inversion de plusieurs trajets.
c--modif le 30 avril 1997 permettant l'inversion de lnMo
c  on a dg/dlnM = M*dg/dM  = gp = Dgr dans ce pgm. JJ a verifie
c  qu'on a bien Dgr = gp.
c  Il suffit donc de ne pas diviser Dgr par M comme cela etait
c  fait dans les versions precedentes. 
c le passage a l'inversion pour log10Mo permet d'eviter le
c scaling ene 25 qui de P0, Dp0 Pinit( plus de risque d'overflow)
c On supprime ces e25 partout
c---------------------------------------------------------------

      Subroutine LECDGK(P0,CP0,DP0,D0,CD0,Pinit,Dinit,G,n,mnbs,nd,
     *		mp,nvar,Prof,titre,distcor,Coupl,ncouch,nmode,
     *		nfiltr,ienv,Ndata,Sig2r,Sig2i,nbs,nbt,nt)

	parameter(echpha=1.e4)
      Parameter (mmp=200, nnd=200)
      Dimension P0(mp), CP0(mp,mp), DP0(mp), D0(nd), CD0(nd,nd),
     *    Pinit(mp), Dinit(nd), G(nd,mp), Coupl(7,7), Idata(nnd),
     *    nfiltr(*) 
      Dimension  Sig2r(nnd),Sig2i(nnd)
      Dimension  Sigr(nnd),Sigi(nnd), Dgr(nnd),Dgi(nnd), Prof(mmp),
     *               Ndata(nt,nd,nd), Ep(mmp), Xinf(mmp), Xsup(mmp)
      character name*40, text*30, titre*80
      real*8 pi, pim
      real tfiltr,par

      data pi/3.14159265359/

      if(mmp.ne.mp.or.nnd.ne.nd) stop 'erreur de dimension dans LECDGK'

      il1=9
      write(*,*) 'entrez le nom du fichier-modele(tv.n)'
      read(*,'(a)') name
      open(il1, file=name, status='old', err=999)

      write(*,*) 'entrez le nom du fichier-listing (tv.n+1)'
      read(*,'(a)') name
      open(10, file=name)

      read(il1,'(a)') titre
      read(il1,*) nvar, ncouch
      m=nvar*ncouch

c   lecture du modele a priori
c    de maniere coherente avec le remplissage du tableau G
      do 10 i=1,ncouch
10    read(il1,*) Prof(i),(P0(k),k=i,nvar*ncouch,ncouch)
c-----------------------------------------------------------------
c   attention, pour les fichiers en entree de ce programme
c   on devra mettre en premier les parametres correspondant 
c   aux seismes inverses pour love et rayleigh
      write(*,*) 'entrez le nombre de seismes a inverser:'
      read (*,*) nbs
      write(*,*) 'entrez le nombre de seismes inverses pour love'
      write(*,*) 'et rayleigh:'
      read (*,*) nlr
c   nbt est le nombre de trajets  total a inverser
      nbt=nlr+nbs
c----------------------------------------------------------------
c   le fichier tv.n devra contenir autant de moment sismique qu'il
c   y a de seismes.

c   lecture du moment sismique et de son incertitude
      mnbs=m+nbs
	if(mnbs.gt.mp) stop 'dimension mp trop petite'
      do 12 i=m+1,mnbs
12    read(il1,*) P0(i)
c----mise a l'echelle numerique de P0,DP0,Pinit et G (facteur 10**25)
c--12 P0(i)=P0(i)*1.e-25

c   lecture des ecarts-type a priori pour ce modele
      do 20 i=1,ncouch
20    read(il1,*) Prof(i),(DP0(k),k=i,nvar*ncouch,ncouch)
      do 22 i=m+1,mnbs
22    read(il1,*) DP0(i)
c--22-DP0(i)=DP0(i)*1.e-25

c   lecture du modele resultant de l'iteration precedente
      do 30 i=1,ncouch
30    read(il1,*) Prof(i),(Pinit(k),k=i,nvar*ncouch,ncouch)
      do 32 i=m+1,mnbs
32    read(il1,*) Pinit(i)
c--32-Pinit(i)=Pinit(i)*1.e-25

c   calcul de la matrice de covariance a priori
      read(il1,*) distcor
      do 40 i=1,nvar
40    read(il1,*) (Coupl(i,j), j=1,nvar)
      call epprof(Ep,Prof,Xinf,Xsup,nvar,ncouch)
      call makecp0(CP0,DP0,Coupl,Ep,Xinf,Xsup,distcor,
     *              nvar,ncouch,mp)
      ij=1
      do 42 i=m+1,mnbs
      CP0(i,i)=DP0(i)*DP0(i)
c     write(*,*) (i,ij,CP0(i,ij),ij=1,i)
42    continue

c   Lecture des donnees reelles.
c   Ici ienv1 et nmode devront toujours avoir la meme valeur pour
c   les differents trajets.Ceci n'est pas le cas du nombre de periodes
c   qui peut changer entre chaque trajets.
      l=0
      ll=0
      do 82 inbt=1,nbt
      read(il1,*) nmode,nfiltr(inbt),ienv1
c     write(*,*) nmode,nfiltr(inbt),ienv1
      do 90 imode=1,nmode
      do 90 ifiltr=1,nfiltr(inbt)
      read(il1,*) ndat,(Sigr(i),Sigi(i),i=1,ndat)
      Ndata(inbt,imode,ifiltr)=ndat
      do 100 i=1,ndat
      l=l+1
      D0(l)=Sigr(i)
      if(ienv1.eq.1.and.ndat.ne.1)
     *    D0(l)=sqrt(Sigr(i)*Sigr(i)+Sigi(i)*Sigi(i))
      if(ienv1.eq.1.and.ndat.eq.1.and.(Sigi(i).ne.0..or.Sigr(i).ne.0.))
     *		D0(l)=echpha*atan2(Sigi(i),Sigr(i))
100    continue
90     continue
       n=l

c   lecture des erreurs sur les donnees reelles
      do 130 imode=1,nmode
      do 130 ifiltr=1,nfiltr(inbt)
      read(il1,*) ndat,(Sigr(i),i=1,ndat)
      do 140 i=1,ndat
      ll=ll+1
140   CD0(ll,ll)=Sigr(i)*Sigr(i)
130   continue
      if(ll.ne.l)stop'incompatibilite dans CD0'
82     continue

      close(il1)

      il2=9
      write(*,*) 'entrez le nom du fichier-(dgk.n)'

      read(*,'(a)') name
      write(*,*) "name is ", name

      open(il2, file=name, status='old', err=999)

c   lecture des derivees partielles dgk
      
c       k= indice du parametre dans l'inversion
c       n= indice de la donnee dans l'inversion
c       n1tot= nb total de donnees pour le mode en cours
c       ntot = nb total de donnees pour tous les modes

c   remplissage de la matrice des derivees partielles dans le cas 
c   de nbt trajets donc boucle surle nombre de trajets


      ntot=0
      jsei=0
      nlr2 = nlr*2
      l = 0
      ll = 0

      do 58 inbt=1,nbt

      read(il2,*) text,nmod,nfilt,nva,ncouc,ienv
      write(*,*) "*******",text,nmod,nfilt,nva,ncouc,ienv

c      write(*,*) "hello"

      if(nmod.ne.nmode. or .nfilt.ne.nfiltr(inbt). or .nva.ne.nvar. or .
     * ncouc.ne.ncouch. or .ienv.ne.ienv1) 
     * stop ' Fichiers modele et donnees incompatibles'
      if (inbt.gt.nlr2.or.mod(inbt,2).ne.0) jsei=jsei+1
      m=nvar*ncouch
      do 60 imode=1,nmode
      k=0
      inew=1
      do 70 ivar=1,nvar
      do 70 icouch=1,ncouch
      k=k+1
      n1tot=0
      do 990 ifiltr=1,nfiltr(inbt)
c      write(*,*) "hello1"
c PROBLEM ------------------------- \/ ------------------------
        write(*,*) inew
	if(inew.eq.1) then
c--modif du 30 avril 1997 visant a  permettre l'inversion de
c--lnQ. On ne divise plus Dgr/Pinit.

c	lecture des der. part. du moment sismique

        write(*,*) "ndat is " ,ndat, " il2 is ", il2

c      	read(il2,*) jmode
c      	write(*,*) "jmode is ", jmode
      	read(il2,*) jmode,tfiltr,par,jvar,jcouch,
     *    ndat,(idata(i),i=1,ndat), (Dgr(i),Dgi(i),i=1,ndat)

        write(*,*) "ndat is ", ndat

      	write(*,*) jmode,tfiltr,par,jvar,jcouch,
     *    ndat,(idata(i),i=1,ndat), (Dgr(i),Dgi(i),i=1,ndat)

c PROBLEM -------------------------   /\ -------------------------
      write(*,*)" 29 DP M0 LUES"
c     il y aura autant de moments sismiques que de seismes (nbs)
      	do 880 km0=m+1,m+nbs
      	do 880 i=1,ndat
      	n=ntot+n1tot+i
            if (km0.eq.(m+jsei)) then 
c      	g(n,km0)=Dgr(i)/Pinit(km0)
       	g(n,km0)=Dgr(i)*alog(10.)
c           write(*,*)'dp Mo : ',g(n,km0)
            else 
       	g(n,km0)=0.
            endif
c     write(*,*) n,km0,g(n,km0),g(n,km0)
880   continue    
      	endif

c      write(*,*) "hello"

      read(il2,*) jmode,tfiltr,par,jvar,jcouch,
     *            ndat,(idata(i),i=1,ndat), (Dgr(i),Dgi(i),i=1,ndat)
      write(*,*)'29 DP SU LUES m,t,p,v,c',jmode,tfiltr,par,jvar,jcouch
      do 80 i=1,ndat
      n=ntot+n1tot+i
      g(n,k)=Dgr(i)
c     write(*,*) n,k,g(n,k)
80    continue
      n1tot=n1tot+ndat
990   continue
      inew=0
70    continue
	if(km0.ne.k+nbs+1) stop ' Pb avec der. part. du moment sismique'
      ntot=ntot+n1tot
60    continue
      if(inbt.eq.nbt.and.ntot.ne.n) stop 'incompatibilite 
     *dans les dgk'

C     DO 227 nn=1,ntot
C     WRITE(*,'(11e11.4)') (g(nn,j),j=1,k)
C 227 continue

c      write(*,*) "hello"

c     on relit dans dgk.n les donnees reelles, deja lues dans tv.n
      do 230 imode=1,nmode
      do 230 ifiltr=1,nfiltr(inbt) 
      read(il2,*) ndat,(Sigr(i),Sigi(i),i=1,ndat)
      write(*,*)'29 DON REELLES LUES'
      do 240 i=1,ndat
      l=l+1
      Sig2r(l)=Sigr(i)
      Sig2i(l)=Sigi(i)
      Dinit(l)=Sigr(i)
      if(ienv.eq.1.and.ndat.ne.1)
     *    Dinit(l)=sqrt(Sigr(i)*Sigr(i)+Sigi(i)*Sigi(i))
      if(ienv1.eq.1.and.ndat.eq.1.and.(Sigi(i).ne.0..or.Sigr(i).ne.0.))
     *		Dinit(l)=echpha*atan2(Sigi(i),Sigr(i))
c     Test de compatibilite des donnees
      if(ndat.eq.1.and.D0(l).eq.0..and.Dinit(l).eq.0.) goto 240
      test=abs((Dinit(l)-D0(l))/D0(l))
      if(test.gt.1.e-4) then
               write(*,'(10e14.6)')(Dinit(j),j=1,l)
               write(*,'(10e14.6)')(D0(j),j=1,l)
               write(*,*)'test',test  
               stop 'Signal de reference incorrect ?'
               endif
240   continue
230   continue
c     lecture des donnees synthetiques
      do 110 imode=1,nmode
      do 110 ifiltr=1,nfiltr(inbt) 
      read(il2,*) ndat,(Sigr(i),Sigi(i),i=1,ndat)
      write(*,*)'29 DON SYNTHETIQUES LUES'
      do 120 i=1,ndat
      ll=ll+1
      Dinit(ll)=Sigr(i)

c      write(*,*) "hello"

      if(ienv.eq.1.and.ndat.ne.1)
     *    Dinit(ll)=sqrt(Sigr(i)*Sigr(i)+Sigi(i)*Sigi(i))
      if(ienv1.eq.1.and.ndat.eq.1.and.(Sigi(i).ne.0..or.Sigr(i).ne.0.))
     *		then
		pim=pi*echpha
		Dinit(ll)=echpha*atan2(Sigi(i),Sigr(i))
		if(Dinit(ll)-D0(ll).gt.pim) Dinit(ll)=Dinit(ll)-2*pim
		if(D0(ll)-Dinit(ll).gt.pim) Dinit(ll)=Dinit(ll)+2*pim
		dphas=(D0(ll)-Dinit(ll))*180/pim
		if(abs(dphas).gt.90.)write(*,*)'Achtung: le',
     *		' dephasage peut etre incorrect de 2.k.pi. Verifier.'
		write(*,*) 'Dephasage donnee',l,': ', dphas, ' degres'
		endif
120    continue
110    continue

58     continue
      if(ll.ne.l.and.l.ne.n)stop'incompatibilite dans sigsynt 
     *du fichier donnees'
      write(*,*) n,mnbs,nvar,ncouch,nbs

      close(il2)
      return
  999 stop 'erreur de lecture fichier'
      write(*,*) "end sub"
      end
c----------------------------------------------------------------------
      subroutine epprof(ep,prof,xinf,xsup,nvar,ncouch)
*     passage des profondeurs aux epaisseurs
      parameter (mmp=200)

      dimension ep(mmp),prof(mmp),xinf(mmp),xsup(mmp)

      ep(1)=prof(2)-prof(1)
      xinf(1)=prof(1)-ep(1)/2.
      xsup(1)=(prof(1)+prof(2))/2.

      do 1 i=2,ncouch-1
      ep(i)=(prof(i+1)-prof(i-1))/2.
      xinf(i)=(prof(i-1)+prof(i))/2.
      xsup(i)=(prof(i)+prof(i+1))/2.
    1 continue

      ep(ncouch)=prof(ncouch)-prof(ncouch-1)
      xinf(ncouch)=(prof(ncouch-1)+prof(ncouch))/2.
      xsup(ncouch)=prof(ncouch)+ep(ncouch)/2.

      do 2 i=ncouch+1,ncouch*nvar
      ep(i)=ep(i-ncouch)
      prof(i)=prof(i-ncouch)
      xinf(i)=xinf(i-ncouch)
      xsup(i)=xsup(i-ncouch)
    2 continue

***   do 3 i=1,ncouch*nvar
***   write(8,'(10f8.2)') ep(i),xinf(i),prof(i),xsup(i)
*** 3 continue

      return
      end
c---------------------------------------------------------------------
      Subroutine makecp0(CP0,SP0,coupl,ep,xinf,xsup,alcorr,
     *                   nvar,ncouch,mp)
      parameter(mmp=200)

c     Ce programme calcule la matrice de covariance
c     pour le modele discretise
c     a partir de la covariance continue definie par des
c     ecarts-type continus et une correlation gaussienne
c     de largeur  alcorr .

      dimension CP0(mp,mp),SP0(mp),coupl(7,7),ep(mp)
      dimension xinf(mp),xsup(mp)

         iwarn=0

         m=nvar*ncouch

         do 1 i=1,m

         if(alcorr.lt.ep(i)) iwarn=1

	 ivar=(i-1)/ncouch+1

         do 2 j=1,i
	 jvar=(j-1)/ncouch+1
         coeff=gauint(xinf,xsup,i,j,alcorr,mp,m)
         coeff=coeff*coupl(ivar,jvar)
         CP0(i,j)=SP0(i)*SP0(j)*coeff/(ep(i)*ep(j))
         CP0(j,i)=CP0(i,j)
    2    continue

    1 continue

         if(iwarn.eq.1) write(*,*) 'Attention, certaines couches sont 
     *plus epaisses que la longueur de correlation: on est loin des hypo
     *theses de discretisation d''un modele continu'

      return
      end
c---------------------------------------------------------------------
      real*4 function gauint(dxinf,dxsup,i,j,dalcorr,mp,m)
      parameter(mmp=200)

c     calcul de l'integrale double pour zi1<zi<zi2 et zj1<zj<zj2
c     de la gaussienne exp( -(zj-zi)**2/(2*alcorr**2))

      implicit real*8 (a-h,o-z)
c----modif eric 05/02/2002 pour compilation Linux------------
c     real*4 dxinf,dxsup,dalcorr,gauint
      real*4 dxinf,dxsup,dalcorr
c------------------------------------------------------------
      dimension dxinf(mp), dxsup(mp)
      alcorr=dalcorr
      pi=3.14159265359
      sq2=dsqrt(2.d0)
      zi1=dxinf(i)
      zi2=dxsup(i)
      zj1=dxinf(j)
      zj2=dxsup(j)

c     int(u1,u2)(exp(-u**2/(2*L**2)))
c             = L*sqrt(pi/2)*(erf(u2/L/sq(2))-erf(u1/L/sq(2)))

c     int(z,infini)(erfc(t)) = -z*erfc(z) + exp(-z**2)/sqrt(pi)

      a11=(zj1-zi1)/(alcorr*sq2)
      a21=(zj2-zi1)/(alcorr*sq2)
      a12=(zj1-zi2)/(alcorr*sq2)
      a22=(zj2-zi2)/(alcorr*sq2)
      g11=(dsqrt(pi)*a11*derfc(a11)-dexp(-a11*a11))
      g12=(dsqrt(pi)*a12*derfc(a12)-dexp(-a12*a12))
      g21=(dsqrt(pi)*a21*derfc(a21)-dexp(-a21*a21))
      g22=(dsqrt(pi)*a22*derfc(a22)-dexp(-a22*a22))
      gauint=alcorr*alcorr*(g11-g12-g21+g22)


      return
      end
