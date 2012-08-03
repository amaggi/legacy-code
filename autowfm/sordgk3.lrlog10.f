c----------------------------------------------------------------------
C      @(#)sordgk.lrlog10.f (derive de sordgk3.f1.3 5/26/87).

      Subroutine SORDGK(P,CP,Resol,Residu,P0,CP0,DP0,D0,CD0,Dinit,
     *     n,m,nd,mp,titre,nvar, Prof,distcor,Coupl,ncouch,nmode,
     *                  nfiltr,ienv,Ndata,Sig2r,Sig2i,nbs,nbt,nt)

*     Ce sous-programme affiche les resultats de l'inversion 
*      Tarantola et Valette (1982) dans un format compatible
*      avec la lecture de l'iteration suivante.
c----------------------------------------------------------------------
c   modif du 30 avril 1997
c   Version adaptee a l'inversion multi-sismogramme ainsi
c   qu'a l'inversion de log Mo (plus de scaling par e25).
c   modif du 1er Septembre 97
c   Calcul un qui enveloppe par periode et une perturbation
c   quadratique du modele pour les besoins de l'automatisation.
c----------------------------------------------------------------------

      parameter(mmp=200)

      Dimension P(mp),CP(mp,mp),Resol(mp,mp),Residu(nd),
     *          P0(mp),CP0(mp,mp),DP0(mp), D0(nd), CD0(nd,nd),
     *          Coupl(7,7),Ndata(nt,nd,nd),Dinit(nd), Prof(mp)
      Dimension Sig2r(nd), Sig2i(nd) 
	dimension ep(mmp),xinf(mmp),xsup(mmp)
	dimension ne(20),nebid(20),sum1env(20),sum2env(20)
	dimension nfiltr(*)
      real quim

      character iout*30, titre*80, form*80

	if(m.ne.nvar*ncouch+nbs)stop'sordgk2: Pb avec les mom. sismiques'

      write(*,*) 'nom du fichier de sortie pour dessin (lu11) ?'
      read(*,'(a)') iout
      write(*,*) 'nombre de couches a dessiner ?'
      read(*,*) ncoulim
      open(11,file=iout)
      write(11,'(a)') titre(1:40)//iout(1:10)
      write(11,*) ' 1'
      write(11,*) ncoulim, nvar
      do 33 i=1,ncouch
33    write(11,'(10e14.6)') Prof(i),(P(k),k=i,nvar*ncouch,ncouch)

c     Residu quadratique moyen avant inversion
      sum1=0.
	nv=0
      do 77 i=1,n
      resbid=D0(i)-Dinit(i)
	if(resbid.ne.0.) nv=nv+1
   77 sum1=sum1+ resbid*resbid/CD0(i,i)
      if (nv.ne.0.) then
      sum1=sqrt(sum1/nv)
      else 
      sum1=0.
      endif

c     Residu quadratique moyen apres inversion
      sum2=0.
	nw=0
      do 78 i=1,n
	if(Residu(i).ne.0.) nw=nw+1
   78 sum2=sum2+Residu(i)*Residu(i)/CD0(i,i)
	if(nw.ne.nv) write(0,*) "Pb dans le calcul du RMS"
      if (nw.ne.0.) then
      sum2=sqrt(sum2/nw)
      else 
      sum2=0.
      endif


      write(11,*) ' RMS normalise apres inversion:', sum2
c     Residus + Erreurs-donnees
      write(11,*) 1
      write(11,*) n
c Residu exact avant inversion, estime apres inversion, ecart-type a priori
      do 66 i=1,n
      write(11,'(i4,10e14.6)') i,D0(i)-Dinit(i),Residu(i),sqrt(CD0(i,i))
   66 continue

      write(11,*) ' Modele inverse'
      write(11,'(10f8.4)') (P(i), i=1,m)

      write(11,*) ' Erreurs-modele'
      write(11,'(10f8.4)') ( sqrt(CP(i,i)), i=1,m)

      write(11,*) 'Residu exact avant inversion:', sum1,
     *            ' Estime apres:',sum2

      write(10,'(a)') titre(1:40)//iout(1:10)
      write(10,*) nvar,ncouch
      do 1 i=1,ncouch
1     write(10,'(10e14.6)') Prof(i),(P0(k),k=i,nvar*ncouch,ncouch)

c------------------------------------------------------------------
c-- modif 30 avril 1997 on inverse pour log10Mo 
c	ecriture de log10Mo (apres correction d'echelle)
c     et ce dans le cas de nbs seismes
      do 11 i=nvar*ncouch+1,m 
11    write(10,*) P0(i)

      do 2 i=1,ncouch
2     write(10,'(10e14.6)')
     *             Prof(i),(DP0(k),k=i,nvar*ncouch,ncouch)

c	ecriture de l'erreur sur log moment sismique a priori
      do 22 i=nvar*ncouch+1,m
22    write(10,*) DP0(i)
c------------------------------------------------------------------

      do 3 i=1,ncouch
3     write(10,'(10e14.6)') Prof(i),(P(k),k=i,nvar*ncouch,ncouch)

c	ecriture du moment sismique obtenu a cette iteration
      do 34 i=nvar*ncouch+1,m
34    write(10,*) P(i),'      Moment sismique'

      write(10,*) distcor,'     distance de correlation'
	write(form,'(a,i3,a)') '(',nvar,'f8.4,a,i3)'
      do 4 i=1,nvar
4     write(10,form) (coupl(i,j),j=1,nvar),'        couplages',i
c--------------------------------------------------------------------
c   reecritures dans le fichier tv.n des donnees reelles pour chaque 
c   trajets.

      n1=0
      n2=0
      do 44 inbt=1,nbt
      write(10,*) nmode,nfiltr(inbt),ienv, '       nmode, nfiltr, ienv'
      do 5 imode=1,nmode
      do 5 ifiltr=1,nfiltr(inbt)
      ndat=ndata(inbt,imode,ifiltr)
      write(10,'(i3,(10e14.6))') ndat,(Sig2r(i),Sig2i(i),i=n1+1,n1+ndat)
      n1=n1+ndata(inbt,imode,ifiltr)
5     continue
      do 6 imode=1,nmode
      do 6 ifiltr=1,nfiltr(inbt)
      ndat=ndata(inbt,imode,ifiltr)
      write(10,'(i3,(10e14.6))') ndat,(sqrt(CD0(i,i)),i=n2+1,n2+ndat)
      n2=n2+ndata(inbt,imode,ifiltr)
6     continue
44    continue

c--------------------------------------------------------------------
c     Ici on affiche le tableau d0(i) calcule dans lecdgk.
c     on introduit dans la boucle le calcul du residu par type
c     de donnee(enveloppes ou phase).
c     ne(ifiltr) est le nombre de donnees d'enveloppes pour la periode ifiltr
c     et np est le nombre de donnees de phase.
      n1=0
      nfiltrmax=0
      np=0
      sum2pha = 0
      do 45 inbt=1,nbt
      write(10,*) ' Donnees inversees'
      do 15 imode=1,nmode
      do 15 ifiltr=1,nfiltr(inbt)
      if((inbt.eq.1).and.(imode.eq.1)) then
         nebid(ifiltr)=0
         ne(ifiltr)=0
         sum1env(ifiltr)=0
         sum2env(ifiltr)=0
      endif
      if (nfiltr(inbt).gt.nfiltrmax) nfiltrmax=nfiltr(inbt)
      ndat=ndata(inbt,imode,ifiltr)
      write(10,'(i3,(10e14.6))') ndat,(D0(i),i=n1+1,n1+ndat)
      if (ndat.ge.3)then
       do 14,i=1,ndat
       nr=n1+i
       reresbid=D0(nr)-Dinit(nr)
       if(reresbid.ne.0.) nebid(ifiltr)=nebid(ifiltr)+1
       sum1env(ifiltr)=sum1env(ifiltr)+reresbid*reresbid/CD0(nr,nr)
       if(Residu(nr).ne.0.) ne(ifiltr)=ne(ifiltr)+1
14     sum2env(ifiltr)=sum2env(ifiltr)+Residu(nr)*Residu(nr)/CD0(nr,nr)
      elseif (ndat.eq.1) then
       nr=nr+1
       if(Residu(nr).ne.0.) np=np+1
       sum2pha=sum2pha+Residu(nr)*Residu(nr)/CD0(nr,nr)
      endif
      n1=n1+ndata(inbt,imode,ifiltr)
15     continue

45    continue

      do 16 ifiltr=1,nfiltrmax
      if (nebid(ifiltr).ne.0) then 
        sum1env(ifiltr)=sqrt(sum1env(ifiltr)/nebid(ifiltr))
      else
        sum1env(ifiltr)=0
      endif
      if (ne(ifiltr).ne.0) then 
        sum2env(ifiltr)=sqrt(sum2env(ifiltr)/ne(ifiltr))
      else
        sum2env(ifiltr)=0
      endif
16    continue

      if (np.ne.0) sum2pha=sqrt(sum2pha/np)
      else sum2pha = 0
c--------------------------------------------------------------------
      write(10,*) ' Modele inverse'
      write(10,'(10f8.4)') (P(i), i=1,m)

      write(10,*) ' Erreurs-modele'
      write(10,'(10f8.4)') ( sqrt(CP(i,i)), i=1,m)
c ----------------------------------------------------------------------
c   Le format d'ecriture des residus a ete chamge de 10f8.4 a 10f10.4
c   pour les besoins d'un essai. il n'est pas sur que ce nouveaux format
c   sera compatible Rceux imposes en lecture dans les autres pgms.
      write(10,*) ' Residus'
      write(10,'(10f10.4)') (Residu(i), i=1,n)
c----------------------------------------------------------------------
      trace=0.
      do 23 i=1,m
   23 trace=trace+Resol(i,i)

c     Perturbation quadratique du modele (quim)
c     P est le modele inverse P(i) le modele a priori
c     soit le modele de l'ite de depart.
c  modif eric 4 fevr 98
c  dans le cas ou nvar=3 on inverse Vs Xi et Q
c  en love rayleigh simultanement. 
c  Sigma(log10Q) a priori est alors choisi tres large
c  pour eviter les pb d'incompatibilite love/Rayleigh
c  Il en resulte que le Q inverse a plus de mal a converger
c  et peut faire varier considerablement le khim entre 2 iterations.
c  Puisque il n'est pas interprete on ne le prendra pas en
c  compte dans le calcul de khim et seul les 2 premiers
c  parametres (Vs et Xi) seront utilises. Ainsi dans le pgm
c  d'inversion automatique seule la convergence de Vs et Xi sera
c  verifiee a travers le khim.
c  On ne tient pas compte non plus de la perturbation de Mo pour 
c  le calcul de ce khim.
      nvarbid=nvar
      if (nvar.eq.3) then
        write(*,*)'ATTENTION: on suppose que les 3 parametres inverses'
        write(*,*)'sont Vs,Xi et log10Q et que l\'ecart-type a priori  '
        write(*,*)'sur log10Q est large=>on calcule la perturbation   '
        write(*,*)'quadratique du modele seulement pour Vs et Xi car  '
        write(*,*)'log10Q ne sera pas interprete.                     '
        nvarbid=2
      endif
      quim=0.
      do 24 i=1,nvarbid*ncouch
   24 quim=quim+ (P(i)-P0(i))*(P(i)-P0(i))/CP(i,i)
      quim=sqrt(quim/m)

      write(10,*) trace,' degres de liberte'
c     write(10,*) 'Residu exact avant inversion:', sum1,
c    *            ' Estime apres:',sum2,' khim:',quim
      write(10,1000)sum1,sum2,quim
1000  format(' Residu exact avant inversion:',f9.6,
     *' Estime apres:',f9.6,' khim:',f9.6)
c     write(10,*) 'Residu ex. av. inv. enveloppes T1: ', sum1env(1),
c    *            ' Estime apres:',sum2env(1)
      write(10,1001)sum1env(1),sum2env(1)
1001  format(' Residu  ex. av. inv. enveloppes T1:',f9.6,
     *' Estime apres:',f9.6)
c     write(10,*) ' Residu estime phase:         ', sum2pha
      write(10,1002)sum2pha
1002  format(' Residu estime phase: ',f9.6)

      write(*,*) ' entrez 1 pour avoir la correlation'
      write(*,*) ' (dans le fichier-dessin)'
      read(*,*) ibid
      if(ibid.eq.1) then

      write(11,'(a)') ' Correlation a priori'
      do 7 i=1,m
      write(11,'(10e14.6)') (CP0(i,j)/sqrt(CP0(i,i)*CP0(j,j)), j=1,m)
    7 continue

      write(11,'(a)') ' Correlation a posteriori'
      do 8 i=1,m
      write(11,'(10e14.6)') (CP(i,j)/sqrt(CP(i,i)*CP(j,j)), j=1,m)
    8 continue
      end if

      write(*,*) ' Entrez 1 pour avoir la resolution'
      write(*,*) ' (sans dimension, dans le fichier-dessin)'
      read(*,*) ibid
      if(ibid.eq.1) then

c     On calcule la resolution en relatif, car plus facilement
c     interpretable lorsqu'on a des parametres physiques
c     differents (pb d'unites, par exemple)

	call epprof(ep,Prof,xinf,xsup,nvar,ncouch)
c-------------------------------------------------------- 
c initialement,la resolution etait calculee de i=1,m ce qui
c induit une erreur pour le moment sismique car ep(j) est
c definit de 1 a nvar*ncouch et m de 1 a nvar*ncouch+nbs
c ep(j) prenait donc la valeur 0 dans les nbs dernieres
c valeurs de m, Resol(i,j) une valeur infinie d'ou
c l'apparition d'un message d'erreur.C'est pourquoi
c les boucles vont maintenant de 1 a m-nbs.
c--------------------------------------------------------
      write(11,'(a)') ' Resolution relative'
      do 9 i=1,m-nbs
      write(11,'(10e14.6)')
     &	(Resol(i,j)*sqrt(CP0(j,j)/CP0(i,i))/ep(j), j=1,m-nbs)
    9 continue
      write(11,'(f6.2,a)') trace,' degres de liberte'

      write(11,'(a)') ' Resolution absolue'
      do 10 i=1,m-nbs
      write(11,'(10e14.6)') (Resol(i,j), j=1,m-nbs)
   10 continue
      write(11,'(f6.2,a)') trace,' degres de liberte'
      end if

      return
      end
