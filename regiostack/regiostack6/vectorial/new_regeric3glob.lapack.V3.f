c Use frt on vpp and frtvpp on covpp
c frt -Wv,-Of -Ne,cosdel new_regeric.f -o regeric
c For vectorization listings do:
c frt -Wv,-md -Ps -Z listing new_regeric.f -o regeric
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c     continuous regionalization                                       c
c    Version vectorisable prevue pour tourner sur le VPP300 de l'anu   c
c    les variables theta(300000),phi(300000), s2az(300000),            c
c    c2az(300000),fcov(300000)                                         c
c    sont stockees dans un seul tableau de dimension                   c
c    (nombre de points par trajet)*(nombre de trajets)                 c
c    Avec 3000 trajets et un step de 0.5 degres, en admettant que pour c
c    une tomo regionale la longueur moyenne de trajet est de 6000 km   c
c    on peut dire (a la louche) qu'il nous faut une dimension de       c
c    120*3000 soit 360 0000. C'est tres grand, on se contente de       c
c    300 000 pour l'Australie en stopant le pgm si cette dimension est c
c    insuffisante.        
c    Modif eric 9/7/2001                                              c
c    dimensions augmentees pour la tomo Globale
c    nbre de trajets maxi monte a 36000, pas 5 degres, longueur de
c    trajet moyen 10000 km soit (20 * 36 000 = 1 080 000)
c    On garde  2 500 000 pour etre sur d'etre large 
c    Ces tableaux intervienent aussi dans la routine grille_vpp        c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     references: montagner(1986, annales geophysicae, 4: b3, 283-294).c
c     modifications pour antarctique mars 1991
c     sortie de l anisotropie(juillet 91)
c     calcul de khi(aout 91)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  modifs eric.
c     modifie pour sortir les parametres dan1 et dan2 de la 
c     subroutine grille(eric avril 95).
c-----------------------------------------------
c   Dernieres modif apportees: On supprime la possibilite d'inverser
c   une seconde fois en augmentant l'erreur sur les donnees mal fitees.
c   Les derivees partielles des parametres anisotropes sont corrigees.
c   ON SORT Viso A1 A2 (Smith et Dahlen) de ce pgm alors que les parametres
c   inverses sont les (-A*vitesse a priori/vitesse inversee**2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     parameter (NDIM=200,NTRAJ=36000,NVECT=1500000)
      parameter (NDIM=200,NTRAJ=18000,NVECT=750000)
c     parameter (NLAT=180,NLON=360)
      parameter (NLAT=90,NLON=180)
      parameter (NTDIM=2*NDIM+8)
c     dimension ux(NTRAJ),bcd(NTRAJ,2),ulf(NLAT,NLON),az(181,NTRAJ)
      dimension ux(NTRAJ),ulf(NLAT,NLON),az(181,NTRAJ)
      dimension t(200),delt(NTRAJ),dsi(NTRAJ),kdi(NTRAJ)
      dimension uld(NTRAJ),ed(NTRAJ),w1(5,NTRAJ),ndat(200)
      dimension fcov(NVECT),s1(181,NTRAJ),s1a(181),dan1(NLAT,NLON)
      dimension cov(NTRAJ,NTRAJ),s(NTRAJ,NTRAJ),ww(NTRAJ),res(NTRAJ)
!      dimension  s2az(NVECT), c2az(NVECT), theta(NVECT), phi(NVECT)
      integer s_path(NTRAJ), e_path(NTRAJ)
      dimension v(NTRAJ),vv(NTRAJ),pac(NLAT,NLON),sig(NLAT,NLON),
     *dan2(NLAT,NLON)
c     dimension s2(NTRAJ,NTRAJ),verif(NTRAJ,NTRAJ)
c     dimension bcd2(NTRAJ,2),jt(200),ind(NTRAJ)
      dimension jt(200),ind(NTRAJ)
      dimension sigan1(NLAT,NLON),sigan2(NLAT,NLON)
      dimension ed1(NTRAJ)
c     modif eric
c     dimension da1(NLAT,NLON),da2(NLAT,NLON)
      dimension resn(NTRAJ)

      character*7 staz(NTRAJ)
      real randerr

      character*32 fileout,fianiout,ficdata
      common/c1/tetai(181,NTRAJ),phii(181,NTRAJ)
      common/c2/s
      common/c3/pac,sig,ulf
      common/c4/az
      common/c5/dan1,dan2
      common/c6/xo,yo,nxo,nyo,pxyo
      common/c61/xinf,yinf,nx,ny,pxy
      common/c7/w(NTDIM,NTRAJ)

      common/dbs/theta(NVECT),phi(NVECT),s2az(NVECT),
     *c2az(NVECT),
     .             s_path, e_path

      real work(NTRAJ)
      integer ipiv(NTRAJ)
      integer info


      pi=3.141592653
      rad=pi/180.


      print *,'Some uncertainties (source focal mechanism or '
      print *,'location) may not be contained in the velocity '
      print *,'error provided in the data file. Enter your '
      print *,'estimation for this error (0.01 km/s has been '
      print *,'found to be a good value for path average Vs'
      print *,'models obtained after waveform modelling) '
      read(*,*)randerr
      print *,randerr
      print *,'enter sampling step along the paths'
      read(*,*)d0
      print *,d0
      print *,'enter number of paths'
      read(*,*)ntra
      print *,ntra
      d1=d0*rad
      print *,' enter xo,yo,nxo,nyo,pxyo'
      read(*,*) xo,yo,nxo,nyo,pxyo
      print *,xo,yo,nxo,nyo,pxyo
c  la dimension de ww=5 fois celle fe ux
c  dimensions of arrays (x0,y0 :initial colatitude and longitude)
c  nxo,nyo :nb of points in latitude,longitude
c  pxyo :step in lat,long
c  these values can be different from xinf yinf nx ny pxy because
c  it is possible to compute the final model in an area smaller
c  than the region covered by the paths (zooming)
c
      print *,'enter xinf,yinf,nx,ny,pxy'
      read(*,*) xinf,yinf,nx,ny,pxy
      print *,xinf,yinf,nx,ny,pxy
c
c xinf,yinf: colatitude and longitude of the origin point
c            located at the north-western corner of the grid
c
c
c  xo,yo                        nyo
c      o------------------------------------------>
c      |
c      |                         ny           
c      |   xinf,yinfo--------------------------->
c      |            |
c      |            |
c      |       nx   |
c      |            |
c      |            |
c      v            v
c
      print *,'enter nper,nbper,indices periodes a traiter'
      read(*,*)nper,nbper,(jt(i),i=1,nbper)
      print *,nper,nbper,(jt(i),i=1,nbper)
      read(*,*)(t(i),i=1,nper)
      print *,(t(i),i=1,nper)
      do 9 i=1,ntra
      do 9 j=1,NTDIM
   9  w(j,i)=10.
c
c   reading of dataset in a formatted file data
c   tam-agr-bca-                       forinstance
c   after data are read. the are ordered in an array w(32,ntra)
c   w(1,i),w(2,i) : m1(teta1,phi1)   w(3,i),w(4,i):m2(teta2,phi2)
c                    (epicenter)                    (station)
c   w(5,i): phipl-phi(average azse) azimuth with respect to the direction 
c   of the plate velocity
c   not used in this version
c   w(7....NDIM+6,i)  group or phase velocities measured at different periods.
c   w(23...NTDIM,i)  errors on group or phase velocities.
c   in that case, we had 16 periods from 10s to 250s:
c   delt epicentral distance for each path i, dsi incremental step on
c   each path i. kdi number of sampled points along each path i.
c   w data array.
c   tetai(colat(mi),i) phii(long(mi),i)
c
      print *,'enter data file name'
      read(5,'(a32)') ficdata
      print 999,ficdata
      open(10,file=ficdata,form='formatted')
      call lecdatV3(randerr,nper,staz,delt,kdi,dsi,ntra,d1,t)
c
c  in lecdat, teta must be between 0. and teta and phi(radians)
c  covariance matrix for the whole paths is computed at each period
c  in order enable correlation length to be period dependent.
c  however, as this computation is quite long,it can be
c  computed at once if we forget this period dependency.
c  dcor, dcoran correlation lengths.
 999  format(a32)
      print *,'enter final velocity file'
      read(5,999)fileout
      print 999,fileout
      open(7,file=fileout)
      print *,'enter final anisotropy file'
      read(5,999)fianiout
      print 999,fianiout
      open(9,file=fianiout)
      print *,'enter dcor,dcoran (in km) for velocity and anisotropy'
      read(*,*)dcor,dcoran
      print *,dcor,dcoran
      dcor=dcor*rad/111.1949
      dcoran=dcoran*rad/111.1949
c
      k1=1
c-----------------------------------------------------------------------
c      loop on periods
c-----------------------------------------------------------------------
      do 5 ktt=1,nbper
      ntt=jt(ktt)
      k2=ntt
      ifois=0
c  in vinit,one reads starting velocity and anisotropy
c  at a given period, data are sorted because the number
c  of data can be period dependent
      print *,' '
      do 55 jk=k1,k2
      call vinit(sigma,anis1,anis2,sigan,ut0,tt)
 55   k1=k2+1
      write(7,*)tt
      write(9,*)tt
      k=0
      do 6 j=1,ntra
         nk=ntt+6
         nkk=ntt+NDIM+6
         if(w(nkk,j).gt.2.) go to 6
         k=k+1
         ind(k)=j
         uld(k)=w(nk,j)
         ed(k)=w(nkk,j)
         ed1(k)=ed(k)/uld(k)**2
         kdi(k)=kdi(j)
         dsi(k)=dsi(j)
         delt(k)=delt(j)
c        bcd2(k,1)=bcd(j,1)
c        bcd2(k,2)=bcd(j,2)
         staz21(k)=staz1(j)
         staz22(k)=staz2(j)
         kk1=kdi(j)
         do 7 ki=1,kk1
            az(ki,k)=az(ki,j)
            tetai(ki,k)=tetai(ki,j)
            phii(ki,k)=phii(ki,j)
    7    continue
         w1(1,k)=w(1,j)
         w1(2,k)=w(2,j)
         w1(3,k)=w(3,j)
         w1(4,k)=w(4,j)
         w1(5,k)=-uld(k)*cos(2.*w(5,j))
    6 continue
      ndat(ntt)=k
      print *,' '
      print 2007,t(ntt),k,dcor
 2007 format('periode,nda,dcor',f6.2,i6,f6.2)
c-------------------------------------------------------
c modif eric
c on calcul les moyennes et ecarts-types en lenteur
c car on inverse en lenteur.
      print *,'selected data '
c     print 2003,(uld(j),j=1,k)
      do 4 it=1,k
  4   ux(it)=1./uld(it)
      print 2003,(ux(i),i=1,k)
c     call me(ux,xm,xe,k)
      call me(uld,xm,xe,k)
c-------------------------------------------------------
CJJL      print 1100,xm,xe
      write(*,*) 'moyenne et ecart-type de uld (1/donnees)',xm,xe
      varmoy=xe*xe
      write(*,*)'variance-donnees p.r. a la moyenne des mesures:',varmoy
      write(*,*) 'modele initial= ',ut0,
     &   '           moyenne des mesures= ',1/xm

      x0=xe
 1100 format(' initial values : xm= ',f7.3,' xe= ',f7.4)
      nda=ndat(ntt)

      ri2dcor2 = 1.0/(2*(dcor**2))
      ri2dcoran2 = 1.0/(2*(dcoran**2))
      sigan2riut02 = (sigan**2)*(1/ut0)*(1/ut0)
      sigma2 = sigma**2
      
!      open(unit=13,file='looplog')
!      write(13,*) 'ri2dcor2= ', ri2dcor2
!      write(13,*) 'ri2dcoran2= ',  ri2dcoran2
!      write(13,*) 'sigan2riut02=',sigan2riut02
!      write(13,*) 'sigma2 =',sigma2
      ik = 0
      do i = 1,nda
         ik = ik+1
         s_path(i) = ik
         ki=1
         tetai(ki,i)=w1(1,i)
         phii(ki,i)=w1(2,i)
         theta(ik) = w1(1,i)!tetai(ki,i)
         phi(ik) = w1(2,i)!phii(ki,i)
         do ki = 2,kdi(i)
            ik = ik+1
            s2az(ik) = sin(2.*az(ki,i))
            c2az(ik) = cos(2.*az(ki,i))
            theta(ik) = tetai(ki,i)
            phi(ik) = phii(ki,i)
         enddo
         e_path(i) = ik
      enddo

      write(*,*) 'total points = ',ik
      if( ik > 1500000) stop 399

      do i=1,nda
!         if(mod(i,100) .eq. 0) then
!            write(13,*) i, ' of ',nda, ' vector length ', kdi(i)
!            call flush(13)
!         endif

c
c loop on paths j_<i then double loop on every
c point between m1(i),m2(i), and m1(j), m2(j)
c
         kki=kdi(i)
         kkj=kdi(j)
         do ki=1,kki
            ik = s_path(i) + ki - 1
          !  s1(ki)=0.
                                !      s1a(ki)=0.
            s2azkii = s2az(ik)
            c2azkii = c2az(ik)
            ctik=cos(theta(ik))
            stik=sin(theta(ik))

            do jk = 1,e_path(i)
               
!               ctik=cos(theta(ik))
               ctjk=cos(theta(jk))
!               stik=sin(theta(ik))
               stjk=sin(theta(jk))
               sdmimj=ctik*ctjk+(stik*stjk)*cos(phi(ik)-phi(jk))
               sdmimj = max( min( sdmimj, 1.0), -1.0 )
               sdmimj=acos(sdmimj)
               sdmimj=sdmimj**2
               ex=sdmimj*ri2dcor2 !/(2*(dcor**2))
c-------------------------------------------------------------------
c modif eric
c  Ici c'est dcoran et non dcor
c     exa=sdmimj/(2*(dcor**2))
               exa=sdmimj*ri2dcoran2 !/(2*(dcoran**2))
c-------------------------------------------------------------------
c              if(ex.gt.25.) then
               if(ex.gt.4.59) then
                  fcov(jk)=0.
               else
                  fcova1 =exp(-exa)*sigan2riut02 
                  fcova2 =fcova1*s2azkii*s2az(jk)
                  fcova1 =fcova1*c2azkii*c2az(jk)
                  fcov(jk)=sigma2*exp(-ex)+fcova1+fcova2
               endif
            enddo
            
            do j = 1,i
               s1(ki,j) = 0.5*fcov(s_path(j))
               do jk = s_path(j)+1,e_path(j)-1
                  s1(ki,j)=s1(ki,j)+fcov(jk)
               enddo
               s1(ki,j) = (s1(ki,j)+0.5*fcov(e_path(j)))*dsi(j)
            enddo
            
         enddo

         do j = 1,i
            call trapez(kdi(i),s1(1,j),dsi(i),co)
            cov(i,j)=(co)/(delt(i)*delt(j))
            cov(j,i)=cov(i,j)
            s(i,j)=cov(i,j)
            s(j,i)=s(i,j)
c
c  at this level, the double integral dsi dsj*exp(-sin(delt)**2
c /2*dcor**2   has been computed.
c
            if (i.eq.j) s(i,j)=(ed(i)**2)+s(i,j)
         enddo
c  this integral is  the hessian of the system.
      enddo
c     print *,'cov(i,j)'
c     do 102 i=1,nda
c102  print 2003,(cov(i,j),j=1,nda)
c calculation of s-1
c     do i=1,nda
c        do j=1,nda
c           s2(i,j)=s(i,j)
c           verif(i,j)=0.
c        enddo
c     enddo
c     e=0.

c      call ma22b(s,15000,nda,ww,e)
      call SGETRF(nda,nda,s,NTRAJ,ipiv,info)
      call SGETRI(nda,s,NTRAJ,ipiv,work,NTRAJ,info)
c       call SPOTRF('U',nda,s,NTRAJ,info)
c       call SPOTRI('U',nda,s,NTRAJ,info)       

c     do i=1,nda
c        do j=1,nda
c           do k=1,nda
c              verif(i,j)=verif(i,j)+s2(i,k)*s(k,j)
c           enddo
c        enddo
c     enddo
      
c     print *,'verif',nda
c     aa=0.
c     do i=1,nda
c        aa=aa+verif(i,i)
c     enddo
c     print*, 'verif',aa
      

c    forward problem resolution.
      do 130 i=1,nda
         v(i)=0.
         call pbdirc(w1(5,i),uld(i),v(i),pac,az,kdi(i),delt(i),dsi(i),i)
  130 continue
      print*,'vecteur v (d-g(p) en lenteur):'
      print 2003,(v(i),i=1,nda)
CJJL  Calcul de la variance des donnees par rapport aux synthetiques
CJJL  du modele initial (et non par rapport a la moyenne des donnees)
      varini=0.
      do 139 i=1,nda
         varini=varini+v(i)*v(i)
  139 continue
      varini=varini/nda
      write(*,*)'variance-donnees p.r. au modele initial:       ',varini
CC modif eric : calcul du qui initial en lenteur 
      call khi (v,ed,qui,nda)
      print*, 'qui initial :',qui
      do 140 i=1,nda
      vv(i)=0.
      do 140 j=1,nda
 140  vv(i)=vv(i)+s(i,j)*v(j)
c     print *,'vvi'
c     print 2012,(vv(i),i=1,nda)
      call me(vv,xm,xe,nda)
c     print 1100,xm,xe
c
      call grille_vpp(vv,dsi,delt,kdi,nda,dcor,dcoran,sigma,sigan,uld,
     1w1,sigan1,sigan2,ut0)
c
      do 94 i=1,nda
      res(i)=0.
      re=res(i)
      call pbdirc(w1(5,i),uld(i),re,ulf,az,kdi(i),delt(i),dsi(i),i)
c   modif eric
c   On calcul le residus en lenteur (inversion 1/Q)
c   et non en Q
c     res(i)=re*(ux(i)**2)
c     resn(i)=sqrt((res(i)/ed1(i))**2)
      res(i)=re
      resn(i)=sqrt((res(i)/ed(i))**2)
  94  continue
	print*,'residus normalises par trajets (et par couche):'
      print 1007,(staz21(i),staz22(i),resn(i),i=1,nda)
c     print 1007,((bcd2(i,j),j=1,2),resn(i),i=1,nda)
      icomp=0
c------------------------------------------------------------------
c  modif eric
c     do 96 i=1,nda
c     esup=abs(res(i))
c     if(esup.lt.0.1)go to 96
c     ed(i)=5.*ed(i)
c     ed1(i)=2.*ed1(i)
c     icomp=icomp+1
c 96  continue
c------------------------------------------------------------------
      do 95 i=1,nda
      do 95 j=1,nda
      cov(i,j)=0.
  95  s(i,j)=0.
c------------------------------------------------------------------
c  modif eric
c     if(icomp.eq.0.or.ifois.eq.1)go to 97
c     print *,'on recommence'
c     ifois=ifois+1
c     go to 11
c 97  call me(res,xm,xe,nda)
c------------------------------------------------------------------
c--x0 et xe ont ete passe en lenteur, normalement redvar esten 1/Q.
      call me(res,xm,xe,nda)
CJJL  redvar=(x0-xe)/x0
CJJL
CJJL  print 1100,xm,xe
      write(*,*) 'moyenne et ecart-type de res (residus fin)',xm,xe
      varres=xm*xm+xe*xe
      write(*,*)'variance-donnees p.r. a la moyenne des mesures:',varmoy
      write(*,*)'variance-donnees p.r. au modele final:         ',varres
      redvar =(varmoy-varres)/varmoy
      redvar2=(varini-varres)/varini

      print *,'red. variance(%val moy des donnees en lenteur)',redvar
      print *,'red. variance(% modele initial)     en lenteur',redvar2
c     call khi(res,ed1,qui,nda)
c------------------------------------------------------------------
c modif eric : ed est l'erreur en 1/Q, ed1 celle en vitesse
c On doit donc calculer un qui en 1/Q
c pour khim : sigma1 est passe en lenteur dans vinit. On le garde
c comme ca.
      call khi(res,ed,qui,nda)
CJJL      print*,'qui  (ecart-type sur les donnees)=',qui
      print*,'qui  (res. quadr. normalise des donnees en lenteur)=',qui
c     calcul de quim
      quim=0.
c     sigma1=sigma*(ut0**2)
      do i=1,nxo
        do j=1,nyo
c          quim=quim+((ulf(i,j)-pac(i,j))**2)/sigma1**2
           quim=quim+(((1./ulf(i,j))-(1./pac(i,j)))**2)/sigma**2
        enddo
      enddo
      quim=sqrt(quim/(nxo*nyo))
CJJL      print*,'khim (ecart-type sur le modele)=',quim
c     print*,'khim*nxo*nyo=',quim,'    khim=',quim/(nxo*nyo)
      print*,'khim (pert.quadr. normal. du modele en lenteur)= ',quim
c------------------------------------------------------------------
c     ecriture du resultat
      do 127 i=1,nxo
 127  write(7,2014)(ulf(i,j),j=1,nyo)
      do 126 i=1,nxo
 126  write(7,2014)(sig(i,j),j=1,nyo)
c-------------------------------------------------------------------
c modif eric Dec 1995
c Les derivees partieles des parametres anisotropes sont (-cos2teta)/L*ut0
c avec  ut0=vitesse initiale(cf grille). Ceci implique(voir these)que les 
c parametres inverses sont dan=(A1*ut0)/(ulf*ulf) ulf est la vitesse inversee
c et non la lenteur car ulf est converti en vitesse dans la routine grille.
c dans ce cas les parametres Ai (coefficients de smiths et Dahlen)sont: 

      do 1210 i=1,nxo
      do 1200 j=1,nyo
      dan1(i,j)=dan1(i,j)*(ulf(i,j)*ulf(i,j))/ut0
      dan2(i,j)=dan2(i,j)*(ulf(i,j)*ulf(i,j))/ut0
1200   continue
1210  continue
 
c apres cette conversion, le pgm sort les coefficients bruts de Smith et Dahlen
c tels que donnes dans la formule C= C0+A1cos2phi+A2sin2phi
c-------------------------------------------------------------------
c     modif eric.                          
c     ecriture des parametres anisotropes.
      do 125 i=1,nxo
 125  write(9,2014)(dan1(i,j),j=1,nyo)
      do 124 i=1,nxo
 124  write(9,2014)(dan2(i,j),j=1,nyo)
      do 123 i=1,nxo
 123  write(9,2014)(sigan1(i,j),j=1,nyo)
      do 122 i=1,nxo
 122  write(9,2014)(sigan2(i,j),j=1,nyo)
c------------------------------------------------------------------
  5   continue

 1007 format(6(a4,a4,1x,f5.2,1x))
 2003 format(10f7.4)
 2009 format(20f6.4)
 2012 format(10f9.2)
 2014 format(20(f10.6,1x))
      stop
      end
c------------------------------------------------------------------
      subroutine coordo(tetam1,phim1,tetam2,phim2,tt,pht,deltm,azm)
c coordo computes coordinates of points located on a great circle
c  between 2 points m1 and m2. deltm is the epicentral distance
      pi=3.141592653
      s=sin(deltm)
      c=cos(deltm)
      cot=cosdel(tetam1,phim1,tetam2,phim2)/sindel(tetam1,phim1,
     *tetam2,phim2)
      alpha=c-s*cot
      beta=s/sindel(tetam1,phim1,tetam2,phim2)
      ct=alpha*cos(tetam1)+beta*cos(tetam2)
c     write(*,*)'DELTM',deltm,' S ',s,' C',c,' COT',cot,
c    *' ALPHA',alpha,' BETA',beta,' CT',ct
      if(abs(ct).gt.1)then
        write(*,*)'     abs(CT) > 1,ct= ',ct
c       ct=sign(1,ct)
        if(ct.lt.0)ct=-1.
        if(ct.gt.0)ct=1.
        write(*,*)'     new ct= ',ct
      endif
c MODIF ERIC-JJ 6/09/2001 (tt=colat point courant)
      if (ct.eq.1.0) then
c       CAS DU PASSAGE POLE NORD
        tt=0.
        pht=phim1
        azm=0.
      else if(ct.eq.-1.0)then
c       CAS DU PASSAGE POLE SUD
        tt=180.
        pht=phim1
        azm=0.
      else
        tt=acos(ct)
        st=sqrt(1-ct**2)
        st1=sin(tetam1)
        st2=sin(tetam2)
        sp1=sin(phim1)
        sp2=sin(phim2)
        ct1=cos(tetam1)
        ct2=cos(tetam2)
        cp1=cos(phim1)
        cp2=cos(phim2)
        cpht=(alpha*st1*cp1+beta*st2*cp2)/st
        if(cpht.ge.1.)cpht=1.
        if(cpht.lt.-1.)cpht=-1.
        pht=acos(cpht)
        spht=(alpha*st1*sp1+beta*st2*sp2)
        if(spht.lt.0.)pht=-pht
        om1=sp1*st1*ct2-ct1*sp2*st2
        om2=ct1*cp2*st2-ct2*cp1*st1
        om3=st2*st1*sin(phim2-phim1)
        vn=spht*om1-cpht*om2
        ve=-om1*ct*cpht-om2*ct*spht+om3*st
        azm=atan2(ve,vn)
      endif
      if(azm.lt.0.)azm=2.*pi+azm
      rad=0.017453293
      te=tt/rad
      p=pht/rad
      azz=azm/rad
c       write(*,*)'TE P AZZ ',te,p,azz
c     print 2020,te,p,azz,spht
c2020 format(4f8.2)
      amap=amax1(phim1,phim2)
      amip=amin1(phim1,phim2)
c     if(pht.gt.amap.or.pht.lt.amip) print 2021,pht,amap,amip
 2021 format('mauvais trajet:pht,amap,amip',3f8.5)
      return
      end
c----------------------------------------------------------
      subroutine trapez(npoint,f,ds,sint)
      dimension f(*)
      sint=0.
      nmax=npoint-1
      do 1 ns=2,nmax
  1   sint=sint+2.*f(ns)
      sint=sint+f(1)+f(npoint)
      sint=sint*ds*0.5
      return
      end
c------------------------------------------------------------
c     lecdatV3 version for new intomodesVs file: no station list
c     at the start
c------------------------------------------------------------
      subroutine lecdatV3(randerr,nper,staz,delt,kdi,dsi,ndat,d1,t)
c     parameter (NDIM=200,NTRAJ=36000,NVECT=1500000)
      parameter (NDIM=200,NTRAJ=18000,NVECT=750000)
      parameter (NTDIM=2*NDIM+8)
      dimension t(200),ug(200),edat(200),az(181,NTRAJ)
      dimension delt(*),dsi(*),kdi(*)
      dimension deltei(181),wlat(600),wlon(600)

      character*7 staz(NTRAJ)

      common/c1/tetai(181,NTRAJ),phii(181,NTRAJ)
      common/c4/az
      common/c7/w(NTDIM,NTRAJ)

      if(nper.gt.NDIM) stop 'trop de periodes lues'
      rad=0.017453293

c     This subroutine reads for each epicenter path
c     the coordinates of the epicenters and stations,
c     the seismic velocities (path-average SV,
c     group or phase velocities) and their errors.
c     w(1,i),w(2,i)=lat lon of event 
c     w(3,i),w(4,i)=lat lon of station 

      nt=nper
      do 1 i=1,ndat
        read(10,'(a7)')staz(i)
c       write(*,'(a10)')staz

        read(10,*)w(1,i),w(2,i),w(3,i),w(4,i)
        if(w(2,i).le.-180.)then
           w(2,i)=w(2,i)+360.
        else if(w(2,i).gt.180.)then
           w(2,i)=w(2,i)-360.
        endif
        if(w(4,i).le.-180.)then
           w(4,i)=w(4,i)+360.
        else if(w(4,i).gt.180.)then
           w(4,i)=w(4,i)-360.
        endif

        w(1,i)=(90.-w(1,i))*rad
        w(2,i)=w(2,i)*rad
        in0=1
        fi=45.
        w(3,i)=(90.-w(3,i))*rad
        w(4,i)=w(4,i)*rad
        w(5,i)=fi*rad
        w(6,i)=nt
c     ug : group or phase velocity data for this path
c     edat error on ug
        read(10,*)(ug(j),j=1,nt)
        read(10,*)(edat(j),j=1,nt)
c     print 1999,(ug(j),j=1,nt)
c     print 1999,(edat(j),j=1,nt)
c1999 format(11f7.3)

        do 2 k=1,nt
          kk=k+6
          w(kk,i)=1/ug(k)
          kkk=k+NDIM+6
c         modif eric may 2003
c         Uncertainties in the origin time
c         in the source mecanism and possible mislocations
c         in the epicenter may not be reflected in the path
c         average a posteriori errors (this would be the case 
c         for waveform modelling).
c         Assuming the errors are normally distributed we add
c         them to the error in the data file :
          w(kkk,i)=sqrt(edat(k)**2+randerr**2)
          w(kkk,i)=w(kkk,i)/(ug(k)**2)
2       continue

        in=1
c
c       calculation of coordinates between m1 and m2
c       one computes the number of steps between m1i and m2i and dsi
c       is recalculated such that dsi*kdi=delt(i)
c
        delt(i)=acos(cosdel(w(1,i),w(2,i),w(3,i),w(4,i)))
        kdi(i)=int(delt(i)/d1)
        dsi(i)=delt(i)/kdi(i)
        ki=1
        kdi(i)=kdi(i)+1
        kd1=kdi(i)
        dds=dsi(i)/rad
        ddel=delt(i)/rad
        print 2005,kd1,dds,delt(i),ddel
 2005   format('nbr of pts,ds,delta(rad),delta(o)',i3,3f8.3)
        tm1=w(1,i)/rad
        pm1=w(2,i)/rad
        tm2=w(3,i)/rad
        pm2=w(4,i)/rad
        print 2006,tm1,pm1,tm2,pm2
 2006   format('t1 p1  t2 p2',4f8.2)

        do 4 ki=1,kd1
          gi=float(ki-1)
          deltei(ki)=gi*dsi(i)
          az(ki,i)=0.
          call coordo(w(1,i),w(2,i),w(3,i),w(4,i),tetai(ki,i),
     &                phii(ki,i),deltei(ki),az(ki,i))
          if(i.eq.5)then
c MODIF
c           write(*,*)'TEST',tetai(ki,i),phii(ki,i),deltei(ki),az(ki,i)
          endif
  4     continue
  1   continue
      return
      end
c------------------------------------------------------------
      subroutine lecdat(nper,staz1,staz2,delt,kdi,dsi,ndat,d1,t)
c     parameter (NDIM=200,NTRAJ=36000,NVECT=1500000)
      parameter (NDIM=200,NTRAJ=18000,NVECT=750000)
      parameter (NTDIM=2*NDIM+8)
c     dimension bcdb(600),bcd(NTRAJ,2)
      dimension t(200),ug(200),edat(200),az(181,NTRAJ)
      dimension delt(*),dsi(*),kdi(*)
      dimension deltei(181),wlat(600),wlon(600)

      character staz*10,nom1*4,nom2*6
      character*4 nomsta(600),staz1(NTRAJ)
      character*6 staz2(NTRAJ)

      common/c1/tetai(181,NTRAJ),phii(181,NTRAJ)
      common/c4/az
      common/c7/w(NTDIM,NTRAJ)
c   in this subroutine, we read the dataset(phase or group velocities)
c   the coordinates of the points along the paths are computed.
CCC	print *,"NDIM, NTDIM=", NDIM, NTDIM
      if(nper.gt.NDIM) stop 'trop de periodes lues'
      rad=0.017453293
      print *,'nsta: number of stations'
      read(10,*)nsta
      print *,nsta
      do 88 i=1,nsta
        read(10,1998)nomsta(i),wlat(i),wlon(i)
c       read(10,1998)bcdb(i),wlat(i),wlon(i)
c       MODIF eric pour avoir des lon entre 0 et 360 degres
c       if(wlon(i).lt.0.)wlon(i)=wlon(i)+360
      print 1998,nomsta(i),wlat(i),wlon(i)
  88  continue
c1998 format(a3,f8.3,f9.3)
 1998 format(a4,2x,2f9.4)

      nt=nper
      do 1 i=1,ndat
c     read(10,2000)(bcd(i,j),j=1,2)
c     print 2000,(bcd(i,j),j=1,2)
c2000 format(a4,a6)
      read(10,'(a10)')staz
        write(*,'(a10)')staz

        ir=index(staz,'.')
        staz1(i)=staz(:ir-1)
        staz2(i)=staz(ir+1:)
        write(*,'(a4,1x,a6)')'BCD ',staz1(i),staz2(i)
        icheck=0

      do 89 j=1,nsta
        nom1=nomsta(j)
        nom2=staz1(i)
        if(nom1(:lnblnk(nom1)).eq.
     *     nom2(:lnblnk(nom2))) then
           w(1,i)=wlat(j)
           w(2,i)=wlon(j)
           icheck=icheck+1
        endif
  89  continue
      if(icheck.ne.1) stop' station non trouvee '
      read(10,*)alate,alone
c     MODIF eric pour avoir des lon entre 0 et 360 degres
c     if(alone.lt.0.)alone=alone+360
c     print 2002,alate,alone
      w(1,i)=(90.-w(1,i))*rad
      w(2,i)=w(2,i)*rad
      in0=1
      fi=45.
      w(3,i)=(90.-alate)*rad
      w(4,i)=alone*rad
      w(5,i)=fi*rad
      w(6,i)=nt
      read(10,*)(ug(j),j=1,nt)
c  ug : group or phase velocity data for this path
c     print 1999,(ug(j),j=1,nt)
c  edat error on ug
c1999 format(11f7.3)
      read(10,*)(edat(j),j=1,nt)
c     print 1999,(edat(j),j=1,nt)
      do 2 k=1,nt
      kk=k+6
      w(kk,i)=1/ug(k)
      kkk=k+NDIM+6
      w(kkk,i)=edat(k)/(ug(k)**2)
  2   continue
      in=1
c2002 format(2f8.3,3f7.2)
c     print 2008,(w(k,i),k=7,NDIM+6)
 2008 format('1/ug',13f6.3)
c
c  calculation of coordinates between m1 and m2
c  one computes the number of steps between m1i and m2i and dsi
c is recalculated such that dsi*kdi=delt(i)
c
      delt(i)=acos(cosdel(w(1,i),w(2,i),w(3,i),w(4,i)))
      kdi(i)=int(delt(i)/d1)
      dsi(i)=delt(i)/kdi(i)
      ki=1
      kdi(i)=kdi(i)+1
      kd1=kdi(i)
      dds=dsi(i)/rad
      ddel=delt(i)/rad
      print 2005,kd1,dds,delt(i),ddel
 2005 format('nbr of pts,ds,delta(rad),delta(o)',i3,3f8.3)
      tm1=w(1,i)/rad
      pm1=w(2,i)/rad
      tm2=w(3,i)/rad
      pm2=w(4,i)/rad
      print 2006,tm1,pm1,tm2,pm2
 2006 format('t1 p1  t2 p2',4f8.2)
      do 4 ki=1,kd1
      gi=float(ki-1)
      deltei(ki)=gi*dsi(i)
      az(ki,i)=0.
      call coordo(w(1,i),w(2,i),w(3,i),w(4,i),tetai(ki,i),
     1phii(ki,i),deltei(ki),az(ki,i))
      if(i.eq.5)then
c MODIF
c     write(*,*)'TEST',tetai(ki,i),phii(ki,i),deltei(ki),az(ki,i)
      endif
  4   continue
  1   continue
      return
      end
c----------------------------------------------------------
      subroutine vinit(sigma,anis1,anis2,sigan,ut0,tt)
c     parameter (NLAT=180,NLON=360)
      parameter (NLAT=90,NLON=180)
      dimension pac(NLAT,NLON),sig(NLAT,NLON),dan1(NLAT,NLON),
     *dan2(NLAT,NLON)
      dimension ulf(NLAT,NLON)
      common/c3/pac,sig,ulf
      common/c5/dan1,dan2
      common/c6/xo,yo,nxo,nyo,pxyo
c     sigma:erreur a priori sur les parametres
      print *,'enter tt,ut0,sigma,anis1,anis2,sigan'
      read(*,*)tt,ut0,sigma,anis1,anis2,sigan
      print *,tt,ut0,sigma,anis1,anis2,sigan
      do 3 i=1,nxo
      do 3 j=1,nyo
      pac(i,j)=ut0
      sig(i,j)=sigma
      dan1(i,j)=anis1
  3   dan2(i,j)=anis2
      sigma=sigma/(ut0**2)
      print*,'sigma'  
      write(*,*) sigma
      return
      end
c----------------------------------------------------------
      subroutine pbdirc(phian,uld,v,pac,az,ki1,delt,dsi,i)
c     parameter (NDIM=200,NTRAJ=36000,NVECT=1500000)
      parameter (NDIM=200,NTRAJ=18000,NVECT=750000)
c     parameter (NLAT=180,NLON=360)
      parameter (NLAT=90,NLON=180)
      dimension pac(NLAT,NLON),ulo(181)
      dimension az(181,NTRAJ),dan1(NLAT,NLON),dan2(NLAT,NLON)
      common/c1/tetai(181,NTRAJ),phii(181,NTRAJ)
      common/c5/dan1,dan2
      common/c6/xo,yo,nxo,nyo,pxyo
c     calculation by integration of:
c     1/u=1/delt dsi/u(1+dan1*cos2phi+dan2*sin2phi)
c
      rad=0.017453293
      v=0.
c     print *,'pbdirc'
c     print 990,ki1,i,delt,dsi,uld,anis1
 990  format('ki1,i,delt,dsi,uld,anis',2i3,4f7.4)
      do 103 kki=1,ki1
      cs=cos(2.*az(kki,i))
      sn=sin(2.*az(kki,i))
      tet=tetai(kki,i)/rad
      ph=phii(kki,i)/rad 
c  MODIF ERIC 19 SEPT---------------------------
c     if (ph.gt.180.) ph=ph-360.
      if (ph.ge.180.) ph=ph-360.
c  END MODIF ERIC 19 SEPT---------------------------
      ct=(tet-(xo-pxyo))/pxyo
      cp=(ph-(yo-pxyo))/pxyo
c  MODIF JJ ERIC 19 SEPT 2001---------------------------
      ict=int(ct)
      if(ict.le.0)ict=1
      if(ict.gt.nxo-1)ict=nxo-1
      ict1=ict+1
      icp=int(cp)
      if(icp.le.0)icp=1
      icp1=icp+1
      if(icp.eq.nyo)icp1=1
c  END MODIF JJ ERIC 19 SEPT 2001---------------------------
c--MODIF traitemenet des poles Nord et Sud
c Eric et JJ 6/09/2001
      if((abs(tet).lt.1.e-7).or.((abs(tet-180.)).lt.1.e-7)) then
         um=pac(ict,icp)*(1.+dan1(ict,icp)*cs+dan2(ict,icp)*sn)
c        write(*,*)'POLE ',tet,um,pac(ict,icp),dan1(ict,icp),
c    *   dan2(ict,icp),cs,sn
      else
c        write(*,*)'NOPOLE '
        tei=pxyo*ict+(xo-pxyo)
        fi=pxyo*icp+(yo-pxyo)
c       if (i.eq.5)write(*,*)'TEST2 ICT ICP PXY0',ict,icp,pxyo
        d1=de(tet,ph,tei,fi)
        d2=de(tet,ph,tei+pxyo,fi)
        d3=de(tet,ph,tei,fi+pxyo)
c       write(*,*)'D3 tet ph tei fi+pxyo ',d3,tet,ph,tei,fi+pxyo
        d4=de(tet,ph,tei+pxyo,fi+pxyo)
c---MODIF
c     if(i.eq.5)write(*,*)'D1234',d1,d2,d3,d4
c       print 996,tet,ph,tei,fi
c       print 998,d1,d2,d3,d4,pac(ict,icp),pac(ict1,icp1)
        pc1=pac(ict,icp)*(1.+dan1(ict,icp)*cs+dan2(ict,icp)*sn)
        pc2=pac(ict1,icp)*(1.+dan1(ict1,icp)*cs+dan2(ict1,icp)*sn)
        pc3=pac(ict,icp1)*(1.+dan1(ict,icp1)*cs+dan2(ict,icp1)*sn)
        pc4=pac(ict1,icp1)*(1.+dan1(ict1,icp1)*cs+dan2(ict1,icp1)*sn)
c       print 996,pc1,pc2,pc3,pc4
        sd=d1*d2*d3+d2*d3*d4+d3*d4*d1+d4*d1*d2
        um=pc1*d2*d3*d4+pc2*d1*d3*d4+pc3*d1*d2*d4+pc4*d1*d2*d3
c---MODIF
c       if(i.eq.5)write(*,*)'UM,SD',um,sd
        um=um/sd
      endif
c     print 997,um
 997  format('um',f7.3)
      ulo(kki)=1./um
 103  continue
 996  format('tet,ph,tei,fi',4f7.2)
 998  format('d1,d2,d3,d4,pac1,pac4',6f7.4)
c     print 1000,(ulo(kki),kki=1,ki1)
 1000 format('ulo',10f6.3)
      call trapez(ki1,ulo,dsi,ss)
      v=-ss/delt
      v=uld+v
      return
      end
c-------------------------------------------------------------------------
      subroutine grille_vpp(vv,dsi,delt,kdi,nda,dcor,dcoran,sigma,
     1                      sigan,uld,w1,sigan1,sigan2,ut0)
c     parameter (NDIM=200,NTRAJ=36000,NVECT=1500000)
      parameter (NDIM=200,NTRAJ=18000,NVECT=750000)
c     parameter (NLAT=180,NLON=360)
      parameter (NLAT=90,NLON=180)
      parameter (NTDIM=2*NDIM+8)
      dimension vv(*),a(NLAT,NLON),pac(NLAT,NLON)
      dimension sig(NLAT,NLON),s(NTRAJ,NTRAJ),dan1(NLAT,NLON),
     *w1(5,NTRAJ)
      dimension dsi(*),kdi(*),delt(*),sia1(NTRAJ),sia2(NTRAJ),uld(NTRAJ)
      dimension fcoa2(NVECT),dan2(NLAT,NLON)
c     dimension phiao(NLAT,NLON),alpha(NLAT,NLON)
      dimension fco(NVECT),si(NTRAJ),dul(NLAT,NLON),iul(NLAT,NLON),
     .           ulf(NLAT,NLON)
      dimension az(181,NTRAJ),fcoa1(NVECT)
      dimension sigan1(NLAT,NLON),sigan2(NLAT,NLON)
      dimension aa1(NLAT,NLON),aa2(NLAT,NLON)
c     dimension da1(NLAT,NLON),da2(NLAT,NLON)
      dimension sigan2delt(NTRAJ), sigma2delt(NTRAJ)
      integer s_path(NTRAJ), e_path(NTRAJ)

      common/c1/tetai(181,NTRAJ),phii(181,NTRAJ)
      common/c2/s
      common/c3/pac,sig,ulf
      common/c4/az
      common/c5/dan1,dan2
      common/c6/xo,yo,nxo,nyo,pxyo
      common/c61/xinf,yinf,nx,ny,pxy
      common/c7/w(NTDIM,NTRAJ)

      common/dbs/theta(NVECT), phi(NVECT),s2az(NVECT),
     *c2az(NVECT),
     .                s_path, e_path

c calcul de 1/u-1/u0
      rad=0.017453293
      radi=0.5/rad
      rpxyo = 1.0/pxyo
c     print *,'grille'
      do i=1,nxo
         do j=1,nyo
            a(i,j)=0.
            ulf(i,j)=pac(i,j)
         enddo
      enddo
      it0=int((xinf-(xo-pxyo))*rpxyo)
      ip0=int((yinf-(yo-pxyo))*rpxyo)
      ipxy=int(pxy*rpxyo)
      itn=it0+(nx-1)*ipxy
      ipn=ip0+(ny-1)*ipxy

      do ij=1,nda
         sigan2delt(ij) =(sigan**2)/delt(ij)
         sigma2delt(ij) = sigma**2/delt(ij)
      enddo
      ri2dcor2 = 1.0/(2*(dcor**2))
      ri2dcoran2 = 1.0/(2*(dcoran**2))
      rmiut0 = -1/ut0
      sigma2 = sigma**2
!     write(13,*) 'In grille'

      do ii=1,nx
!         write(13,*) ii, ' of ',nx
!         call flush(13)
         tx=xinf+(ii-1)*pxy
         tt=tx*rad
         ctt = cos(tt)
         stt = sin(tt)
         it1=int((tx-(xo-pxyo))*rpxyo)
         do j=1,ny
            dul(ii,j)=0.
            an1=0.
            an2=0.
         !   iul(ii,j)=0.
            py=yinf+(j-1)*pxy
            ip1=int((py-(yo-pxyo))*rpxyo)
        !    ulf(it1,ip1)=0.
            pp=py*rad
 
            do ik = 1,e_path(nda)
               ctik=cos(theta(ik))
               stik=sin(theta(ik))
               sdmimj=ctik*ctt+stik*stt*cos(phi(ik)-pp)
!                  sdmimj=cosdel(tetai(ki,ij),phii(ki,ij),tt,pp)
               sdmimj = max( min( sdmimj, 1.0), -1.0 )
               sdmimj=acos(sdmimj)
               sdmimj=sdmimj**2
               ex=sdmimj*ri2dcor2 !/(2*(dcor**2))
               exa=sdmimj*ri2dcoran2 !/(2*(dcoran**2))
c              if(ex.lt.25.) then
               if(ex.lt.4.59) then
                  fco(ik)=exp(-ex)
c-------------------------------------------------------------------
c modif eric
c Dans les derivees partieles des parametres anisotropes c'est 1/ut0
c (avec ut0=vitesse initiale) qui intervient et non uld(=1/vitesse
c  des donnees)
c On remplace uld(ij) par 1/ut0
c     fcoa=-uld(ij)*exp(-exa)
c-------------------------------------------------------------------
                  fcoa = rmiut0*exp(-exa)
                  fcoa1(ik)=fcoa*c2az(ik)
                  fcoa2(ik)=fcoa*s2az(ik)
               else
                  fco(ik)=0.
                  fcoa1(ik)=0.
                  fcoa2(ik)=0.
               endif
            enddo
 
            do ij=1,nda
               si(ij)   = 0.5*fco(s_path(ij))
               sia1(ij) = 0.5*fcoa1(s_path(ij))
               sia2(ij) = 0.5*fcoa2(s_path(ij))
               do jk = s_path(ij)+1,e_path(ij)-1
                  si(ij)   = si(ij)   + fco(jk)
                  sia1(ij) = sia1(ij) + fcoa1(jk)
                  sia2(ij) = sia2(ij) + fcoa2(jk)
               enddo
               si(ij)   = (si(ij)   + 0.5*fco(e_path(ij)) )*dsi(ij)
               sia1(ij) = (sia1(ij) + 0.5*fcoa1(e_path(ij)))*dsi(ij)
               sia2(ij) = (sia2(ij) + 0.5*fcoa2(e_path(ij)))*dsi(ij)
 
               si(ij)   = si(ij)*sigma2delt(ij)
               sia1(ij) = sia1(ij)*sigan2delt(ij)
               sia2(ij) = sia2(ij)*sigan2delt(ij)
               an1 = an1+vv(ij)*sia1(ij)
               an2 = an2+vv(ij)*sia2(ij)
               dul(ii,j) = dul(ii,j)+vv(ij)*si(ij)
            enddo
 
!            iul(ii,j)=int(1000.*dul(ii,j))
            dan1(it1,ip1)=dan1(it1,ip1)+an1
            dan2(it1,ip1)=dan2(it1,ip1)+an2
            ulf(it1,ip1) = pac(it1,ip1)/(dul(ii,j)*pac(it1,ip1)+1.0)
!            ulf(it1,ip1)=dul(ii,j)+1./pac(it1,ip1)
!            ulf(it1,ip1)=1./ulf(it1,ip1)
       
c     alpha(it1,ip1)=sqrt((dan1(it1,ip1)**2)+(dan2(it1,ip1)**2))
c     phiao(it1,ip1)=radi*atan2(dan2(it1,ip1),dan1(it1,ip1))
 
c     calculation of the final error in it1,ip1
c
            a(it1,ip1)=0.
            aa1(it1,ip1)=0.
            aa2(it1,ip1)=0.
            do  l=1,nda
               do  k=1,nda
                  aa1(it1,ip1)=aa1(it1,ip1)+s(k,l)*sia1(k)*sia1(l)
                  aa2(it1,ip1)=aa2(it1,ip1)+s(k,l)*sia2(k)*sia2(l)
                  a(it1,ip1)=a(it1,ip1)+(s(k,l)*si(k)*si(l))
               enddo
            enddo
            siga1=sigan**2-aa1(it1,ip1)
            siga2=sigan**2-aa2(it1,ip1)
            sigan1(it1,ip1)=sqrt(abs(siga1))
            sigan2(it1,ip1)=sqrt(abs(siga2))
            sigm=sigma**2-a(it1,ip1)
            sigm=abs(sigm)
            sigm=sqrt(sigm)
            sig(it1,ip1)=sigm*(ulf(it1,ip1)**2)
         enddo
      enddo
 
c--------------------------------------------------------
c  modif eric
c     print *,'anisotropy'
c     do 120 i=it0,itn,ipxy
c     do 118 j=ip0,ipn,ipxy
c     da1(i,j)=dan1(i,j)*100.
c118  da2(i,j)=dan2(i,j)*100.
c     print 2014,(da1(i,j),j=ip0,ipn,ipxy)
c     print 2014,(da2(i,j),j=ip0,ipn,ipxy)
c120  continue
c     do 131 i=1,nxo
c131  write (9,2014)(da1(i,j),j=1,nyo)
c     do 331 i=1,nxo
c331  write (9,2014)(da2(i,j),j=1,nyo)
c     do 431 i=1,nxo
c431  write(9,2014)(sigan1(i,j),j=1,nyo)
c     do 531 i=1,nxo
c531  write (9,2014)(sigan2(i,j),j=1,nyo)
c--------------------------------------------------------
 2012 format(18f8.3)
 2013 format('kd1,delt,dcor,dcoran',i3,3f7.4)
 2020 format(30i4)
 2014 format(20f6.3)
 2021 format(30f4.2)
      return
      end
c---------------------------------------------------------
      subroutine me(x,xm,xe,n)
      dimension x(n)
      xm=0.
      xe=0.
      do 900 ii=1,n
      xm=xm+x(ii)
 900  xe=xe+x(ii)**2
      xm=xm/n
      xe=(xe/n-xm**2)
      xe=sqrt(xe)
      return
      end
c---------------------------------------------------------
      function cosdel(tetam1,phim1,tetam2,phim2)
      cosdel=0.
      c1=cos(tetam1)
      c2=cos(tetam2)
      s1=sin(tetam1)
      s2=sin(tetam2)
      cosdel=c1*c2+(s1*s2)*cos(phim1-phim2)
      return
      end
c----------------------------------------------------------
      function sindel(tetam1,phim1,tetam2,phim2)
      co=cosdel(tetam1,phim1,tetam2,phim2)
      un=1.0000000000
      a=abs(un-co**2)
      sindel=sqrt(a)
      return
      end
c------------------------------------------------------------
      function de(t1,p1,t2,p2)
      rad=0.017453293
      t=(t1+t2)*rad/2.
      de=((t2-t1)**2)+((p2-p1)**2)*(sin(t)**2)
      de=sqrt(de)
      return
      end
c---------------------------------------------------------
      subroutine ma22b(a,ia,n,w,e)
      dimension a(*),w(*)
      call a22m(a,a,ia,n,w,e)
      return
      end
c----------------------------------------------------------
      subroutine a22m(a,b,ia,n,w,e)
      dimension a(ia,*),w(*),b(*)
      l=0
      eps=0.5e-5
c     print 401
401   format('a22m')
c     call time
      do 1 j=1,n
      do 1 i=1,j
      l=l+1
      b(l)=a(i,j)
1     continue
c     print 402
402   format('sinv')
c     call time
      call sinv(a,n,eps,ierr)
c     print 402
      e=ierr
      ii=n
      jj=n
20    ii=jj
24    a(ii,jj)=b(l)
      l=l-1
      ii=ii-1
      if(ii)22,22,24
22    jj=jj-1
      if(jj)30,30,20
30    do 40  j=1,n
      do 40 i=1,j
40    a(j,i)=a(i,j)
c     print 401
c     call time
      return
      end
c----------------------------------------------------------
      subroutine sinv(a,n,eps,ier)
      real*8 din,work
      dimension a(*)


c     routine modifiee par David (2 sept 98)
c     permet d'augmenter la vitesse de 10%
c     sur le VPP


      call mfsd(a,n,eps,ier)
c     call time

      if (ier .lt. 0) return

      ipiv=n*(n+1)/2
      ind=ipiv

c      call FTRACE_REGION_BEGIN("boucle1")     
      do i=1,n
         din=1.d0/dble(a(ipiv))
         a(ipiv)=din
         min=n
         kend=i-1
         lanf=n-kend
         if(kend.gt.0) then
            j=ind
            do k=1,kend
               work=0.d0
               min=min-1
               lhor=ipiv
               lver=j
               do  l=lanf,min
                  lver=lver+1
                  lhor=lhor+l
                  work=work+dble(a(lver)*a(lhor))
               enddo
               a(j)=-work*din
               j=j-min
            enddo
         endif
        ipiv=ipiv-min
         ind=ind-1
      enddo
c      CALL FTRACE_REGION_END("boucle1")        
c     call time

c      call FTRACE_REGION_BEGIN("boucle2")    
      do  i=1,n
         ipiv=ipiv+i
         j=ipiv
         do k=i,n
            work=0.d0
!DBS            lhor=j
            do l = 0,n-k  !DBS l=k,n
               lhor = j + l*k + l*(l-1)/2
               lver=lhor+k-i
               work=work+dble(a(lhor)*a(lver))
!DBS               lhor=lhor+l
            enddo
            a(j)=work
            j=j+k
         enddo
      enddo
c      CALL FTRACE_REGION_END("boucle2")   

      return
      end
c---------------------------------------------------
      subroutine mfsd(a,n,eps,ier)
      real*8 dpiv,dsum
      dimension a(*)
      if(n-1) 12,1,1
1     ier=0
      kpiv=0
      do 11 k=1,n
      kpiv=kpiv+k
      ind=kpiv
      lend=k-1
      tol=abs(eps*a(kpiv))
      do 11 i=k,n
      dsum=0.d0
      if(lend) 2,4,2
2     do 3 l=1,lend
      lanf=kpiv-l
      lind=ind-l
3     dsum=dsum+dble(a(lanf)*a(lind))
4     dsum=dble(a(ind))-dsum
      if(i-k) 10,5,10
5     if(sngl(dsum)-tol) 6,6,9
6     if(dsum) 12,12,7
7     if(ier) 8,8,9
8     ier=k-1
9     dpiv=dsqrt(dsum)
      a(kpiv)=dpiv
      dpiv=1.d0/dpiv
      goto 11
10    a(ind)=dsum*dpiv
11    ind=ind+i
      if(ier.ne.0) print 21,ier
      return
12    ier=-1
      print 22
21    format (' ret',i4)
22    format (' break')
      return
      end
c------------------------------------------------
      subroutine khi(x,ed,qui,n)
      dimension x(n),ed(n)
      qui=0.
      do 1 i=1,n
      xed=x(i)/ed(i)
      qui=qui+xed**2
   1  continue
      qui=qui/float(n)
      qui=sqrt(qui)
      return
      end
c-----------------------------------------------------------------
