c  Pgm "continuous regionalization"  of Montagner(1986, annales
c  geophysicae, 4: b3, 283-294) optimized by
c  E. Debayle and M. Sambridge
c  for application to large datasets (up to 100 000 paths)
c  The storage of data arrays is optimized, the computation 
c  CmGt and GmGt has been entirely rewritten (see routines
c  grillestack, CmGtstack and GCmGtstack) to speed up the 
c  computation and (GCmGt +Cd)-1 is not inverted directly
c  but instead, the system (GCmGt +Cd)-1*(d-Gmo) is solved
c  using a conjugate gradient method for sparse system
c (see numerical recipies pp 77)
c This version is the same as regiostack5 excepted that
c problems at the north and south poles have been removed
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  other modifs eric.
c  modifie pour sortir les parametres dan1 et dan2 de la 
c  subroutine grillestack(eric avril 95).
c 
c  Les derivees partielles des parametres anisotropes sont corrigees.
c  ON SORT Viso A1 A2 (Smith et Dahlen) de ce pgm alors que 
c  les parametres
c  inverses sont les (-A*vitesse a priori/vitesse inversee**2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        parameter (NTRAJ=10000,NSTEP=181,NDEPTH=50)
c       parameter (NLAT=360,NLON=720)
        parameter (NLAT=90,NLON=180)
        parameter (NTDEPTH=2*NDEPTH+8,KMAX=NLAT*NLON)
        parameter (NTR=NTRAJ/2)
        parameter (NMAXS=NTRAJ*3500)

      DOUBLE PRECISION v(NTRAJ),vv(NTRAJ), uld(NTRAJ),res(NTRAJ)
      dimension ux(NTRAJ),delt(NTRAJ),dsi(NTRAJ),kdi(NTRAJ)
      dimension ed(NTRAJ)
      dimension ed1(NTRAJ),resn(NTRAJ),cov(NTRAJ)
      character*4 staz1(NTRAJ),staz21(NTRAJ)
      character*8 staz2(NTRAJ),staz22(NTRAJ)

      dimension iray(KMAX,NTR),iss(KMAX,NTR),ise(KMAX,NTR)

      dimension dul(KMAX),an1(KMAX),an2(KMAX)
      dimension teta(KMAX),phi(KMAX)
      dimension nray(KMAX)
      logical mask(KMAX)

      DOUBLE PRECISION sa,err,tol,re

      dimension w1(5,NTRAJ)
      dimension az(NSTEP,NTRAJ)

      dimension t(NDEPTH),ndat(NDEPTH),jt(NDEPTH)

      dimension dan1(NLAT,NLON),dan2(NLAT,NLON)
      dimension pac(NLAT,NLON),ulf(NLAT,NLON)

      character*32 fileout,fianiout,ficdata
      common/c1/tetai(NSTEP,NTRAJ),phii(NSTEP,NTRAJ)
      common/c3/pac,ulf
      common/c4/az
      common/c5/dan1,dan2
      common/c6/xo,yo,nxo,nyo,pxyo
      common/c61/xinf,yinf,nx,ny,pxy
      common/c7/w(NTDEPTH,NTRAJ)
      common /mat/sa(NMAXS),ija(NMAXS)

      pi=3.141592653
      rad=pi/180.
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
      do 9 j=1,NTDEPTH
   9  w(j,i)=10.
c
c   reading of dataset in a formatted file data
c   After data are read they are ordered in an array w(NTDEPTH,ntra)
c   w(1,i),w(2,i) : m1(teta1,phi1)   w(3,i),w(4,i):m2(teta2,phi2)
c                    (epicenter)                    (station)
c   w(5,i): phipl-phi(average azse) azimuth with respect to the direction 
c   of the plate velocity (not used in this version)
c   w(7....NDEPTH+6,i)  seismic velocities measured at different periods or depths.
c   w(NDEPTH+7...NTDEPTH,i)  errors on seismic velocities.
c   delt epicentral distance for each path i, dsi incremental step on
c   each path i. kdi number of sampled points along each path i.
c   tetai(colat(mi),i) phii(long(mi),i)
c
      print *,'enter data file name'
      read(5,'(a32)') ficdata
      print 999,ficdata
      open(10,file=ficdata,form='formatted')
       call lecdat(nper,staz1,staz2,delt,kdi,dsi,ntra,d1,t)

c  In lecdat, teta must be between 0. and teta and phi(radians)
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

       k1=1
c      loop on periods
       do 5 ktt=1,nbper

       ntt=jt(ktt)
       k2=ntt
c  in vinit,one reads starting velocity and anisotropy
c  at a given period, data are sorted because the number
c  of data can be period dependent
      print *,' '
      do 55 jk=k1,k2
        call vinit(sigma,sigan,ut0,tt)
 55   k1=k2+1
      write(7,*)tt
      write(9,*)tt
      k=0
      do 6 j=1,ntra
      nk=ntt+6
      nkk=ntt+NDEPTH+6
      if(w(nkk,j).gt.2.) go to 6
      k=k+1
      uld(k)=w(nk,j)
      ed(k)=w(nkk,j)
      ed1(k)=ed(k)/uld(k)**2
      kdi(k)=kdi(j)
      dsi(k)=dsi(j)
      delt(k)=delt(j)
      staz21(k)=staz1(j)
      staz22(k)=staz2(j)
      kk1=kdi(j)
      do 7 ki=1,kk1
        az(ki,k)=az(ki,j)
        tetai(ki,k)=tetai(ki,j)
  7   phii(ki,k)=phii(ki,j)
      w1(1,k)=w(1,j)
      w1(2,k)=w(2,j)
      w1(3,k)=w(3,j)
      w1(4,k)=w(4,j)
      w1(5,k)=-uld(k)*cos(2.*w(5,j))
  6   continue
      ndat(ntt)=k
      print *,' '
      print 2007,t(ntt),k,dcor
 2007 format('periode,nda,dcor',f6.2,i6,f6.2)
c-------------------------------------------------------
c modif eric
c on calcul les moyennes et ecarts-types en lenteur
c car on inverse en lenteur.
      print *,'selected data '
      do 4 it=1,k
  4   ux(it)=1./uld(it)
      print 2003,(ux(i),i=1,k)
      call me(uld,xm,xe,k)
c-------------------------------------------------------
      write(*,*) 'moyenne et ecart-type de uld (1/donnees)',xm,xe
      varmoy=xe*xe
      write(*,*)'variance-donnees p.r. a la moyenne des mesures:',varmoy
      write(*,*) 'modele initial= ',ut0,
     &   '           moyenne des mesures= ',1/xm
      x0=xe
 1100 format(' initial values : xm= ',f7.3,' xe= ',f7.4)
      nda=ndat(ntt)

c---------------------------------------------------------------
cComputation of (GCmGt +Cdo) (E. Debayle march 2002)

      do 50 j=1,KMAX
            mask(j)=.TRUE.
            phi(j)=0.
            teta(j)=0.
            nray(j)=0.
c        do 50 jj=1,NTR
c           iss(j,jj)=0.
c           ise(j,jj)=0.
c           iray(j,jj)=0.
 50   continue
 
      do 54 ii=1,nx
         do 54 jj=1,ny
          j2=ny*(ii-1)+jj
          teta(j2)=xinf+(ii-1)*pxy+pxy/2.
          teta(j2)=teta(j2)*rad
          phi(j2)=yinf+(jj-1)*pxy+pxy/2.
          phi(j2)=phi(j2)*rad
 54   continue

      jold=0
      do 70 i=1,nda
        cov(i)=0
        tetai(1,i)=w1(1,i)
        phii(1,i)=w1(2,i)
        kki=kdi(i)

        do 70 ki=1,kki

           tx=tetai(ki,i)/rad
           py=phii(ki,i)/rad
c modif 29/07/2002 :
c En tomo regionale lorsque la region
c chevauche la longitude +180 degres
c (i.e. Australie)
c il faut remettre py entre 0 et 360
c            if(py.lt.0) py=py+360.

           it2=int((tx-(xinf-pxy))/pxy)
           ip2=int((py-(yinf-pxy))/pxy)

           j2=ny*(it2-1)+ip2

           if(i.eq.nda.and.ki.eq.kki.and.j2.eq.jold)ise(j2,nray(j2))=ki

           if(j2.ne.jold.or.ki.eq.1) then
             nray(j2)=nray(j2)+1
             if(nray(j2).gt.NTR)then
                 write(*,*)'J2 ',j2,' nray(j2) ',nray(j2)
                 write(*,*)'NTR ',NTR
                 stop' NTR must be increased'
             endif
             iray(j2,nray(j2))=i
             iss(j2,nray(j2))=ki
             if(i.eq.nda.and.ki.eq.kki)then
                 write(*,*) 'DERNIER POINT DERNIER TRAJET'
                 ise(j2,nray(j2))=ki
             endif
             if(ki.ne.1)then
                  ise(jold,nray(jold))=ki-1
             else if(ki.eq.1.and.jold.ne.0.and.j2.ne.jold)then
                  ise(jold,nray(jold))=kdi(i-1)
             else if(ki.eq.1.and.jold.ne.0.and.j2.eq.jold)then
                  ise(jold,nray(jold)-1)=kdi(i-1)
             endif
             jold=j2
           endif

 70    continue

      fpxys2=pxy/2.
      nxny=nx*ny
      nxm=nx-1
      nym=ny-1
      nyp=ny+1
c     nypxy=ny*int(pxy)
      nypxy=ny*pxy
      fnys2=float(ny)/2.
      nynxm=ny*nxm
      pimpxy=180.-pxy
      dcorlim=2.648*dcor
      dcor2=dcor**2
      dcoran2=dcoran**2
      sigma2=sigma**2
      sigan2=sigan**2
      fcoagcgt=(1/ut0)**2*sigan2
      fcoacgt=-(1/ut0)*sigan2

       k=nda+1
       ija(1)=nda+2
       do 100 i=1,nda
         tetai(1,i)=w1(1,i)
         phii(1,i)=w1(2,i)
         kki=kdi(i)
         do 110 ki=1,kki
         call GCmGtstack(iray,iss,ise,cov,dsi,kdi,teta,phi,nray,mask
     *,dcorlim,dcor2,dcoran2,sigma2,fcoagcgt,nxny,nym,nyp,nypxy
     *,fnys2,nynxm,pimpxy,fpxys2,i,ki,kki)
  110 continue

c we store the matrix (GCmGt+Cdo) in row indexed sparse
c storage mode (see num. recipies p 71)
c we also reconstruct the full matrix (which is symetric 
c from its upper triangle part) (when if(ij.lt.i)  test is verified)

       do 105 ij=1,nda
          if(ij.lt.i) then
             do it=ija(ij),ija(ij+1)-1
               if(ija(it).eq.i)then
                 k=k+1  
                 sa(k)=sa(it)
                 ija(k)=ij
               endif
             enddo
          else if(ij.eq.i) then
             sa(ij)=ed(i)**2+cov(i)/(delt(i)*delt(ij))
          else if(abs(cov(ij)).gt.10e-10)then
             k=k+1
             sa(k)=cov(ij)/(delt(i)*delt(ij))
             ija(k)=ij
          endif
          cov(ij)=0
 105   continue
       ija(i+1)=k+1
       if((k+1).gt.NMAXS) pause' NMAXS too small in main code'
 100   continue
       write (*,*)'KSAMAX ',k+1,' NMAXS ',NMAXS

c      write(*,*)' ija(k) '
c      write(*,'(10(1x,i6))') (ija(ii),ii=1,k+1)
c      write(*,*) 'sa(k) '
c      write(*,*) (sa(ii),ii=1,k+1)

c--------------------------------------------------------------
c   forward problem resolution.
c   Vector v contains (d-g(p))in slowness

      do 130 i=1,nda
      v(i)=0.
      vv(i)=0.
      call pbdirc(uld(i),v(i),pac,az,kdi(i),delt(i),dsi(i),i)
 130  continue
      print*,'vecteur v (d-g(p) en lenteur):'
      print 2003,(v(i),i=1,nda)
cJJL  Calcul de la variance des donnees par rapport aux synthetiques
cJJL  du modele initial (et non par rapport a la moyenne des donnees)
      varini=0.
      do 139 i=1,nda
      varini=varini+v(i)*v(i)
139   continue
      varini=varini/nda
      write(*,*)'variance-donnees p.r. au modele initial:       ',varini
c  modif eric : calcul du qui initial en lenteur 
      call khi (v,ed,qui,nda)
      print*, 'qui initial :',qui

c  solving the linear system (GCmGt +Cdo)-1*(d-Gmo) using
c  conjugate gradient

      itol=3
      itmax=1000
      tol= 1e-12
      iter=0
      err=0

      write(*,*)' INPUT LINBCG '
      write(*,*) 'NDA ',nda,' ITOL ',itol,' TOL ',tol,' ITMAX ',itmax
      write(*,*) 'ITER ',iter,' ERR ',err
c     print *,'vv(i)'
c     write(*,*) (vv(i),i=0,nda)
      call linbcg(nda,v,vv,itol,tol,itmax,iter,err)

 2012 format(10(f9.6,1x))

      call me(vv,xm,xe,nda)
      print 1100,xm,xe
c###################################################################### 
c  This part of the code include the previous grillestack subroutine
c  which call the CmGtstack2 code
c NB :
c  repere 1 : xo,yo,nxo,nyo,pxyo
c  repere 2 : xinf,yinf,nx,ny,pxy

      do 9950 j=1,KMAX
            mask(j)=.TRUE.
            dul(j)=0.
            an1(j)=0.
            an2(j)=0.
9950  continue

      do 99110 ij=1,nda
        kd1=kdi(ij)
        do 99120 ki=1,kd1
        call CmGtstack2(dul,an1,an2,vv(ij),dsi(ij),teta,phi,mask,
     *dcorlim,dcor2,dcoran2,sigma2,fcoacgt,nxny,nym,nyp,nypxy,
     *fnys2,fpxys2,nynxm,pimpxy,delt(ij),ij,ki,kd1)

99120   continue
99110 continue

      do 99140 ii=1,nx
        tx=xinf+(ii-1)*pxy+fpxys2
        it1=int((tx-(xo-pxyo))/pxyo)

        do 99140 jj=1,ny
          py=yinf+(jj-1)*pxy+fpxys2
          ip1=int((py-(yo-pxyo))/pxyo)
          if(it1.gt.NLAT.or.it1.le.0.or.ip1.gt.NLON.or.ip1.le.0)then
              write(*,*)'IT1 ',it1,'IP1',ip1
              stop 'pb dimension'
          endif
          ulf(it1,ip1)=1./pac(it1,ip1)
          j2=ny*(ii-1)+jj

          dan1(it1,ip1)=dan1(it1,ip1)+an1(j2)
          dan2(it1,ip1)=dan2(it1,ip1)+an2(j2)
          ulf(it1,ip1)=ulf(it1,ip1)+dul(j2)
          ulf(it1,ip1)=1./ulf(it1,ip1)

c     The calculation of the final a posteriori covariance
c     matrix is skipped as it would require to solve m times
c     (m beeing the number of model parameter) the
c     system (GCmGt +Cdo)-1*(each column of GCm)
c     using the conjugate gradient method

99140  continue
 
c###################################################################### 
      do 94 i=1,nda
      res(i)=0.
      re=res(i)
      call pbdirc(uld(i),re,ulf,az,kdi(i),delt(i),dsi(i),i)
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
      icomp=0
c------------------------------------------------------------------
c--x0 et xe ont ete passe en lenteur, normalement redvar esten 1/Q.

      call me(res,xm,xe,nda)

      write(*,*) 'moyenne et ecart-type de res (residus fin)',xm,xe
      varres=xm*xm+xe*xe
      write(*,*)'variance-donnees p.r. a la moyenne des mesures:',varmoy
      write(*,*)'variance-donnees p.r. au modele final:         ',varres
      redvar =(varmoy-varres)/varmoy
      redvar2=(varini-varres)/varini

      print *,'red. variance(%val moy des donnees en lenteur)',redvar
      print *,'red. variance(% modele initial)     en lenteur',redvar2
c------------------------------------------------------------------
c modif eric : ed est l'erreur en 1/V, ed1 celle en vitesse
c On doit donc calculer un qui en 1/V
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
 126  write(7,2014)(sigma,j=1,nyo)
c-------------------------------------------------------------------
c modif eric Dec 1995
c Les derivees partieles des parametres anisotropes sont (-cos2teta)/L*ut0
c avec  ut0=vitesse initiale(cf grillestack). Ceci implique(voir these)que les 
c parametres inverses sont dan=(A1*ut0)/(ulf*ulf) ulf est la vitesse inversee
c et non la lenteur car ulf est converti en vitesse dans la routine grillestack.
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
 123  write(9,2014)(sigan,j=1,nyo)
      do 122 i=1,nxo
 122  write(9,2014)(sigan2,j=1,nyo)
c------------------------------------------------------------------
  5   continue

 1007 format(6(a4,a8,1x,f5.2,1x))
 2003 format(10f7.4)
 2009 format(20f6.4)
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
c-modif eric 27/10/2002--------------------------------
      if (abs(cosdel(tetam1,phim1,tetam2,phim2)).gt.1)
     *stop 'coordo : abs(cosdel)>1'
      if (abs(sindel(tetam1,phim1,tetam2,phim2)).gt.1.or.
     *sindel(tetam1,phim1,tetam2,phim2).eq.0)
     *stop 'coordo :  sindel=0 or abs(sindel)>1'
c------------------------------------------------------
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
      subroutine lecdat(nper,staz1,staz2,delt,kdi,dsi,ndat,d1,t)
      
        parameter (NTRAJ=10000,NSTEP=181,NDEPTH=50)
        parameter (NSTAMAX=400)
c       parameter (NLAT=360,NLON=720)
        parameter (NLAT=90,NLON=180)
        parameter (NTDEPTH=2*NDEPTH+8,KMAX=NLAT*NLON)

      dimension delt(NTRAJ),dsi(NTRAJ),kdi(NTRAJ)
      dimension wlat(NSTAMAX),wlon(NSTAMAX)
      dimension t(NDEPTH),ug(NDEPTH),edat(NDEPTH)
c     dimension bcd(NTRAJ,2),bcdb(NSTAMAX)
      dimension az(NSTEP,NTRAJ)
      dimension deltei(NSTEP)
      character staz*13
      character*4 nomsta(NSTAMAX),staz1(NTRAJ),nom1,nom2
      character*8 staz2(NTRAJ)

      common/c1/tetai(NSTEP,NTRAJ),phii(NSTEP,NTRAJ)
      common/c4/az
      common/c7/w(NTDEPTH,NTRAJ)

c   in this subroutine, we read the dataset(phase or group velocities)
c   the coordinates of the points along the paths are computed.

      if(nper.gt.NDEPTH) stop 'trop de periodes lues'
      rad=0.017453293
      print *,'nsta: number of stations'
      read(10,*)nsta
      print *,nsta

      do 88 i=1,nsta
         read(10,1998)nomsta(i),wlat(i),wlon(i)
c        read(10,1998)bcdb(i),wlat(i),wlon(i)
c        MODIF eric pour avoir des lon entre 0 et 360 degres
c        if(wlon(i).lt.0.)wlon(i)=wlon(i)+360
         print 1998,nomsta(i),wlat(i),wlon(i)
  88  continue
 1998 format(a4,2x,2f9.4)

      nt=nper
      do 1 i=1,ndat
c       read(10,2000)(bcd(i,j),j=1,2)
c       print 2000,(bcd(i,j),j=1,2)
c2000   format(a4,a6)
        read(10,'(a13)')staz
        write(*,'(a13)')staz

        ir=index(staz,'.')
        staz1(i)=staz(:ir-1)
        staz2(i)=staz(ir+1:)
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
  89    continue
        if(icheck.ne.1) stop' station non trouvee '
        read(10,*)alate,alone
c       write(*,*)'ICI',alate,alone,w(1,i),w(2,i)
c       write(*,*)'Lat lon epicenter, station',alate,alone,w(1,i),w(2,i)
c       MODIF eric pour avoir des lon entre 0 et 360 degres
c       if(alone.lt.0.)alone=alone+360
c       print 2002,alate,alone
        w(1,i)=(90.-w(1,i))*rad
        w(2,i)=w(2,i)*rad
        in0=1
        fi=45.
        w(3,i)=(90.-alate)*rad
        w(4,i)=alone*rad
        w(5,i)=fi*rad
        w(6,i)=nt
        read(10,*)(ug(j),j=1,nt)
c       ug : group or phase velocity data for this path
        print 1999,(ug(j),j=1,nt)
c       edat error on ug
 1999   format(11f7.3)
        read(10,*)(edat(j),j=1,nt)
        do 2 k=1,nt
          kk=k+6
          w(kk,i)=1/ug(k)
          kkk=k+NDEPTH+6
          w(kkk,i)=edat(k)/(ug(k)**2)
  2     continue
        in=1
c2002   format(2f8.3,3f7.2)
c       print 2008,(w(k,i),k=7,NDEPTH+6)
 2008   format('1/ug',13f6.3)
c
c       calculation of coordinates between m1 and m2
c       one computes the number of steps between m1i and m2i and dsi
c       is recalculated such that dsi*kdi=delt(i)
c
c modif eric 27/10/2002---------------------------
        delt(i)=cosdel(w(1,i),w(2,i),w(3,i),w(4,i))
        if (abs(delt(i)).gt.1)stop'abs(delt(i))>1'
c-------------------------------------------------
        delt(i)=acos(delt(i))
c       delt(i)=acos(cosdel(w(1,i),w(2,i),w(3,i),w(4,i)))
        kdi(i)=int(delt(i)/d1)
        dsi(i)=delt(i)/kdi(i)
        ki=1
        kdi(i)=kdi(i)+1
        kd1=kdi(i)
        if(kd1.gt.NSTEP) then
          write(*,*)'Path ',i,'kd1= ',kd1,' NSTEP= ',NSTEP
          stop'kd1 greater than NSTEP'
        endif
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
     *phii(ki,i),deltei(ki),az(ki,i))
  4       continue
  1     continue
      return
      end
c----------------------------------------------------------
      subroutine vinit(sigma,sigan,ut0,tt)

c       parameter (NLAT=360,NLON=720)
        parameter (NLAT=90,NLON=180)

      dimension pac(NLAT,NLON)
c     dimension pac(NLAT,NLON),sig(NLAT,NLON)
      dimension ulf(NLAT,NLON),dan1(NLAT,NLON),dan2(NLAT,NLON)

      common/c3/pac,ulf
c     common/c3/pac,sig,ulf
      common/c5/dan1,dan2
      common/c6/xo,yo,nxo,nyo,pxyo
c     sigma:erreur a priori sur les parametres
      print *,'enter tt,ut0,sigma,anis1,anis2,sigan'
      read(*,*)tt,ut0,sigma,anis1,anis2,sigan
      print *,tt,ut0,sigma,anis1,anis2,sigan
      do 3 i=1,nxo
      do 3 j=1,nyo
      pac(i,j)=ut0
c     sig(i,j)=sigma
      dan1(i,j)=anis1
  3   dan2(i,j)=anis2
      sigma=sigma/(ut0**2)
      print*,'sigma'  
      write(*,*) sigma
      return
      end
c----------------------------------------------------------
      subroutine pbdirc(uld,v,pac,az,ki1,delt,dsi,i)
c     subroutine pbdirc(uld,v,ki1,delt,dsi,i)

        parameter (NTRAJ=10000,NSTEP=181)
c       parameter (NLAT=360,NLON=720)
        parameter (NLAT=90,NLON=180)
        parameter (KMAX=NLAT*NLON)

      dimension pac(NLAT,NLON),dan1(NLAT,NLON),dan2(NLAT,NLON)
c     dimension sig(NLAT,NLON),ulf(NLAT,NLON)

      dimension ulo(NSTEP)

      dimension az(NSTEP,NTRAJ)

      DOUBLE PRECISION v,uld

      common/c1/tetai(NSTEP,NTRAJ),phii(NSTEP,NTRAJ)
c     common/c3/pac,sig,ulf
c     common/c4/az
      common/c5/dan1,dan2
      common/c6/xo,yo,nxo,nyo,pxyo
c     calculation by integration of:
c     1/u=1/delt dsi/u(1+dan1*cos2phi+dan2*sin2phi)
c
      rad=0.017453293
      v=0.
c     print *,'pbdirc'
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
      subroutine me(x,xm,xe,n)
c     dimension x(n)
      DOUBLE PRECISION x(n)
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
      subroutine khi(x,ed,qui,n)
c     dimension x(n),ed(n)
      dimension ed(n)
      DOUBLE PRECISION x(n)
      qui=0.
      do 1 i=1,n
      xed=x(i)/ed(i)
      qui=qui+xed**2
   1  continue
      qui=qui/float(n)
      qui=sqrt(qui)
      return
      end
c--------------------------------------------------------------
      SUBROUTINE CmGtstack2(dul,an1,an2,vv,dsi,teta,phi,mask,
     *dcorlim,dcor2,dcoran2,sigma2,fcoacgt,nxny,nym,nyp,nypxy,
     *fnys2,fpxys2,nynxm,pimpxy,delt,ij,ki,kd1)

        parameter (NTRAJ=10000,NSTEP=181)
c       parameter (NLAT=360,NLON=720)
        parameter (NLAT=90,NLON=180)
        parameter (KMAX=NLAT*NLON)
 
        DOUBLE PRECISION vv

        dimension dul(KMAX),an1(KMAX),an2(KMAX)
        dimension phi(KMAX),teta(KMAX)

        dimension az(NSTEP,NTRAJ)

        dimension indice(8),list(KMAX)
        logical   mask(KMAX)

        common/c1/tetai(NSTEP,NTRAJ),phii(NSTEP,NTRAJ)
        common/c4/az
        common/c6/xo,yo,nxo,nyo,pxyo
        common/c61/xinf,yinf,nx,ny,pxy

        rad=0.017453293
 
c       mask(j) =TRUE => cell j has not been put on stack

c indices are in repere 2 : xinf,yinf,nx,ny,pxy
 
             tx=tetai(ki,ij)/rad
             py=phii(ki,ij)/rad
c En tomo regionale lorsque la region
c chevauche la longitude +180 degres
C (i.e. Australie)
c il faut remettre py entre 0 et 360
c            if(py.lt.0) py=py+360.
 
             it2=int((tx-(xinf-pxy))/pxy)
             ip2=int((py-(yinf-pxy))/pxy)
 
             j2=ny*(it2-1)+ip2

             kk=1
             list(kk)=j2
             mask(j2)=.FALSE.

             call stackinit()
 
100          sdmimj=cosdel(tetai(ki,ij),phii(ki,ij),
     *teta(j2),phi(j2))
             if(mod(j2,ny).eq.0)then
               it2=int(j2/ny)
             else
               it2=int(j2/ny)+1
             endif
             tx2=xinf+(it2-1)*pxy+fpxys2

             if(sdmimj.gt.1) sdmimj=1
             if(sdmimj.lt.-1)sdmimj=-1
             sdmimj=acos(sdmimj)

             if(sdmimj.lt.dcorlim) then
               sdmimj2=sdmimj**2
               fco=sigma2*exp(-sdmimj2/(2*dcor2))
               fcoa=fcoacgt*exp(-sdmimj2/(2*dcoran2))
               fcoa1=fcoa*cos(2.*az(ki,ij))
               fcoa2=fcoa*sin(2.*az(ki,ij))

c when the contribution of the first and last point of
c each great circle is not zero we take half of it.
c This should make the summation equivalent to a
c trapezoidal integration. TO BE CHECKED ON REGIOSTACK3

                if(ki.eq.1.or.ki.eq.kd1) then
                 dul(j2)=dul(j2)+(0.5*fco*dsi/delt)*vv
                 an1(j2)=an1(j2)+(0.5*fcoa1*dsi/delt)*vv
                 an2(j2)=an2(j2)+(0.5*fcoa2*dsi/delt)*vv
                else
                 dul(j2)=dul(j2)+(fco*dsi/delt)*vv
                 an1(j2)=an1(j2)+(fcoa1*dsi/delt)*vv
                 an2(j2)=an2(j2)+(fcoa2*dsi/delt)*vv
                endif
 
c finding the neighbour of j2

                indice(1)=j2-nyp
                indice(2)=j2-ny
                indice(3)=j2-nym
                indice(4)=j2-1
                indice(5)=j2+1
                indice(6)=j2+nym
                indice(7)=j2+ny
                indice(8)=j2+nyp
c eastern boundary
             if (mod(j2,ny).eq.0.and.nypxy.eq.360.) then
                indice(3)=j2-2*ny+1
                indice(5)=j2-nym
                indice(8)=j2+1
             elseif(mod(j2,ny).eq.0.and.nypxy.ne.360.) then
                write(*,*)' ny*pxy =',nypxy
                stop'agrandir la zone a l est'
             endif
c western boundary
             if (mod(j2,ny).eq.1.and.nypxy.eq.360.) then
                indice(1)=j2-1
                indice(4)=j2+nym
                indice(6)=j2+2*ny-1
             elseif(mod(j2,ny).eq.1.and.nypxy.ne.360.) then
                write(*,*)' ny*pxy =',nypxy
                stop'agrandir la zone a l ouest'
             endif
c northern boundary and north pole
             if(j2.le.ny.and.tx2.le.pxy) then
                if(nypxy.ne.360.) then
                    write(*,*)' 2 ny*pxy =',nypxy
                    stop'Path at north pole requires to
     *     define a region covering the 360 deg. of longitude variation'
                endif
                if(mod(ny,2).ne.0) stop'uneven number of column?'
                if(j2.gt.int(fnys2)) then
                  indice(2)=j2-int(fnys2)
                  if(indice(2)-1.gt.0) then
                     indice(1)=indice(2)-1
                  else
                   indice(1)=ny
                  endif
                  indice(3)=indice(2)+1
                endif
                if(j2.le.int(fnys2)) then
                  indice(2)=j2+int(fnys2)
                  if(indice(2)+1.gt.ny) then
                     indice(3)=1
                  else
                     indice(3)=indice(2)+1
                  endif
                  indice(1)=indice(2)-1
                endif
             else if (j2.le.ny.and.tx2.gt.pxy) then
                stop'agrandir la zone au Nord'
             endif
c southern boundary and south pole
             if(j2.gt.nynxm.and.tx2.ge.pimpxy) then
                if(nypxy.ne.360.) then
                    write(*,*)' ny*pxy =',nypxy
                    stop'Path at south pole requires to
     *     define a region covering the 360 deg. of longitude variation'
                endif
                if(mod(ny,2).ne.0) stop'uneven number of column?'
                if(j2.gt.int(fnys2)) then
                  indice(7)=j2-int(fnys2)
                  if(indice(7)-1.gt.nynxm) then
                     indice(6)=indice(7)-1
                  else
                     indice(6)=nxny
                  endif
                  indice(8)=indice(7)+1
                endif
                if(j2.le.int(fnys2)) then
                  indice(7)=j2+int(fnys2)
                  if(indice(7)+1.gt.nxny) then
c                   indice(8)=1
                    indice(8)=nynxm+1
                  else
                    indice(8)=indice(7)+1
                  endif
                  indice(6)=indice(7)-1
                endif
             else if(j2.gt.nynxm.and.tx2.lt.pimpxy) then
                stop 'Agrandir la zone au Sud'
             endif

c------------------------------------------- 
 
                do inb=1,8
                 jcurrent=indice(inb)
                  if(jcurrent.gt.KMAX.or.jcurrent.le.0) 
     *stop 'Dimension array exceeded'
                  if(mask(jcurrent)) then
                     call push(jcurrent)
                     kk=kk+1
                     if(kk.gt.KMAX) stop 'Dimension array list exceeded'
                     list(kk)=jcurrent
                     mask(jcurrent)=.FALSE.
                   endif
                end do
             endif
 
             call stackempty(ibid)
             if(ibid.eq.1) go to 200
 
             call pop(jj)
             j2=jj
 
             goto 100
 
200          call stackflush()
        do k=1,kk
          mask(list(k))=.TRUE.
        end do
       return
       end
c------------------------------------------------
      SUBROUTINE GCmGtstack(iray,iss,ise,cov,dsi,kdi,teta,phi,nray,mask
     *,dcorlim,dcor2,dcoran2,sigma2,fcoagcgt,nxny,nym,nyp
     *,nypxy,fnys2,nynxm,pimpxy,fpxys2,ij,ki,kki)
        parameter (NTRAJ=10000,NSTEP=181)
c       parameter (NLAT=360,NLON=720)
        parameter (NLAT=90,NLON=180)
        parameter (KMAX=NLAT*NLON)
        parameter (NTR=NTRAJ/2)
 
        dimension  dsi(NTRAJ),kdi(NTRAJ),cov(NTRAJ)
 
        dimension phi(KMAX),teta(KMAX)
        dimension nray(KMAX)
 
        dimension iray(KMAX,NTR),iss(KMAX,NTR),ise(KMAX,NTR)

        dimension az(NSTEP,NTRAJ)
 
        dimension indice(8),list(KMAX)
        logical   mask(KMAX)
 
 
        common/c1/tetai(NSTEP,NTRAJ),phii(NSTEP,NTRAJ)
        common/c4/az
        common/c6/xo,yo,nxo,nyo,pxyo
        common/c61/xinf,yinf,nx,ny,pxy
 
        rad=0.017453293
 
c       mask(j) =TRUE => cell j has not been put on stack
 
c indices are in repere 2 : xinf,yinf,nx,ny,pxy
 
             tx=tetai(ki,ij)/rad
             py=phii(ki,ij)/rad
c En tomo regionale lorsque la region
c chevauche la longitude +180 degres
C (i.e. Australie)
c il faut remettre py entre 0 et 360
c            if(py.lt.0) py=py+360.
 
             it2=int((tx-(xinf-pxy))/pxy)
             ip2=int((py-(yinf-pxy))/pxy)
 
             j2=ny*(it2-1)+ip2
 
             kk=1
             list(kk)=j2
             mask(j2)=.FALSE.
 
             call stackinit()
 
100          sdmimj=cosdel(tetai(ki,ij),phii(ki,ij),
     *teta(j2),phi(j2))
             if(sdmimj.gt.1) sdmimj=1
             if(sdmimj.lt.-1) sdmimj=-1
             sdmimj=acos(sdmimj)

             if(mod(j2,ny).eq.0)then
               it2=int(j2/ny)
             else
               it2=int(j2/ny)+1
             endif
             tx2=xinf+(it2-1)*pxy+fpxys2
 
             if(sdmimj.lt.dcorlim) then

               do 110 ir=1,nray(j2)

             if(iray(j2,ir).ge.ij) then

                 do 112 ik=iss(j2,ir),ise(j2,ir)

                  sdmimj=cosdel(tetai(ki,ij),phii(ki,ij),
     *tetai(ik,iray(j2,ir)),phii(ik,iray(j2,ir)))
                  if(sdmimj.gt.1)sdmimj=1
                  if(sdmimj.lt.-1)sdmimj=-1
                  sdmimj=acos(sdmimj)

                  if(sdmimj.lt.dcorlim) then

                    sdmimj2=sdmimj**2
                    fco=sigma2*exp(-sdmimj2/(2*dcor2))
                    fcoa=fcoagcgt*exp(-sdmimj2/(2*dcoran2))
                    fcoa1=fcoa*cos(2.*az(ik,iray(j2,ir)))*
     *cos(2.*az(ki,ij))
                    fcoa2=fcoa*sin(2.*az(ik,iray(j2,ir)))*
     *sin(2.*az(ki,ij))

c when the contribution of the first and last point of
c each great circle is not zero we take half of it.
c This should make the summation equivalent to a
c trapezoidal integration. TO BE CHECKED ON REGIOSTACK3

                    if(ki.eq.1.or.ki.eq.kki.or.
     &ik.eq.1.or.ik.eq.kdi(iray(j2,ir))) then
                       cov(iray(j2,ir))=cov(iray(j2,ir))+
     &0.5*(fco+fcoa1+fcoa2)*dsi(iray(j2,ir))*dsi(ij)
                       else
                       cov(iray(j2,ir))=cov(iray(j2,ir))+
     &(fco+fcoa1+fcoa2)*dsi(iray(j2,ir))*dsi(ij)
                       endif

                  endif
112    continue

                  endif

110    continue

c finding the neighbour of j2

                indice(1)=j2-nyp
                indice(2)=j2-ny
                indice(3)=j2-nym
                indice(4)=j2-1
                indice(5)=j2+1
                indice(6)=j2+nym
                indice(7)=j2+ny
                indice(8)=j2+nyp
c eastern boundary
             if (mod(j2,ny).eq.0.and.nypxy.eq.360.) then
                indice(3)=j2-2*ny+1
                indice(5)=j2-nym
                indice(8)=j2+1
             elseif(mod(j2,ny).eq.0.and.nypxy.ne.360.) then
                write(*,*)' STOP: ny*pxy =',nypxy
                stop'Agrandir la zone a l est'
             endif
c western boundary
             if (mod(j2,ny).eq.1.and.nypxy.eq.360.) then
                indice(1)=j2-1
                indice(4)=j2+nym
                indice(6)=j2+2*ny-1
             elseif(mod(j2,ny).eq.1.and.nypxy.ne.360.) then
                write(*,*)' STOP: ny*pxy =',nypxy
                stop'Agrandir la zone a l ouest'
             endif
c northern boundary and north pole
             if(j2.le.ny.and.tx2.le.pxy) then
                if(nypxy.ne.360.) then
                   write(*,*)' STOP: ny*pxy =',nypxy
                   stop 'Path at north pole requires to
     *     define a region covering the 360 deg. of longitude variation'
                endif
                if(mod(ny,2).ne.0) stop'uneven number of column?'
                if(j2.gt.int(fnys2)) then
                  indice(2)=j2-int(fnys2)
                  if(indice(2)-1.gt.0) then
                     indice(1)=indice(2)-1
                  else
                   indice(1)=ny
                  endif
                  indice(3)=indice(2)+1
                endif
                if(j2.le.int(fnys2)) then
                  indice(2)=j2+int(fnys2)
                  if(indice(2)+1.gt.ny) then
                     indice(3)=1
                  else
                     indice(3)=indice(2)+1
                  endif
                  indice(1)=indice(2)-1
                endif
             else if (j2.le.ny.and.tx2.gt.pxy) then
                write(*,*)'J2 TX2 PXY',j2,tx2,pxy
                stop'Agrandir la zone au Nord'
             endif
c southern boundary and south pole
             if(j2.gt.nynxm.and.tx2.ge.pimpxy) then
                if(nypxy.ne.360.) then
                    write(*,*)' STOP: ny*pxy =',nypxy
                    stop'Path at south pole requires to
     *     define a region covering the 360 deg. of longitude variation'
                endif
                if(mod(ny,2).ne.0) stop'uneven number of column?'
                if(j2.gt.int(fnys2)) then
                  indice(7)=j2-int(fnys2)
                  if(indice(7)-1.gt.nynxm) then
                     indice(6)=indice(7)-1
                  else
                     indice(6)=nxny
                  endif
                  indice(8)=indice(7)+1
                endif
                if(j2.le.int(fnys2)) then
                  indice(7)=j2+int(fnys2)
                  if(indice(7)+1.gt.nxny) then
c                   indice(8)=1
                    indice(8)=nynxm+1
                  else
                    indice(8)=indice(7)+1
                  endif
                  indice(6)=indice(7)-1
                endif
             else if(j2.gt.nynxm.and.tx2.lt.pimpxy) then
                stop'Agrandir la zone au Sud'
             endif

c------------------------------------------- 
 
                do 130 inb=1,8
                 jcurrent=indice(inb)
                  if(jcurrent.gt.KMAX.or.jcurrent.le.0) then
                     write(*,*) 'jcurrent INB tx2',jcurrent,inb,tx2
                     stop'Dimension array exceeded GCmGt'
                  endif 
                  if(mask(jcurrent)) then
                     call push(jcurrent)
                     kk=kk+1
                     if(kk.gt.KMAX) stop 
     *'Dimension array list exceeded GCmGt'
                     list(kk)=jcurrent
                     mask(jcurrent)=.FALSE.
                   endif
130                continue
             endif
 
             call stackempty(ibid)
             if(ibid.eq.1) go to 200
 
             call pop(jj)
             j2=jj
 
             goto 100
 
200          call stackflush()
        do k=1,kk
          mask(list(k))=.TRUE.
        end do
       return
       end
c-----------------------------------------------------------------
      SUBROUTINE linbcg(n,b,x,itol,tol,itmax,iter,err)
      INTEGER iter,itmax,itol,n,NTRAJ
      DOUBLE PRECISION err,tol,b(*),x(*),EPS
      PARAMETER (NTRAJ=10000,EPS=1.d-14)
CU    USES atimes,asolve,snrm
      INTEGER j
      DOUBLE PRECISION ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zm1nrm,
     *znrm,p(NTRAJ),pp(NTRAJ),r(NTRAJ),rr(NTRAJ),z(NTRAJ),zz(NTRAJ),snrm
      iter=0
      call atimes(n,x,r,0)
      do 11 j=1,n
        r(j)=b(j)-r(j)
        rr(j)=r(j)
11    continue
C     call atimes(n,r,rr,0)
      znrm=1.d0
      if(itol.eq.1) then
        bnrm=snrm(n,b,itol)
      else if (itol.eq.2) then
        call asolve(n,b,z,0)
        bnrm=snrm(n,z,itol)
      else if (itol.eq.3.or.itol.eq.4) then
        call asolve(n,b,z,0)
        bnrm=snrm(n,z,itol)
        call asolve(n,r,z,0)
        znrm=snrm(n,z,itol)
      else
        pause 'illegal itol in linbcg'
      endif
      call asolve(n,r,z,0)
100   if (iter.le.itmax) then
        iter=iter+1
        zm1nrm=znrm
        call asolve(n,rr,zz,1)
        bknum=0.d0
        do 12 j=1,n
          bknum=bknum+z(j)*rr(j)
12      continue
        if(iter.eq.1) then
          do 13 j=1,n
            p(j)=z(j)
            pp(j)=zz(j)
13        continue
        else
          bk=bknum/bkden
          do 14 j=1,n
            p(j)=bk*p(j)+z(j)
            pp(j)=bk*pp(j)+zz(j)
14        continue
        endif
        bkden=bknum
        call atimes(n,p,z,0)
        akden=0.d0
        do 15 j=1,n
          akden=akden+z(j)*pp(j)
15      continue
        ak=bknum/akden
        call atimes(n,pp,zz,1)
        do 16 j=1,n
          x(j)=x(j)+ak*p(j)
          r(j)=r(j)-ak*z(j)
          rr(j)=rr(j)-ak*zz(j)
16      continue
        call asolve(n,r,z,0)
        if(itol.eq.1.or.itol.eq.2)then
          znrm=1.d0
          err=snrm(n,r,itol)/bnrm
        else if(itol.eq.3.or.itol.eq.4)then
          znrm=snrm(n,z,itol)
          if(abs(zm1nrm-znrm).gt.EPS*znrm) then
            dxnrm=abs(ak)*snrm(n,p,itol)
            err=znrm/abs(zm1nrm-znrm)*dxnrm
          else
            err=znrm/bnrm
            goto 100
          endif
          xnrm=snrm(n,x,itol)
          if(err.le.0.5d0*xnrm) then
            err=err/xnrm
          else
            err=znrm/bnrm
            goto 100
          endif
        endif
        write (*,*) ' iter=',iter,' err=',err
      if(err.gt.tol) goto 100
      endif
c modif eric 20/10/20002
      IF(err.gt.tol) STOP 'err gt than tol, the desired 
     *convergence tolerance'
      return
      END
c-----------------------------------------------------------------
C  (C) Copr. 1986-92 Numerical Recipes Software Y"1p3-1.
      SUBROUTINE atimes(n,x,r,itrnsp)
      INTEGER n,itrnsp,ija,NMAXS
      DOUBLE PRECISION x(n),r(n),sa
      PARAMETER (NTRAJ=10000,NMAXS=NTRAJ*3500)
      COMMON /mat/ sa(NMAXS),ija(NMAXS)
CU    USES dsprsax,dsprstx
      if (itrnsp.eq.0) then
        call dsprsax(sa,ija,x,r,n)
      else
        call dsprstx(sa,ija,x,r,n)
      endif
      return
      END
c-----------------------------------------------------------------
C  (C) Copr. 1986-92 Numerical Recipes Software Y"1p3-1.
      SUBROUTINE asolve(n,b,x,itrnsp)
      INTEGER n,itrnsp,ija,NMAXS,i
      DOUBLE PRECISION x(n),b(n),sa
      PARAMETER (NTRAJ=10000,NMAXS=NTRAJ*3500)
      COMMON /mat/ sa(NMAXS),ija(NMAXS)
      do 11 i=1,n
        x(i)=b(i)/sa(i)
11    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software Y"1p3-1.
c-----------------------------------------------------------------
      SUBROUTINE dsprstx(sa,ija,x,b,n)
      INTEGER n,ija(*)
      DOUBLE PRECISION b(n),sa(*),x(n)
      INTEGER i,j,k
      if (ija(1).ne.n+2) pause 'mismatched vector and matrix in sprstx'
      do 11 i=1,n
        b(i)=sa(i)*x(i)
11    continue
      do 13 i=1,n
        do 12 k=ija(i),ija(i+1)-1
          j=ija(k)
          b(j)=b(j)+sa(k)*x(i)
12      continue
13    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software Y"1p3-1.
c-----------------------------------------------------------------
      SUBROUTINE dsprsax(sa,ija,x,b,n)
      INTEGER n,ija(*)
      DOUBLE PRECISION b(n),sa(*),x(n)
      INTEGER i,k
      if (ija(1).ne.n+2) pause 'mismatched vector and matrix in sprsax'
      do 12 i=1,n
        b(i)=sa(i)*x(i)
        do 11 k=ija(i),ija(i+1)-1
          b(i)=b(i)+sa(k)*x(ija(k))
11      continue
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software Y"1p3-1.
c-----------------------------------------------------------------
      FUNCTION snrm(n,sx,itol)
      INTEGER n,itol,i,isamax
      DOUBLE PRECISION sx(n),snrm
      if (itol.le.3)then
        snrm=0.
        do 11 i=1,n
          snrm=snrm+sx(i)**2
11      continue
        snrm=sqrt(snrm)
      else
        isamax=1
        do 12 i=1,n
          if(abs(sx(i)).gt.abs(sx(isamax))) isamax=i
12      continue
        snrm=abs(sx(isamax))
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software Y"1p3-1.
c-----------------------------------------------------------------
