c  Pgm for "continuous regionalization" using massive datsets
c  using massive dataset after Debayle and Sambridge, JGR (2003).
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Compared to the initial approach described 
c  in Montagner(1986, Annales Geophysicae),
c  the new scheme incorporates :
c  -geometrical algorithms which dramatically increases
c   computational efficiency (routines CmGtstack and GCmGtstack).
c   These routines have been //lized under OMP so that this
c   code can be used on a single PC processor or a // machine.
c  -a conjugate gradient method for sparse system (see numerical 
c   recipies pp 77) to avoid a direct inversion of (GCmGt +Cd)-1
c   (we solve instead (GCmGt +Cd)-1*(d-Gmo)). 
c  -an optimized storage of the (GCmGt +Cd) -1 arrays under compact form
c  -dynamic allocation of all arrays.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Previous modifs eric.
c 
c  Les derivees partielles des parametres anisotropes sont corrigees.
c  ON SORT Viso A1 A2 (Smith et Dahlen) de ce pgm alors que 
c  les parametres
c  inverses sont les (-A*vitesse a priori/vitesse inversee**2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      INTEGER, PARAMETER :: stack_size = 1000
      INTEGER :: store(stack_size),ipos
 
      integer OMP_GET_NUM_THREADS
      integer OMP_GET_THREAD_NUM
      REAL*8,DIMENSION(:),allocatable::sa
      INTEGER,DIMENSION(:),allocatable::ija

      REAL,DIMENSION(:,:),allocatable::tetai
      REAL,DIMENSION(:,:),allocatable::phii
      REAL,DIMENSION(:,:),allocatable::az
      REAL,DIMENSION(:,:),allocatable::dan1
      REAL,DIMENSION(:,:),allocatable::dan2
      REAL,DIMENSION(:,:),allocatable::pac
      REAL,DIMENSION(:,:),allocatable::ulf
      REAL,DIMENSION(:,:),allocatable::w
      REAL,DIMENSION(:),allocatable::deltei
      REAL,DIMENSION(:),allocatable::dul
      REAL,DIMENSION(:),allocatable::an1
      REAL,DIMENSION(:),allocatable::an2
      REAL,DIMENSION(:),allocatable::dul_tmp
      REAL,DIMENSION(:),allocatable::an1_tmp
      REAL,DIMENSION(:),allocatable::an2_tmp
      REAL,DIMENSION(:),allocatable::teta
      REAL,DIMENSION(:),allocatable::phi
      REAL*8,DIMENSION(:),allocatable::v
      REAL*8,DIMENSION(:),allocatable::vv
      REAL*8,DIMENSION(:),allocatable::uld
      REAL*8,DIMENSION(:),allocatable::res
      REAL,DIMENSION(:),allocatable::delt
      REAL,DIMENSION(:),allocatable::dsi
      REAL,DIMENSION(:),allocatable::ed
      REAL,DIMENSION(:),allocatable::resn
      REAL,DIMENSION(:),allocatable::cov
      REAL,DIMENSION(:),allocatable::t

      INTEGER,DIMENSION(:,:),allocatable::iray
      INTEGER,DIMENSION(:,:),allocatable::iss
      INTEGER,DIMENSION(:,:),allocatable::ise
      INTEGER,DIMENSION(:),allocatable::nray
      INTEGER,DIMENSION(:),allocatable::kdi
      INTEGER,DIMENSION(:),allocatable::jt

      CHARACTER*7,DIMENSION(:),allocatable::staz

      LOGICAL,DIMENSION(:),allocatable::mask

      real*8 err,tol,re

      character*32 fileout,fianiout,ficdata
      common/c6/xo,yo,nxo,nyo,pxyo
      common/c61/xinf,yinf,nx,ny,pxy
      common/c62/nxny,nxm,nym,nyp,nypxy,nynxm,fpxys2,fnys2,pimpxy
      common/c63/dcorlim,dcor2,dcoran2,sigma2,sigan2,fcoagcgt,fcoacgt
      common/c9/rad

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
c  xo,yo :initial colatitude and longitude
c  nxo,nyo :nb of points in latitude,longitude
c  pxyo :step in lat,long
      xo=0.
      yo=-180.
      print *,' enter pxyo such as mod(360,pxyo)=0 and mod(180,pxyo)=0'
      read(*,*) pxyo
      if(amod(360.,pxyo).ne.0)then
         write(*,*)' amod(360.,pxyo) must be equal to 0 '
         stop
      endif
      if(amod(180.,pxyo).ne.0)then
         write(*,*)' amod(180.,pxyo) must be equal to 0 '
         stop
      endif
        nxo=int(180/pxyo)
        nyo=int(360/pxyo)
      write(*,*)' xo,yo,nxo,nyo,pxyo ',xo,yo,nxo,nyo,pxyo
      xinf=xo
      yinf=yo
      nx=nxo
      ny=nyo
      pxy=pxyo
c  these values can be different from xinf yinf nx ny pxy because
c  it is possible to compute the final model in an area smaller
c  than the region covered by the paths (zooming)
c  Eric May 2003 :not necessary here because only the ''influence
c  zone'' around the path is considered in the computation so
c  that the whole Earth can be adopted both for xo,yo,nxo,nyo,pxyo
c  and xinf,yinf,nx,ny,pxy
c  xinf,yinf: colatitude and longitude of the origin point
c            located at the north-western corner of the grid
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
      print *,'enter number of depths/periods present'
      print *,'in input file'
      read(*,*)nper
 
c  Dynamic allocation vector NTRA--------------------
        ALLOCATE(t(NPER),stat=ierr)
        if(ierr/=0) then
          print *,"Erreur a L\'allocation de t"
          stop 4
        endif
        ALLOCATE(jt(NPER),stat=ierr)
        if(ierr/=0) then
          print *,"Erreur a L\'allocation de jt"
          stop 4
        endif
        ALLOCATE(v(NTRA),stat=ierr)
        if(ierr/=0) then
          print *,"Erreur a L\'allocation de v"
          stop 4
        endif
        ALLOCATE(vv(NTRA),stat=ierr)
        if(ierr/=0) then
          print *,"Erreur a L\'allocation de vv"
          stop 4
        endif
        ALLOCATE(uld(NTRA),stat=ierr)
        if(ierr/=0) then
          print *,"Erreur a L\'allocation de uld"
          stop 4
        endif
        ALLOCATE(res(NTRA),stat=ierr)
        if(ierr/=0) then
          print *,"Erreur a L\'allocation de res"
          stop 4
        endif
        ALLOCATE(delt(NTRA),stat=ierr)
        if(ierr/=0) then
          print *,"Erreur a L\'allocation de delt"
          stop 4
        endif
        ALLOCATE(dsi(NTRA),stat=ierr)
        if(ierr/=0) then
          print *,"Erreur a L\'allocation de dsi"
          stop 4
        endif
        ALLOCATE(ed(NTRA),stat=ierr)
        if(ierr/=0) then
          print *,"Erreur a L\'allocation de ed"
          stop 4
        endif
        ALLOCATE(resn(NTRA),stat=ierr)
        if(ierr/=0) then
          print *,"Erreur a L\'allocation de resn"
          stop 4
        endif
        ALLOCATE(kdi(NTRA),stat=ierr)
        if(ierr/=0) then
          print *,"Erreur a L\'allocation de kdi"
          stop 4
        endif
        ALLOCATE(staz(NTRA),stat=ierr)
        if(ierr/=0) then
          print *,"Erreur a L\'allocation de staz"
          stop 4
        endif
c  End dynamic allocation ntra--------------------------
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
      print *,'enter nbper,indices periodes a traiter'
      read(*,*)nbper,(jt(i),i=1,nbper)
      print *,nper,nbper,(jt(i),i=1,nbper)
      read(*,*)(t(i),i=1,nper)
      print *,(t(i),i=1,nper)

      print *,'enter data file name'
      read(5,'(a32)') ficdata
      print 999,ficdata
      open(10,file=ficdata,form='formatted')
 999  format(a32)

c  Dynamic allocation other arrays-----------------------
c  Note that  NF and NTR are linked to the geometry
c  and may need to be adapted by the user.
         
        if(NTRA.lt.6000) then
           NF=NTRA
        else
           NF=6000+int((NTRA-6000)/2)
        endif
c          NF=6000

        if(NTRA.lt.1000) then
           NTR=NTRA
        else if(NTRA.gt.1000.and.NTRA.le.5000) then
           NTR=int(NTRA/2)
        else if(NTRA.gt.5000.and.NTRA.le.10000) then
           NTR=int(NTRA/5)
        else if(NTRA.gt.10000.and.NTRA.le.20000) then
           NTR=int(NTRA/10)
        else if(NTRA.gt.20000) then
           NTR=int(NTRA/15)
        endif
c          NTR=1450

        KMAX=nx*ny
        NTDEPTH=2*nper+6
        ALLOCATE(w(NTDEPTH,NTRA),stat=ierr)
        if(ierr/=0) then
          print *,"Erreur a L\'allocation de w"
          stop 4
        endif

      call lecdatV3(randerr,nper,staz,delt,kdi,dsi,ntra,d1,
     *t,NSTEP,w,NTDEPTH)

        write(*,*)'KMAX, NTR, NSTEP ',KMAX,NTR,NSTEP
        ALLOCATE(dul(KMAX),stat=ierr)
        if(ierr/=0) then
          print *,"Erreur a L\'allocation de dul"
          stop 4
        endif
        ALLOCATE(an1(KMAX),stat=ierr)
        if(ierr/=0) then
          print *,"Erreur a L\'allocation de an1"
          stop 4
        endif
        ALLOCATE(an2(KMAX),stat=ierr)
        if(ierr/=0) then
          print *,"Erreur a L\'allocation de an2"
          stop 4
        endif
        ALLOCATE(teta(KMAX),stat=ierr)
        if(ierr/=0) then
          print *,"Erreur a L\'allocation de teta"
          stop 4
        endif
        ALLOCATE(phi(KMAX),stat=ierr)
        if(ierr/=0) then
          print *,"Erreur a L\'allocation de phi"
          stop 4
        endif
        ALLOCATE(nray(KMAX),stat=ierr)
        if(ierr/=0) then
          print *,"Erreur a L\'allocation de nray"
          stop 4
        endif
        ALLOCATE(dan1(NX,NY),stat=ierr)
        if(ierr/=0) then
          print *,"Erreur a L\'allocation de dan1"
          stop 4
        endif
        ALLOCATE(dan2(NX,NY),stat=ierr)
        if(ierr/=0) then
          print *,"Erreur a L\'allocation de dan2"
          stop 4
        endif
        ALLOCATE(pac(NX,NY),stat=ierr)
        if(ierr/=0) then
          print *,"Erreur a L\'allocation de pac"
          stop 4
        endif
        ALLOCATE(ulf(NX,NY),stat=ierr)
        if(ierr/=0) then
          print *,"Erreur a L\'allocation de ulf"
          stop 4
        endif
        ALLOCATE(iray(KMAX,NTR),stat=ierr)
        if(ierr/=0) then
          print *,"Erreur a L\'allocation de iray"
          stop 4
        endif
        ALLOCATE(iss(KMAX,NTR),stat=ierr)
        if(ierr/=0) then
          print *,"Erreur a L\'allocation de iss"
          stop 4
        endif
        ALLOCATE(ise(KMAX,NTR),stat=ierr)
        if(ierr/=0) then
          print *,"Erreur a L\'allocation de ise"
          stop 4
        endif
        ALLOCATE(tetai(NSTEP,NTRA),stat=ierr)
        if(ierr/=0) then
          print *,"Erreur a L\'allocation de tetai"
          stop 4
        endif
        ALLOCATE(phii(NSTEP,NTRA),stat=ierr)
        if(ierr/=0) then
          print *,"Erreur a L\'allocation de phii"
          stop 4
        endif
        ALLOCATE(az(NSTEP,NTRA),stat=ierr)
        if(ierr/=0) then
          print *,"Erreur a L\'allocation de az"
          stop 4
        endif        
        ALLOCATE(SA(NTRA*NF),stat=ierr)
        if(ierr/=0) then
          print *,"Erreur a L\'allocation de SA"
          stop 4
        endif
        ALLOCATE(IJA(NTRA*NF),stat=ierr)
        if(ierr/=0) then
          print *,"Erreur a L\'allocation de IJA"
          stop 4
        endif

        ALLOCATE(deltei(NSTEP),stat=ierr)
        if(ierr/=0) then
          print *,"Erreur a L\'allocation de deltei"
          stop 4
        endif
c  End dynamic allocation other arrays-----------------

        do 4 i=1,ntra
        kd1=kdi(i)
        dds=dsi(i)/rad
        ddel=delt(i)/rad
c FOR DEBUGGING
c       print 2005,kd1,dds,delt(i),ddel
c2005   format('nbr of pts,ds,delta(rad),delta(o)',i3,3f8.3)
        tm1=w(1,i)/rad
        pm1=w(2,i)/rad
        tm2=w(3,i)/rad
        pm2=w(4,i)/rad
c FOR DEBUGGING
c       print 2006,tm1,pm1,tm2,pm2
c2006   format('t1 p1  t2 p2',4f8.2)
        do 4 ki=1,kd1
          gi=float(ki-1)
          deltei(ki)=gi*dsi(i)
          az(ki,i)=0.
          call coordo(w(1,i),w(2,i),w(3,i),w(4,i),tetai(ki,i),
     *phii(ki,i),deltei(ki),az(ki,i))
  4       continue

      print *,'enter final velocity file'
      read(5,999)fileout
      print 999,fileout
      open(7,file=fileout)
      print *,'enter final anisotropy file'
      read(5,999)fianiout
      print 999,fianiout
      open(9,file=fianiout)
c  dcor, dcoran correlation lengths.
      print *,'enter dcor,dcoran (in km) for velocity and anisotropy'
      read(*,*)dcor,dcoran
      print *,dcor,dcoran
      dcor=dcor*rad/111.1949
      dcoran=dcoran*rad/111.1949

c      loop on periods
       k1=1
       do 5 ktt=1,nbper

       ntt=jt(ktt)
       k2=ntt
c  in vinit,one reads starting velocity and anisotropy
c  at a given period, data are sorted because the number
c  of data can be period dependent
      print *,' '
      do 55 jk=k1,k2
        call vinit(pac,ulf,dan1,dan2,nx,sigma,sigan,ut0,tt)
 55   k1=k2+1
      write(7,*)tt
      write(9,*)tt
      do 6 j=1,ntra
      nk=ntt+6
c     nkk=ntt+NDEPTH+6
      nkk=ntt+nper+6
      uld(j)=w(nk,j)
      ed(j)=w(nkk,j)
c     write(*,*)'J ULD W ED W ',j,uld(j),w(nk,j),ed(j),w(nkk,j)
      w(5,j)=(-uld(j))*cos(2.*w(5,j))
  6   continue
      print *,' '
      print 2007,t(ntt),ntra,dcor
 2007 format('periode,nda,dcor',f6.2,i6,f6.2)
c-------------------------------------------------------
c modif eric
c on calcul les moyennes et ecarts-types en lenteur
c car on inverse en lenteur.
c FOR DEBUGGING
c     print *,'selected data '
c     print 2003,(1./uld(i),i=1,ntra)
      call me(uld,xm,xe,ntra)
c-------------------------------------------------------
      write(*,*) 'moyenne et ecart-type de uld (1/donnees)',xm,xe
      varmoy=xe*xe
      write(*,*)'variance-donnees p.r. a la moyenne des mesures:',varmoy
      write(*,*) 'modele initial= ',ut0,
     &   '           moyenne des mesures= ',1/xm
      x0=xe
 1100 format(' initial values : xm= ',f7.3,' xe= ',f9.5)

c---------------------------------------------------------------
cComputation of (GCmGt +Cdo) (E. Debayle march 2002)
c NB Eric, may 2003
c The computation of GCmGt has not been yet extracted from the loop
c on the depth (or period) because
c Cm=sig(M1)*sig(M2)**2*exp(1/2(dis(M1,M2)/Lcorr)**2)
c and the standard deviation sig(M1) is assumed to be constant 
c with depth in velocity which means for slowness a dependance
c on the a priori velocity at each depth (sig(M1)/Vs(z)**2).
c Assuming Vs independent of depth is equivalent to assume large
c sig(M1) with depth.
c As the anisotropic S is computed following
c S= (Cd + (GCmGt)_vs + (GCmGt)_A1 +  (GCmGt)_A2))
c we could remove the isotropic and anisotropic standard deviations
c from the computation of the integrals present in the continuous
c expression of S.
c However this would require to store the 3 matrix (GCmGt)_vs 
c (GCmGt)_A1 and (GCmGt)_A2 which each have the size of arrays sa and ija.
c This is too expensive in memory for large systems.

      do 50 j=1,KMAX
            phi(j)=0.
            teta(j)=0.
            nray(j)=0
c        do 50 jj=1,NTR
c           iss(j,jj)=0.
c           ise(j,jj)=0.
c           iray(j,jj)=0.
 50   continue
 
      do 54 ii=1,nx
         do 54 jj=1,ny
          j2=ny*(ii-1)+jj
          if(j2.gt.KMAX) then
            write(*,*)'NX NY ii jj j2',nx,ny,ii,jj,j2
          endif
          teta(j2)=xinf+(ii-1)*pxy+pxy/2.
          teta(j2)=teta(j2)*rad
          phi(j2)=yinf+(jj-1)*pxy+pxy/2.
          phi(j2)=phi(j2)*rad
 54   continue

      jold=0
      ntrmax=0
      do 70 i=1,ntra
        tetai(1,i)=w(1,i)
        phii(1,i)=w(2,i)
        kki=kdi(i)
        tetai(kki,i)=w(3,i)
        phii(kki,i)=w(4,i)

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

           if(i.eq.ntra.and.ki.eq.kki.and.j2.eq.jold)ise(j2,nray(j2))=ki

           if(j2.ne.jold.or.ki.eq.1) then
             nray(j2)=nray(j2)+1
             if(nray(j2).gt.ntrmax)ntrmax=nray(j2)
             if(nray(j2).gt.NTR)then
                 write(*,*)'J2 ',j2,' nray(j2) ',nray(j2)
                 write(*,*)'NTR ',NTR
                 stop' NTR must be increased'
             endif
             iray(j2,nray(j2))=i
             iss(j2,nray(j2))=ki
             if(i.eq.ntra.and.ki.eq.kki)then
                 write(*,*) 'Last point of last path'
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
       write(*,*)'Number max of paths in a cell ',ntrmax
       write(*,*)'Number used for arrays allocation ',NTR

      fpxys2=pxy/2.
      nxny=nx*ny
      nxm=nx-1
      nym=ny-1
      nyp=ny+1
      nypxy=ny*pxy
      fnys2=float(ny)/2.
      nynxm=ny*nxm
      pimpxy=180.-pxy
      dcorlim=2.648*dcor
c     dcorlim=3.03*dcor
      dcor2=dcor**2
      dcoran2=dcoran**2
      sigma2=sigma**2
      sigan2=sigan**2
      fcoagcgt=(1/ut0)**2*sigan2
      fcoacgt=-(1/ut0)*sigan2

       k=ntra+1
       ija(1)=ntra+2
!$omp parallel private(j,mask,cov)
       ALLOCATE(cov(NTRA),stat=ierr)
       if(ierr/=0) then
         print *,"Erreur a L\'allocation de cov"
         stop 4
       endif
       ALLOCATE(mask(KMAX),stat=ierr)
       if(ierr/=0) then
         print *,"Erreur a L\'allocation de mask"
         stop 4
       endif
       do j=1,KMAX
         mask(j)=.TRUE.
       enddo
       do j=1,ntra
         cov(j)=0.
       enddo
!$    print *,'nbre de T:',OMP_GET_NUM_THREADS()
!$omp do ordered schedule(dynamic,1) private(i,kki,ki,ij,it,jj,
!$omp+ jcurrent,store,ipos)
cam!$omp+ jcurrent,store,ipos) firstprivate(mask,cov)
       do 100 i=1,ntra
         kki=kdi(i)
         jj=0
         do 110 ki=1,kki
         call GCmGtstack(KMAX,iray,iss,ise,tetai,phii,az,nstep,cov,
     *mask,i,ki,teta,phi,dsi,kdi,nray,jj,jcurrent,store,ipos)
  110 continue

c we store the matrix (GCmGt+Cdo) in row indexed sparse
c storage mode (see num. recipies p 71)
c we also reconstruct the full matrix (which is symetric 
c from its upper triangle part) (when if(ij.lt.i)  test is verified)

c!$       write(*,*) 'AV ORDERED ',i,' rang ',OMP_GET_THREAD_NUM()
!$omp ordered

c      do 105 ij=1,ntra
       do ij=1,ntra
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
          else if(abs(cov(ij)).gt.1e-9)then
             k=k+1
             sa(k)=cov(ij)/(delt(i)*delt(ij))
             ija(k)=ij
          endif
          cov(ij)=0.
      enddo
c105   continue
       ija(i+1)=k+1
c!$     write (*,*)'IKSAMAX ',k+1,' rang ',OMP_GET_THREAD_NUM()
!$omp end ordered
 100   continue
!$omp end do
      DEALLOCATE(mask,cov)
!$omp end parallel

       write (*,*)'Size arrays SA IJA ',k+1
     *,' allocated size is ',NTRA*NF
      if((k+1).gt.NTRA*NF) then
               write (*,*)'Size arrays SA IJA too small '
               write (*,*)'value of NF must be increased '
               write (*,*)'in main code '
               stop
       endif

c--------------------------------------------------------------
c   forward problem resolution.
c   Vector v contains (d-g(p))in slowness

      do 130 i=1,ntra
      v(i)=0.
      vv(i)=0.
      call pbdirc(pac,dan1,dan2,nx,tetai,phii,az,NSTEP,uld(i),v(i),
     *kdi(i),delt(i),dsi(i),i)
c     write(*,*)'i uld v kdi delt dsi'
c     write(*,*)i,uld(i),v(i),kdi(i),delt(i),dsi(i)
 130  continue
c FOR DEBUGGING
c     print*,'vecteur v (d-g(p) en lenteur):'
c     print 2003,(v(i),i=1,ntra)
cJJL  Calcul de la variance des donnees par rapport aux synthetiques
cJJL  du modele initial (et non par rapport a la moyenne des donnees)
      varini=0.
      do 139 i=1,ntra
      varini=varini+v(i)*v(i)
139   continue
      varini=varini/ntra
      write(*,*)'variance-donnees p.r. au modele initial:       ',varini
c  modif eric : calcul du qui initial en lenteur 
      call khi(v,ed,qui,ntra)
      print*, 'qui initial :',qui

c  solving the linear system (GCmGt +Cdo)-1*(d-Gmo) using
c  conjugate gradient

      itol=3
      itmax=1000
      tol= 1e-12
      iter=0
      err=0

      write(*,*)' INPUT LINBCG '
      write(*,*) 'NDA ',ntra,' ITOL ',itol,' TOL ',tol,' ITMAX ',itmax
      write(*,*) 'ITER ',iter,' ERR ',err
c     print *,'vv(i)'
c     write(*,*) (vv(i),i=1,ntra)
      call linbcg(sa,ija,ntra,v,vv,itol,tol,itmax,iter,err)

 2012 format(10(f9.6,1x))

      call me(vv,xm,xe,ntra)
      print 1100,xm,xe
c###################################################################### 
c  Here we compute m=mO+CmGt(GCmGt+Cdo)-1(d-Gmo)
c  solution of the inverse pb using the CmGtstack2 routine
c NB :
c  repere 1 : xo,yo,nxo,nyo,pxyo
c  repere 2 : xinf,yinf,nx,ny,pxy

      do 9950 j=1,KMAX
            dul(j)=0.
            an1(j)=0.
            an2(j)=0.
9950  continue
!$omp parallel private(j,dul_tmp,an1_tmp,an2_tmp,mask)
      ALLOCATE(dul_tmp(KMAX),stat=ierr)
      if(ierr/=0) then
          print *,"Erreur a L\'allocation de dul"
          stop 4
      endif
      ALLOCATE(an1_tmp(KMAX),stat=ierr)
      if(ierr/=0) then
          print *,"Erreur a L\'allocation de an1"
          stop 4
      endif
      ALLOCATE(an2_tmp(KMAX),stat=ierr)
      if(ierr/=0) then
          print *,"Erreur a L\'allocation de an2"
          stop 4
      endif
      ALLOCATE(mask(KMAX),stat=ierr)
      if(ierr/=0) then
        print *,"Erreur a L\'allocation de mask"
        stop 4
      endif
      do  j=1,KMAX
            mask(j)=.TRUE.
            dul_tmp(j)=0.
            an1_tmp(j)=0.
            an2_tmp(j)=0.
      enddo
!$omp do private(ij,ki,dk1,store,ipos)
cam!$omp+ firstprivate(mask)
      do 99110 ij=1,ntra
        kd1=kdi(ij)
        do 99120 ki=1,kd1
        call CmGtstack2(KMAX,tetai,phii,az,nstep,vv(ij),dsi(ij),mask,
     *delt(ij),ij,ki,kd1,store,ipos,dul_tmp,an1_tmp,an2_tmp,teta,phi)
99120   continue
99110 continue
!$omp end do
!$omp critical
      do  j=1,KMAX
       dul(j)=dul(j)+dul_tmp(j)
       an1(j)=an1(j)+an1_tmp(j)
       an2(j)=an2(j)+an2_tmp(j)
      enddo
!$omp end critical
      DEALLOCATE(dul_tmp,an1_tmp,an2_tmp,mask)
!$omp end parallel

      do 99140 ii=1,nx
        tx=xinf+(ii-1)*pxy+fpxys2
        it1=int((tx-(xo-pxyo))/pxyo)

        do 99140 jj=1,ny
          py=yinf+(jj-1)*pxy+fpxys2
          ip1=int((py-(yo-pxyo))/pxyo)
          if(it1.gt.NX.or.it1.le.0.or.ip1.gt.NY.or.ip1.le.0)then
              write(*,*)'IT1 ',it1,'IP1',ip1
              stop 'pb dimension'
          endif
          j2=ny*(ii-1)+jj

          dan1(it1,ip1)=dan1(it1,ip1)+an1(j2)
          dan2(it1,ip1)=dan2(it1,ip1)+an2(j2)
          ulf(it1,ip1)=dul(j2)+1./pac(it1,ip1)
          ulf(it1,ip1)=1./ulf(it1,ip1)

c     The calculation of the final a posteriori covariance
c     matrix is skipped as it would require to solve m times
c     (m beeing the number of model parameter) the
c     system (GCmGt +Cdo)-1*(each column of GCm)
c     using the conjugate gradient method

99140  continue
 
c###################################################################### 
      do 94 i=1,ntra
      re=0.
      call pbdirc(ulf,dan1,dan2,nx,tetai,phii,az,NSTEP,uld(i),re,
     *kdi(i),delt(i),dsi(i),i)
c   modif eric
c   On calcul le residus en lenteur 
      res(i)=re
      resn(i)=sqrt((res(i)/ed(i))**2)
  94  continue
      print*,'residus normalises par trajets (et par couche):'
      print 1007,(staz(i),resn(i),i=1,ntra)
      icomp=0
c--x0 et xe ont ete passe en lenteur, normalement redvar est en lenteur 

      call me(res,xm,xe,ntra)

      write(*,*) 'moyenne et ecart-type de res (residus fin)',xm,xe
      varres=xm*xm+xe*xe
      write(*,*)'variance-donnees p.r. a la moyenne des mesures:',varmoy
      write(*,*)'variance-donnees p.r. au modele final:         ',varres
      redvar =(varmoy-varres)/varmoy
      redvar2=(varini-varres)/varini

      print *,'red. variance(%val moy des donnees en lenteur)',redvar
      print *,'red. variance(% modele initial)     en lenteur',redvar2
c modif eric : ed est l'erreur en 1/V
c On doit donc calculer un qui en 1/V
c pour khim : sigma1 est passe en lenteur dans vinit. On le garde
c comme ca.
      call khi(res,ed,qui,ntra)
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
      print*,'khim (pert.quadr. normal. du modele en lenteur)= ',quim

c     ecriture du resultat
      do 127 i=1,nxo
 127  write(7,2014)(ulf(i,j),j=1,nyo)
c modif eric Dec 1995
c Les derivees partieles des parametres anisotropes sont (-cos2teta)/L*ut0
c avec  ut0=vitesse initiale. Ceci implique(voir these)que les 
c parametres inverses sont dan=(A1*ut0)/(ulf*ulf) ulf etant a ce stade
c la vitesse inversee
c Dans ce cas les parametres Ai (coefficients de smiths et Dahlen)sont: 

      do 1210 i=1,nxo
      do 1200 j=1,nyo
      dan1(i,j)=dan1(i,j)*(ulf(i,j)*ulf(i,j))/ut0
      dan2(i,j)=dan2(i,j)*(ulf(i,j)*ulf(i,j))/ut0
1200   continue
1210  continue
 
c apres cette conversion, le pgm sort les coefficients bruts de Smith et Dahlen
c tels que donnes dans la formule C= C0+A1cos2phi+A2sin2phi
c     ecriture des parametres anisotropes.
      do 125 i=1,nxo
 125  write(9,2014)(dan1(i,j),j=1,nyo)
      do 124 i=1,nxo
 124  write(9,2014)(dan2(i,j),j=1,nyo)
  5   continue
      close(10)
      close(7)
      close(9)

 1007 format(6(a7,1x,f5.2,1x))
 2003 format(10f7.4)
 2014 format(20(f10.6,1x))
       end
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine trapez(npoint,f,ds,sint)
      real f(*)
      sint=0.
      nmax=npoint-1
      do 1 ns=2,nmax
  1   sint=sint+2.*f(ns)
      sint=sint+f(1)+f(npoint)
      sint=sint*ds*0.5
      return
      end
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine lecdatV3(randerr,nper,staz,
     *delt,kdi,dsi,ndat,d1,t,NSTEP,w,NTDEPTH)
      
      real w(NTDEPTH,*)
      real delt(*),dsi(*)
      integer kdi(*)
      real t(nper),ug(nper),edat(nper)
      character*7 staz(*)

      common/c9/rad

c This subroutine reads for each epicenter path
c the coordinates of the epicenters and stations,
c the seismic velocities (path-average SV,
c group or phase velocities) and their errors.
c w(1,i),w(2,i)=lat lon of event 
c w(3,i),w(4,i)=lat lon of station 
      nt=nper
      NSTEP=1
      do 1 i=1,ndat
        read(10,'(a7)')staz(i)
c       write(*,'(a7)')staz(i)
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
        read(10,*)(ug(j),j=1,nt)
c ug : Vs velocity at each depth or,
c group or phase velocity data for path i.
c print 1999,(ug(j),j=1,nt)
c edat error on ug
        read(10,*)(edat(j),j=1,nt)
        do 2 k=1,nt
          kk=k+6
          w(kk,i)=1/ug(k)
          kkk=k+nper+6
c modif eric may 2003
c Uncertainties in the origin time
c in the source mecanism and possible mislocations
c in the epicenter may not be reflected in the path
c average a posteriori errors (this would be the case 
c for waveform modelling).
c Assuming the errors are normally distributed we add
c them to the error in the data file :
          w(kkk,i)=sqrt(edat(k)**2+randerr**2)
          w(kkk,i)=w(kkk,i)/(ug(k)**2)
  2     continue
        in=1
c2002   format(2f8.3,3f7.2)
c       print 2008,(w(k,i),k=7,nper+6)
c2008   format('1/ug',13f6.3)
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
        if(kdi(i).gt.NSTEP)NSTEP=kdi(i)
  1     continue
      return
      end
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine vinit(pac,ulf,dan1,dan2,nx,sigma,sigan,ut0,tt)

      real pac(nx,*),ulf(nx,*)
      real dan1(nx,*),dan2(nx,*)

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
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine pbdirc(pac,dan1,dan2,nx,tetai,phii,az,NSTEP,uld,v,ki1,
     *delt,dsi,i)

      REAL,DIMENSION(:),allocatable::ulo
      real pac(nx,*),dan1(nx,*),dan2(nx,*)

      real tetai(NSTEP,*),phii(NSTEP,*),az(NSTEP,*)

      real*8 v,uld

      common/c6/xo,yo,nxo,nyo,pxyo
      common/c9/rad
c     calculation by integration of:
c     1/u=1/delt dsi/u(1+dan1*cos2phi+dan2*sin2phi)
c
      ALLOCATE(ulo(NSTEP),stat=ierr)
        if(ierr/=0) then
          print *,"Erreur a L\'allocation de ulo"
          stop 4
      endif
      v=0.
c     print *,'pbdirc'
      do 103 kki=1,ki1
      cs=cos(2.*az(kki,i))
      sn=sin(2.*az(kki,i))
      tet=tetai(kki,i)/rad
      ph=phii(kki,i)/rad
c  MODIF ERIC 19 SEPT (and 6/08/03)--------------
      if (ph.ge.180.) ph=ph-360.
      if (ph.lt.-180.) ph=ph+360.
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
      if(icp.ge.nyo)icp=nyo
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
c       print 996,pc1,pc2,pc3,pc4,ict,icp,ict1,icp1
        sd=d1*d2*d3+d2*d3*d4+d3*d4*d1+d4*d1*d2
        um=pc1*d2*d3*d4+pc2*d1*d3*d4+pc3*d1*d2*d4+pc4*d1*d2*d3
c---MODIF
c       if(i.eq.5)write(*,*)'UM,SD',um,sd
        um=um/sd
      endif
c     print 997,um
 997  format('um',f7.3)
      ulo(kki)=1./um
c     write(*,*)'sd um ',sd,um
 103  continue
 996  format('tet,ph,tei,fi',4f7.2,'ict,icp,ict1,icp1',4(1x,i3))
 998  format('d1,d2,d3,d4,pac1,pac4',6f7.4)
c     print 1000,(ulo(kki),kki=1,ki1)
 1000 format('ulo',10f6.3)
      call trapez(ki1,ulo,dsi,ss)
      v=-ss/delt
c     write(*,*)'delt v',delt,v
      v=uld+v
      return
      end
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine me(x,xm,xe,n)
      real*8 x(n)
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
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      function cosdel(tetam1,phim1,tetam2,phim2)
      cosdel=0.
      c1=cos(tetam1)
      c2=cos(tetam2)
      s1=sin(tetam1)
      s2=sin(tetam2)
      cosdel=c1*c2+(s1*s2)*cos(phim1-phim2)
      return
      end
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      function sindel(tetam1,phim1,tetam2,phim2)
      co=cosdel(tetam1,phim1,tetam2,phim2)
      un=1.0000000000
      a=abs(un-co**2)
      sindel=sqrt(a)
      return
      end
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      function de(t1,p1,t2,p2)
      common/c9/rad
      t=(t1+t2)*rad/2.
      de=((t2-t1)**2)+((p2-p1)**2)*(sin(t)**2)
      de=sqrt(de)
      return
      end
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine khi(x,ed,qui,n)
      real ed(n)
      real*8 x(n)
      qui=0.
      do 1 i=1,n
      xed=x(i)/ed(i)
      qui=qui+xed**2
   1  continue
      qui=qui/float(n)
      qui=sqrt(qui)
      return
      end
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        SUBROUTINE CmGtstack2(KMAX,tetai,phii,az,nstep,vv,dsi,mask
     *,delt,ij,ki,kd1,store,ipos,dul_tmp,an1_tmp,an2_tmp,teta,phi)

        INTEGER, PARAMETER :: stack_size = 1000
        INTEGER :: store(stack_size),ipos
 
        real*8 vv
        real dul_tmp(*),an1_tmp(*),an2_tmp(*)
        real phi(*),teta(*)

        real tetai(nstep,*),phii(nstep,*),az(nstep,*)

        integer indice(8),list(nxny)
        logical   mask(*)

        common/c6/xo,yo,nxo,nyo,pxyo
        common/c61/xinf,yinf,nx,ny,pxy
        common/c62/nxny,nxm,nym,nyp,nypxy,nynxm,fpxys2,fnys2,pimpxy
        common/c63/dcorlim,dcor2,dcoran2,sigma2,sigan2,fcoagcgt,fcoacgt
        common/c9/rad

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

c            call stackinit()

             ipos=0 
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
                 dul_tmp(j2)=dul_tmp(j2)+(0.5*fco*dsi/delt)*vv
                 an1_tmp(j2)=an1_tmp(j2)+(0.5*fcoa1*dsi/delt)*vv
                 an2_tmp(j2)=an2_tmp(j2)+(0.5*fcoa2*dsi/delt)*vv
                else
                 dul_tmp(j2)=dul_tmp(j2)+(fco*dsi/delt)*vv
                 an1_tmp(j2)=an1_tmp(j2)+(fcoa1*dsi/delt)*vv
                 an2_tmp(j2)=an2_tmp(j2)+(fcoa2*dsi/delt)*vv
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
                if((j2-nynxm).gt.int(fnys2)) then
                  indice(7)=j2-int(fnys2)
                  if(indice(7)-1.gt.nynxm) then
                     indice(6)=indice(7)-1
                  else
                     indice(6)=nxny
                  endif
                  indice(8)=indice(7)+1
                endif
                if((j2-nynxm).le.int(fnys2)) then
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
                     call pushf(store,jcurrent,ipos)
                     kk=kk+1
                     if(kk.gt.KMAX) stop 'Dimension array list exceeded'
                     list(kk)=jcurrent
                     mask(jcurrent)=.FALSE.
                   endif
                end do
             endif
 
c            call stackempty(ibid)
c            if(ibid.eq.1) go to 200

             call popf(store,jj,ipos)
                if (jj.gt.0) then
                j2=jj
                elseif (jj.eq.0) then
                goto 200
                endif
 
             goto 100
 
c200          call stackflush()
200           continue
        do k=1,kk
          mask(list(k))=.TRUE.
        end do
       return
       end
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE GCmGtstack(KMAX,iray,iss,ise,tetai,phii,az,nstep,
c    *cov,mask,ij,ki,teta,phi,dsi,kdi,nray,store,ipos)
     *cov,mask,ij,ki,teta,phi,dsi,kdi,nray,jj,jcurrent,store,ipos)

        INTEGER, PARAMETER :: stack_size = 1000
        INTEGER :: store(stack_size),ipos

        real  dsi(*),cov(*)
        integer  kdi(*)

        real phi(*),teta(*)
        integer nray(*)
 
        integer iray(KMAX,*),iss(KMAX,*),ise(KMAX,*)

        real tetai(nstep,*),phii(nstep,*),az(nstep,*) 

        integer indice(8),list(nxny)
        logical   mask(*)
 
        common/c6/xo,yo,nxo,nyo,pxyo
        common/c61/xinf,yinf,nx,ny,pxy
        common/c62/nxny,nxm,nym,nyp,nypxy,nynxm,fpxys2,fnys2,pimpxy
        common/c63/dcorlim,dcor2,dcoran2,sigma2,sigan2,fcoagcgt,fcoacgt
        common/c9/rad

        save sdmimj

        kki=kdi(ij) 
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
 
c            call stackinit()

             ipos=0 
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
                if((j2-nynxm).gt.int(fnys2)) then
                  indice(7)=j2-int(fnys2)
                  if(indice(7)-1.gt.nynxm) then
                     indice(6)=indice(7)-1
                  else
                     indice(6)=nxny
                  endif
                  indice(8)=indice(7)+1
                endif
                if((j2-nynxm).le.int(fnys2)) then
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
                     call pushf(store,jcurrent,ipos)
                     kk=kk+1
                     if(kk.gt.KMAX) stop 
     *'Dimension array list exceeded GCmGt'
                     list(kk)=jcurrent
                     mask(jcurrent)=.FALSE.
                   endif
130                continue
             endif
 
c            call stackempty(ibid)
c            if(ibid.eq.1) go to 200
 
             call popf(store,jj,ipos)
                if (jj.gt.0) then
                j2=jj
                elseif (jj.eq.0) then
                goto 200
                endif
                                                                                                
             goto 100
 
c200          call stackflush()
200           continue
        do k=1,kk
          mask(list(k))=.TRUE.
        end do
       return
       end
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE linbcg(sa,ija,n,b,x,itol,tol,itmax,iter,err)
c     INTEGER iter,itmax,itol,n,NTRAJ,ija(*)
      INTEGER iter,itmax,itol,n,ija(*)
      real*8 err,tol,b(*),x(*),EPS,sa(*)
c     PARAMETER (NTRAJ=6000,EPS=1.d-14)
      PARAMETER (EPS=1.d-14)
CU    USES atimes,asolve,snrm
      INTEGER j
      real*8 ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zm1nrm,
     *znrm,p(n),pp(n),r(n),rr(n),z(n),zz(n),snrm
c    *znrm,p(NTRAJ),pp(NTRAJ),r(NTRAJ),rr(NTRAJ),z(NTRAJ),zz(NTRAJ),snrm
      iter=0
c     call atimes(n,x,r,0)
      call atimes(sa,ija,n,x,r,0)
      do 11 j=1,n
        r(j)=b(j)-r(j)
        rr(j)=r(j)
11    continue
C     call atimes(n,r,rr,0)
C     call atimes(sa,ija,n,r,rr,0)
      znrm=1.d0
      if(itol.eq.1) then
        bnrm=snrm(n,b,itol)
      else if (itol.eq.2) then
c       call asolve(n,b,z,0)
        call asolve(sa,ija,n,b,z,0)
        bnrm=snrm(n,z,itol)
      else if (itol.eq.3.or.itol.eq.4) then
c       call asolve(n,b,z,0)
        call asolve(sa,ija,n,b,z,0)
        bnrm=snrm(n,z,itol)
c       call asolve(n,r,z,0)
        call asolve(sa,ija,n,r,z,0)
        znrm=snrm(n,z,itol)
      else
        pause 'illegal itol in linbcg'
      endif
c     call asolve(n,r,z,0)
      call asolve(sa,ija,n,r,z,0)
100   if (iter.le.itmax) then
        iter=iter+1
        zm1nrm=znrm
c       call asolve(n,rr,zz,1)
        call asolve(sa,ija,n,rr,zz,1)
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
c       call atimes(n,p,z,0)
        call atimes(sa,ija,n,p,z,0)
        akden=0.d0
        do 15 j=1,n
          akden=akden+z(j)*pp(j)
15      continue
        ak=bknum/akden
c       call atimes(n,pp,zz,1)
        call atimes(sa,ija,n,pp,zz,1)
        do 16 j=1,n
          x(j)=x(j)+ak*p(j)
          r(j)=r(j)-ak*z(j)
          rr(j)=rr(j)-ak*zz(j)
16      continue
c       call asolve(n,r,z,0)
        call asolve(sa,ija,n,r,z,0)
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
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  (C) Copr. 1986-92 Numerical Recipes Software Y"1p3-1.
      SUBROUTINE atimes(sa,ija,n,x,r,itrnsp)
c     REAL*8,DIMENSION(:),allocatable::sa
c     INTEGER*4,DIMENSION(:),allocatable::ija
c     INTEGER n,itrnsp,ija(*),NTRAJ
      INTEGER n,itrnsp,ija(*)
      real*8 x(n),r(n),sa(*)
c     PARAMETER (NTRAJ=6000)
c     COMMON /mat/ sa,ija
CU    USES dsprsax,dsprstx
      if (itrnsp.eq.0) then
        call dsprsax(sa,ija,x,r,n)
      else
        call dsprstx(sa,ija,x,r,n)
      endif
      return
      END
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  (C) Copr. 1986-92 Numerical Recipes Software Y"1p3-1.
      SUBROUTINE asolve(sa,ija,n,b,x,itrnsp)
c     REAL*8,DIMENSION(:),allocatable::sa
c     INTEGER*4,DIMENSION(:),allocatable::ija
c     INTEGER n,itrnsp,i,ija(*),NTRAJ
      INTEGER n,itrnsp,i,ija(*)
      real*8 x(n),b(n),sa(*)
c     PARAMETER (NTRAJ=6000)
c     COMMON /mat/ sa,ija
      do 11 i=1,n
        x(i)=b(i)/sa(i)
11    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software Y"1p3-1.
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE dsprstx(sa,ija,x,b,n)
      INTEGER n,ija(*)
      real*8 b(n),sa(*),x(n)
      INTEGER i,j,k
      if (ija(1).ne.n+2) then
          write(*,*)'IJA(1) n+2',ija(1),n+2
          pause 'mismatched vector and matrix in sprstx'
      endif
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
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE dsprsax(sa,ija,x,b,n)
      INTEGER n,ija(*)
      real*8 b(n),sa(*),x(n)
      INTEGER i,k
      if (ija(1).ne.n+2) then
          write(*,*)'IJA(1) n+2',ija(1),n+2
          pause 'mismatched vector and matrix in sprsax'
      endif
      do 12 i=1,n
        b(i)=sa(i)*x(i)
        do 11 k=ija(i),ija(i+1)-1
          b(i)=b(i)+sa(k)*x(ija(k))
11      continue
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software Y"1p3-1.
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      FUNCTION snrm(n,sx,itol)
      INTEGER n,itol,i,isamax
      real*8 sx(n),snrm
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
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!            INTEGER, PARAMETER :: stack_size = 1000
!     INTEGER, SAVE :: store(stack_size), pos = 0
!           CONTAINS
                                                                                                
      SUBROUTINE pushf(store,i,pos)
             IMPLICIT NONE
             INTEGER, PARAMETER :: stack_size = 1000
             INTEGER, INTENT(IN) :: i
             INTEGER :: store(stack_size),pos
c       common/stackk/store(stack_size)
                                                                                                
      IF (pos < stack_size) THEN
       pos = pos + 1; store(pos) = i
      ELSE
       STOP 'Stack Full error'
      END IF
        return
                                                                                                
        END SUBROUTINE pushf
                                                                                                
        SUBROUTINE popf(store,i,pos)
             IMPLICIT NONE
             INTEGER, PARAMETER :: stack_size = 1000
             INTEGER, INTENT(OUT) :: i
             INTEGER :: store(stack_size),pos
c       common/stackk/store(stack_size)
                                                                                                
      IF (pos > 0) THEN
       i = store(pos); pos = pos - 1
      ELSE
        i=0
!      STOP 'Stack Empty error'
      END IF
        return
        END SUBROUTINE popf
c-----------------------------------------------------------------
C replacement for a function available in many fortran libraries
C find the length of the string without trailing blanks
             function lnblnk(string)
             character*(*) string
             nchars=len(string)
             do i=nchars,1,-1
                   if(string(i:i).ne.' ') then
C                  if(string(i:i+1).ne.' ') then
                       lnblnk=i
                       return
                    endif
             enddo
             lnblnk=0
             return
             end
c---------------------------------------------------------
