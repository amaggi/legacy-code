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

        parameter (NTRAJ=40000,NSTEP=350,NDEPTH=50)
        parameter (NLAT=90,NLON=180)
        parameter (NTDEPTH=2*NDEPTH+8,KMAX=NLAT*NLON)
        parameter (NTR=NTRAJ)

c     DOUBLE PRECISION v(NTRAJ),vv(NTRAJ), uld(NTRAJ),res(NTRAJ)
      dimension delt(NTRAJ),dsi(NTRAJ),kdi(NTRAJ)
      character*4 staz1(NTRAJ)
      character*8 staz2(NTRAJ)

c     dimension iray(KMAX,NTR),iss(KMAX,NTR),ise(KMAX,NTR)

      dimension teta(KMAX),phi(KMAX)
      dimension nray(KMAX)


      dimension tetai(NSTEP,NTRAJ),phii(NSTEP,NTRAJ)
      dimension az(NSTEP,NTRAJ)
      dimension w(NTDEPTH,NTRAJ)

      dimension t(NDEPTH)

      character*32 ficdata
c     common/c1/tetai(NSTEP,NTRAJ),phii(NSTEP,NTRAJ)
c     common/c4/az
      common/c6/xo,yo,nxo,nyo,pxyo
      common/c61/xinf,yinf,nx,ny,pxy
c      common/c7/w(NTDEPTH,NTRAJ)

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
c
c
      print *,'enter data file name'
      read(5,'(a32)') ficdata
      print 999,ficdata
      open(10,file=ficdata,form='formatted')
       call lecdat(tetai,phii,az,w,staz1,staz2,delt,kdi,dsi,ntra,d1)
      close(10)

c  In lecdat, teta must be between 0. and teta and phi(radians)
c  covariance matrix for the whole paths is computed at each period
c  in order enable correlation length to be period dependent.
c  however, as this computation is quite long,it can be
c  computed at once if we forget this period dependency.
c  dcor, dcoran correlation lengths.
 999  format(a32)

      open(8,file='densite.xyz')

      do 50 j=1,KMAX
            phi(j)=0.
            teta(j)=0.
            nray(j)=0
 50   continue
 
      do 54 ii=1,nx
         do 54 jj=1,ny
          j2=ny*(ii-1)+jj
          teta(j2)=xinf+(ii-1)*pxy+pxy/2
          teta(j2)=teta(j2)*rad
          phi(j2)=yinf+(jj-1)*pxy+pxy/2
          phi(j2)=phi(j2)*rad
 54   continue

      jold=0
      ntrmax=0
      do 70 i=1,ntra
        tetai(1,i)=w(1,i)
        phii(1,i)=w(2,i)
        kki=kdi(i)

        do 70 ki=1,kki

           tx=tetai(ki,i)/rad
           py=phii(ki,i)/rad

c valable en tomo globale 
            ppy=py-(yinf+pxy*ny)
            if(ppy.le.10e-6.and.ppy.gt.0) then
               write(*,*)'BONJOUR ','TRAJET ',i,' POINT ',ki
               py=yinf+pxy/4.
            else if(ppy.ge.-10e-6.and.ppy.le.0) then
               py=py-pxy/4.
            endif

            ppy=py-yinf
            if(ppy.le.10e-6.and.ppy.gt.0) then
               write(*,*)'BONJOUR ','TRAJET ',i,' POINT ',ki
               py=yinf+pxy/4.
            else if(ppy.ge.-10e-6.and.ppy.le.0) then
               py=(yinf+pxy*ny)-pxy/4.
            endif
c modif 29/07/2002 :
c En tomo regionale lorsque la region
c chevauche la longitude +180 degres
c (i.e. Australie)
c il faut remettre py entre 0 et 360
c            if(py.lt.0) py=py+360.

           it2=int((tx-(xinf-pxy))/pxy)
           ip2=int((py-(yinf-pxy))/pxy)

           j2=ny*(it2-1)+ip2
           if(j2.gt.KMAX.or.j2.lt.1) then
                 write(*,*)'TRAJET ',i,' POINT ',ki
                 write(*,*)'TX ',tx,' PY ',py
                 write(*,*)'py-(yinf+pxy*ny)',py-(yinf+pxy*ny)
                 write(*,*)'J2 ',j2
                 stop' J2 out of range'
           endif


c          if(i.eq.nda.and.ki.eq.kki.and.j2.eq.jold)ise(j2,nray(j2))=ki

           if(j2.ne.jold.or.ki.eq.1) then
             nray(j2)=nray(j2)+1
             if(nray(j2).gt.ntrmax)ntrmax=nray(j2)
             if(nray(j2).gt.NTR)then
                 write(*,*)'TRAJET ',i,'POINT ',ki
                 write(*,*)'J2 ',j2,' nray(j2) ',nray(j2)
                 write(*,*)'TX ',tx,' PY ',py
                 write(*,*)'NTR ',NTR
                 stop' NTR must be increased'
             endif
c            iray(j2,nray(j2))=i
c            iss(j2,nray(j2))=ki
c            if(i.eq.nda.and.ki.eq.kki)then
c                write(*,*) 'DERNIER POINT DERNIER TRAJET'
c                ise(j2,nray(j2))=ki
c            endif
c            if(ki.ne.1)then
c                 ise(jold,nray(jold))=ki-1
c            else if(ki.eq.1.and.jold.ne.0.and.j2.ne.jold)then
c                 ise(jold,nray(jold))=kdi(i-1)
c            else if(ki.eq.1.and.jold.ne.0.and.j2.eq.jold)then
c                 ise(jold,nray(jold)-1)=kdi(i-1)
c            endif
             jold=j2
           endif

 70    continue

      write(*,*)'NTRMAX ',ntrmax
      do 99140 ii=1,nx
        tx=xinf+(ii-1)*pxy+pxy/2.
c       it1=int((tx-(xo-pxyo))/pxyo)

        do 99140 jj=1,ny
          py=yinf+(jj-1)*pxy+pxy/2.
c         ip1=int((py-(yo-pxyo))/pxyo)
          j2=ny*(ii-1)+jj

          write(8,*)90-tx,py,nray(j2)

99140  continue
       close(8)
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
      subroutine lecdat(tetai,phii,az,w,staz1,staz2,delt,kdi,
     *dsi,ndat,d1)
      
        parameter (NTRAJ=40000,NSTEP=350,NDEPTH=50)
        parameter (NSTAMAX=400)
        parameter (NLAT=90,NLON=180)
        parameter (NTDEPTH=2*NDEPTH+8,KMAX=NLAT*NLON)

      dimension delt(NTRAJ),dsi(NTRAJ),kdi(NTRAJ)
      dimension wlat(NSTAMAX),wlon(NSTAMAX)
      dimension t(NDEPTH),ug(NDEPTH),edat(NDEPTH)
      dimension tetai(NSTEP,NTRAJ),phii(NSTEP,NTRAJ)
      dimension az(NSTEP,NTRAJ)
      dimension deltei(NSTEP)
      dimension w(NTDEPTH,NTRAJ)
      character staz*13
      character*4 nomsta(NSTAMAX),staz1(NTRAJ),nom1,nom2
      character*8 staz2(NTRAJ)

c     common/c1/tetai(NSTEP,NTRAJ),phii(NSTEP,NTRAJ)
c     common/c4/az
c     common/c7/w(NTDEPTH,NTRAJ)

c   in this subroutine, we read the dataset(phase or group velocities)
c   the coordinates of the points along the paths are computed.

      rad=0.017453293
      print *,'nsta: number of stations'
      read(10,*)nsta
      print *,nsta

      do 88 i=1,nsta
         read(10,1998)nomsta(i),wlat(i),wlon(i)
         print 1998,nomsta(i),wlat(i),wlon(i)
  88  continue
 1998 format(a4,2x,2f9.4)

c     nt=nper
      nt=1
      do 1 i=1,ndat
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
c       write(*,*)'Lat lon epicenter, station',alate,alone,w(1,i),w(2,i)
c       MODIF eric pour avoir des lon entre 0 et 360 degres
c       if(alone.lt.0.)alone=alone+360
        w(1,i)=(90.-w(1,i))*rad
        w(2,i)=w(2,i)*rad
        in0=1
        fi=45.
        w(3,i)=(90.-alate)*rad
        w(4,i)=alone*rad
        w(5,i)=fi*rad
        w(6,i)=nt
        read(10,*)(ug(j),j=1,nt)
c       write(*,*)(ug(j),j=1,nt)
        read(10,*)(edat(j),j=1,nt)
        do 2 k=1,nt
          kk=k+6
          w(kk,i)=1/ug(k)
          kkk=k+NDEPTH+6
          w(kkk,i)=edat(k)/(ug(k)**2)
  2     continue
        in=1
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
          write(*,*)'TETA1',tetai(1,i),phii(1,i),deltei(1),az(1,i)
          write(*,*)'W',w(1,i)/rad,w(2,i)/rad
  1     continue
      return
      end
c----------------------------------------------------------
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
