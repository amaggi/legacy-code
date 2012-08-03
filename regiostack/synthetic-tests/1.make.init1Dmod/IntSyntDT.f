c sylvana: added perturb. magnitude in input parameters      
c alessia: modified to do checkerboard test
c alessia: modified to take V3 formatted intomodes files
c alessia: modified to take both PREM and 3SMAC input
c alessia: removed cross-section output (can do cross-sections directly from
c          vsfin and anvsfin files)
c alessia: modified to do DT calculations instead of Vs
c
c $Id: IntSyntV3.f,v 1.4 2004/06/16 13:44:14 alessia Exp $
c
      program IntSyntV3

      include 'param_V3.h'

      parameter(NPTSMAX=18000,NCOUPMAX=20, NPERTMAX=100)
      parameter(NTMAX=3)

      real lat, lon

c     Grid and integration steps
      real gridstep, intstep

c     bg model type flag
      integer model_type
c     integer smac50
      real u_dt

c     3D arrays for isotropic velocity and anisotropy
      real c(NPROFMAX,NLIN,NCOL)
      real a1(NPROFMAX,NLIN,NCOL),a2(NPROFMAX,NLIN,NCOL)

c     Intomodes parameters
      real late(NBTRAMAX),lone(NBTRAMAX) 
      real lats(NBTRAMAX),lons(NBTRAMAX) 
      real err(NBTRAMAX,NPROFMAX) 
      character*80 nomtra(NBTRAMAX)
      character*40 nomtrajet

c     used for 1D models
      real sl(NPROFMAX,NPTSMAX)
      real tout(NPROFMAX)
      real cref(NPROFMAX)
      real ssom(NPROFMAX)
      real s(NPTSMAX)

c     1D depth array
      real prof(NPROFMAX)

c     Perturbation type flags
c     integer anis_type, synt_type
      integer synt_type

c     For cylindrical perturbations
      integer n_cyl
      real cyl_lat(NPERTMAX), cyl_lon(NPERTMAX)
      real cyl_mag(NPERTMAX), cyl_rad(NPERTMAX)
      real cyl_minprof(NPERTMAX), cyl_maxprof(NPERTMAX)

c     For plume perturbations
c     integer n_plume
c     real plume_dT(NPERTMAX), plume_rad(NPERTMAX)
c     real plume_lat(NPERTMAX), plume_lon(NPERTMAX)
c     real plume_minprof(NPERTMAX), plume_maxprof(NPERTMAX)
c     real plume_VsTgrad(NPROFMAX,NTMAX) 
c     real plume_profs(NPROFMAX), plume_T(NTMAX)
c     integer n_adb, n_plume_profs

c     For checkerboard perturbations
      real mag, pert_size, profmin, profmax
c     integer check_anis

c     For lithospheric cooling models
c     real age(NLIN,NCOL)
c     integer cooling_type
c     real Tm, kappa, plate_th, dlogVsdT

      integer isor

      common/c1/nbtra,ncouch,nlat,nlon

      isor=12
      pi=3.1415926
      degtorad=pi/180.0
      radtodeg=180.0/pi      


c     ask about number of periods and periods:
      write(*,*) 'Number of periods covered by intomo file:'
      read(*,*) ncouch
      write(*,*) 'List the ',ncouch,' periods:'
      read(*,*) (prof(i),i=1,ncouch)

c     read intomodes file
      write(*,*)'Lecture du fichier intomodes'
      call lecintomo(lats,lons,nomtra,late,lone,err)
      call write_coverage(lats,lons,late,lone,nbtra)

c     what kind of starting model
      write(*,*) 'Starting model: '
      write(*,*) '1. Uniform 0 residual slowness '
      write(*,*) '2. Uniform non zero residual slowness '
      read(*,*)model_type
      if (model_type.eq.2) then
        write(*,*) 'Value of non-zero residual slowness :'
        read(*,*) u_dt
      endif

c     what kind of errors?
      write(*,*) 'Use data errors? [1=yes]'
      read(*,*) ierror
      if (ierror.ne.1) then
        write(*,*) 'Fixed data error value (in s : e.g. 1.0)'
        read(*,*) fixederr
      endif

c     ask for required model step:
      write(*,*)'Geometry: model step in degrees: '
      read(*,*) gridstep
      intstep=gridstep/2.
      nlat=int(180/gridstep)
      nlon=int(360/gridstep)

c     ask for synthetic type
      write(*,*) 'Type of synthetic test desired:'
      write(*,*) '0: No perturbation'
      write(*,*) '1: Checkerboard'
      write(*,*) '2: Cylinders'
c     write(*,*) '3: Temperature defined plume conduits'
c     write(*,*) '4: Anisotropy perturbation only'
c     write(*,*) '5: Lithospheric cooling model'

      read(*,*) synt_type

c     -----------------------------------------------------------------
c     Read more paramters according to each synthetic type
c     -----------------------------------------------------------------
c     -----------------------------------------------------------------
c     cooling model signature
c     -----------------------------------------------------------------
      if (synt_type .eq. 5 ) then
c       write(*,*) 'Type of cooling signature:'
c       write(*,*) '(1) Half space cooling'
c       write(*,*) '(2) Plate model'
c       read(*,*) cooling_type

c       standard parameters for all cooling models
c       write(*,*) 'Mantle temperature Tm (K):'
c       read(*,*) Tm
c       write(*,*) 'Mantle diffusivity Kappa (km^2/Ma):'
c       read(*,*) kappa
c       write(*,*) 'dlovVsdT (%/100K):'
c       read(*,*) dlogVsdT

c       parameters for plate models
c       if (cooling_type .eq. 2) then
c         write(*,*) 'Asymptotic plate thickness (km):'
c         read(*,*) plate_th
c       end if

c       read age file
c       call read_age(age,gridstep)
c     -----------------------------------------------------------------

c     -----------------------------------------------------------------
c     simple anisotropy pattern
c     -----------------------------------------------------------------
      elseif (synt_type .eq. 4 ) then
c       write(*,*) 'Magnitude of azimuthal anisotropy will be 2%.'
c       write(*,*) 'Anisotropy type:'
c       write(*,*) '(1) Uniform NS'
c       write(*,*) '(2) Uniform EW'

c       read(*,*)anis_type
c     -----------------------------------------------------------------

      else if (synt_type .eq. 3 ) then
c     -----------------------------------------------------------------
c     plume conduits
c     -----------------------------------------------------------------
c       print *,'Number of plume perturbations ',
c    *         '(max ', NPERTMAX,')'
c       read(*,*)n_plume
c       if (n_plume.gt.NPERTMAX) then
c         write(*,*) 'Too many perturbations (max = ',NPERTMAX,')'
c         stop
c       endif
c       do 222 i=1,n_plume
c         print *,'Lat, lon, radius(km), mindepth, maxdepth, ',
c    *            'magnitude (delta T in Kelvin) for plume', i
c         write(*,*)'(e.g. 25, -15, 200, 100, 410, 300)'
c         read(*,*)plume_lat(i),plume_lon(i),plume_rad(i),
c    *             plume_minprof(i),plume_maxprof(i),plume_DT(i)
c22     continue
c       hard code values for cammarano file
c       n_adb=3
c       n_plume_profs=27
c       call setup_VsTgrad(plume_VsTgrad,plume_profs,plume_T,
c    *                     n_plume_profs,n_adb)
c     -----------------------------------------------------------------

      elseif (synt_type .eq. 2 ) then
c     -----------------------------------------------------------------
c     cylinders
c     -----------------------------------------------------------------
      print *,'Number of cylindrical perturbations ',
     *         '(max ', NPERTMAX,')'
       read(*,*)n_cyl
       if (n_cyl.gt.NPERTMAX) then
         write(*,*) 'Too many perturbations (max = ',NPERTMAX,')'
         stop
       endif
       do 2 i=1,n_cyl
         print *,'Lat, lon, radius(km), mindepth, maxdepth, ',
     *           'magnitude (s/km) for cylinder', i
         write(*,*)'(e.g. 25, 115, 200, 50, 410, 0.03)'
         read(*,*)cyl_lat(i),cyl_lon(i),cyl_rad(i),cyl_minprof(i),
     *            cyl_maxprof(i),cyl_mag(i)
         if (cyl_lon(i).lt.0) cyl_lon(i) = 360 + cyl_lon(i) 
2      continue
c     -----------------------------------------------------------------

      elseif (synt_type.eq.1) then
c     -----------------------------------------------------------------
c     checkerboard
c     -----------------------------------------------------------------
c       ask for checkerboard parameters
        write(*,*) 'Lateral size of chekerboard (km)'
        write(*,*)'(e.g. 1000) '
        read(*,*) pert_size
        write(*,*)'Prof min et max de la perturbation? '
        write(*,*)'(e.g. 300 450) '
        read(*,*) profmin,profmax
        write(*,*)'Magnitude (s/km) of perturbation? (e.g. 0.03)'
        read(*,*) mag
c       write(*,*)'Anisotropic perturbation? (0=no, 1=yes)'
c       read(*,*) check_anis
c     -----------------------------------------------------------------
      endif


c     -----------------------------------------------------------------
c     END USER INPUT
c     -----------------------------------------------------------------

c-----------------------------------------------
c     SET STARTING MODEL
c-----------------------------------------------
      if (model_type.eq.1) then
        call set_uniform(c,0.0)
      else if (model_type.eq.2) then
        call set_uniform(c,u_dt)
      endif

c-----------------------------------------------
c     ADD PERTURBATIONS
c-----------------------------------------------
      write(*,*) 'Adding perturbations'

c------------------------------------------------------------------
c    start loop over depth layers
c------------------------------------------------------------------
      do 10 jp=1,ncouch
        write(*,*) 'Doing layer at ',prof(jp), 's'

        cref(jp)=0
        do 12 i=1,nlat
        do 12 j=1,nlon
c         find the latitude and longitude of this point
          lat=getlat(i,gridstep)
          lon=getlon(j,gridstep)

c-----------------------------------------------------------------
c         apply the relevant perturbation
c-----------------------------------------------------------------
c-----------------------------------------------------------------
c           cooling models
c-----------------------------------------------------------------
          if(synt_type .eq. 5) then
c           if(age(i,j) .lt. 250) then
c           if (cooling_type .eq. 1) then
c             c(jp,i,j)=halfspcool_vs(age(i,j),c(jp,i,j),prof(jp),
c    &                                dlogVsdT,Tm,kappa)
c           else if (cooling_type .eq. 2) then
c             c(jp,i,j)=platecool_vs(age(i,j),c(jp,i,j),prof(jp), 
c    &                               dlogVsdT,Tm,kappa,plate_th)
c           endif
c           endif

c-----------------------------------------------------------------
c           anisotropy
c-----------------------------------------------------------------
          elseif(synt_type .eq. 4) then
c           if(anis_type.eq.1) then
c             Uniform NS anisotropy 2%
c             a1(jp,i,j)= 0.045
c             a2(jp,i,j)= 0.0
c           elseif(anis_type.eq.2) then
c             Uniform EW anisotropy 2%
c             a1(jp,i,j)= -0.045   
c             a2(jp,i,j)= 0.0
c           endif
          elseif(synt_type .eq. 3) then
c-----------------------------------------------------------------
c           plumes
c-----------------------------------------------------------------
c           no point doing any of this outside our gradient bounds
c           i.e. 100-800 km depth
c           if (prof(jp).ge.100 .and. prof(jp).le.800) then
c             plume_grad=VsTgrad1300(prof(jp),plume_VsTgrad,
c    *                   plume_profs,plume_T,n_plume_profs,n_adb)
c             c(jp,i,j)=vs_plume(c(jp,i,j),lat,lon,prof(jp),plume_lat,
c    *                  plume_lon,plume_rad,plume_maxprof,
c    *                  plume_minprof,plume_DT,n_plume,plume_grad)
c           endif

c         elseif(synt_type .eq. 2) then
c-----------------------------------------------------------------
c           cylinders
c-----------------------------------------------------------------
            c(jp,i,j)=dt_cyl(c(jp,i,j),lat,lon,prof(jp),
     *                        cyl_lat, cyl_lon, cyl_rad,
     *                        cyl_maxprof, cyl_minprof,
     *                        cyl_mag,n_cyl)
          else
c-----------------------------------------------------------------
c           checkerboard
c-----------------------------------------------------------------
            c(jp,i,j)=dtchecker(c(jp,i,j),lat,lon,prof(jp),pert_size,
     *                           profmin, profmax,mag)
c           if (check_anis .eq. 1) then
c             a1(jp,i,j)=0.0
c             a2(jp,i,j)=a2checker(lat,lon,prof(jp),pert_size,
c    *                           profmin, profmax)
c           endif
          endif
c-----------------------------------------------------------------
c         finished applying perturbations to this depth
c-----------------------------------------------------------------
c         add to reference velocity sum
          cref(jp)=cref(jp)+c(jp,i,j)
12      continue
c       reference velocity for depth jp is the mean of the velocities
        cref(jp)=cref(jp)/nlat/nlon

c     do next depth
10    continue
c-----------------------------------------------
c     END PERTURBATIONS
c-----------------------------------------------

c---------------------------------------------------
c     write input anvsfin and vsfin
c---------------------------------------------------
      write(*,*) 'Writing input vsfin / anvsfin files'
      call writefin(prof,c,a1,a2,ncouch,nlat,nlon)

c-----------------------------------------------
c     WRITE INTOMODES FILE
c-----------------------------------------------
      write(*,*) 'Writing intomodes file'
      open(isor,file='SYNTintomo')
      
      do 25 jp=1,ncouch
25    ssom(jp)=0

      do 30 itra=1,nbtra
c       change name of path to indexed station name e.g. KIPP.z14 
c       (guaranteed unique even for clustered paths)
        nomtrajet=nomtra(itra)(:lnblnk(nomtra(itra)))
        write(isor,1002) nomtrajet
        write(isor,*)late(itra),lone(itra),lats(itra),lons(itra)
1002    format(a40)

c       if(itra.eq.802) write(*,*)lats(itra),lons(itra)

        call mygrt(late(itra),lone(itra),lats(itra),lons(itra),
     *             disd,az12,az21,gc)
        delta=disd
        azes=az12
        do while(azes.lt.0.) 
          azes=azes+360.
        enddo
        do while(azes.gt.360.) 
          azes=azes-360.
        enddo
        dpas= intstep
        npas= int(delta/intstep)
        if (npas.eq.0) stop 'Stepping interval smaller than path length'
        npat= npas+1

c  Boucle sur chaque point de l'azimuth epicentre-station
c
        do 60 jcle=1,npat 
          if (jcle.eq.npat) then
c           special case for the last point: make sure it is the actual
c           end of the path
            blat=lats(itra)
            blon=lons(itra)
          else
c           find next lat lon point along the path
            dstcle=dpas*float(jcle-1)*111.195 
            call mygds(late(itra),lone(itra),azes,dstcle,blat,blon)
          endif
c         find the azimuth and backazimuth for the new point
          call mygrt(late(itra),lone(itra),blat,blon,dstcle,azim,
     *               azback,gc)
c         local azimuth of the path is the local backazimuth towards
c         the event minus 180 degrees
          azb=azback-180.
          do while (azb.lt.0.) 
            azb=azb+360.
          enddo
          do while (azb.gt.360.) 
            azb=azb-360.
          enddo

          i=int((90-blat)/gridstep)+1
          if(blon.lt.0.)blon=360.+blon
          j=int(blon/gridstep)+1

          do 70 jp=1,ncouch 
            sl(jp,jcle)=c(jp,i,j)
c           sl(jp,jcle)=slv(jp,jcle)+a1(jp,i,j)*cos(2*degtorad*azb)+
c    *                 a2(jp,i,j)*sin(2*degtorad*azb)
c           if(sl(jp,jcle).eq.0) then
c             write(*,*)'V NUL TRA',itra
c             stop
c           endif
70        continue 
 
60      continue 

c  Il ne reste plus qu'a integrer sur le trajet.
        do 80 jp=1,ncouch
          do 90 jcle=1,npat 
            s(jcle)=sl(jp,jcle)
90        continue
          call trapez(npat,s,dpas,sint)
c         integral over residual slowness IS delay time
          tout(jp)=sint
80      continue

c       add to sum of average slownesses
        do 120 ii=1,ncouch
          ssom(ii)=ssom(ii)+tout(ii)/(dpas*(npat-1))
120     continue

        write(isor,1003) (tout(ii),ii=1,ncouch)
        if(ierror.eq.1) then
          write(isor,1003) (err(itra,ii),ii=1,ncouch)
        else
          write(isor,1003) (fixederr,ii=1,ncouch)
        endif

30    continue

c  Affichage de la moyenne des lenteurs pour l'apriori.

      do 130 ii=1,ncouch
        write(*,*) prof(ii),ssom(ii)/nbtra,' .05 0. 0. 0.005'
130   continue

      close (isor)
1003  format(50f8.4)


      return

      end program



