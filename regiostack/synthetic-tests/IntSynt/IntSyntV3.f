c sylvana: added perturb. magnitude in input parameters      
c alessia: modified to do checkerboard test
c alessia: modified to take V3 formatted intomodes files
c alessia: modified to take both PREM and 3SMAC input
c alessia: removed cross-section output (can do cross-sections directly from
c          vsfin and anvsfin files)
c
c $Id: IntSyntV3.f,v 1.4 2004/06/16 13:44:14 alessia Exp $
c
      program IntSyntV3

      include 'param_V3.h'
      include 'slabs.h'

      parameter(NPTSMAX=18000,NCOUPMAX=20, NPERTMAX=100)
      parameter(NTMAX=3,NSLABMAX=20)
      parameter(PI=3.1415926)

      real lat, lon

c     Grid and integration steps
      real gridstep, intstep

c     bg model type flag
      integer model_type
      integer smac50

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
      real v(NPROFMAX,NPTSMAX)
      real vout(NPROFMAX)
      real cref(NPROFMAX)
      real vsom(NPROFMAX)
      real len(NPTSMAX)

c     1D depth array
      real prof(NPROFMAX)

c     Perturbation type flags
      integer anis_type, synt_type

c     For cylindrical perturbations
      integer n_cyl
      real cyl_lat(NPERTMAX), cyl_lon(NPERTMAX)
      real cyl_mag(NPERTMAX), cyl_rad(NPERTMAX)
      real cyl_minprof(NPERTMAX), cyl_maxprof(NPERTMAX)

c     For slab-like perturbations
      integer n_slabs
      character*256 slab_file(NSLABSMAX)
      real slab_pert

c     For plume perturbations
      integer n_plume
      real plume_dT(NPERTMAX), plume_rad(NPERTMAX)
      real plume_lat(NPERTMAX), plume_lon(NPERTMAX)
      real plume_minprof(NPERTMAX), plume_maxprof(NPERTMAX)
      real plume_VsTgrad(NPROFMAX,NTMAX) 
      real plume_profs(NPROFMAX), plume_T(NTMAX)
      integer n_adb, n_plume_profs

c     For checkerboard perturbations
      real mag, pert_size, profmin, profmax
      integer check_anis

c     For lithospheric cooling models
      real age(NLIN,NCOL)
      integer cooling_type
      real Tm, kappa, plate_th, dlogVsdT
      
c     For cylindrical plumes with anisotropy
      real plate_azimuth, scale_rad, mag_anis
      real plate_angle, theta_rot, anis_angle, anis_azimuth
      
      integer isor

      common/c1/nbtra,ncouch

c----------------------------------------------------------
c    25 JUIN 2001 : CODE ADAPTE A PARTIR DE IntsurtrajetsVs1anis.f
c    en vue du calcul d'un test synthetique complet
c    incluant waveform modelling et inversion 3D.
c
c    Programme de calcul de modeles synthetiques dans PREM + perturbation
c    Il calcul l'integrale des lenteurs sur un trajet donne.
c    les coordonnees des noms de trajets  et des epicentres
c    et stations sont dans un fichier format Montagner.
c
c    Cette version lit 
c        - Un model PREM 1D (des.in) extrapolle en chaque points geographique
c          ou aucune perturbation n'est ajoute 
c        - un fichier intomodesVs 
c          pour la couverture en trajet que l'on veut mettre dans le test.
c    Cette version permet 
c        - d'ajouter une anisotropie synthetique a 3SMAC.
c        - d'ajouter un dVs a PREM pour une region et un interval de profondeur donne 
c    
c----------------------------------------------------------
c     il1=12
      isor=12
      degtorad=pi/180.0
      radtodeg=180.0/pi      

c     read intomodes file
      write(*,*)'Number of depths?'
      read(*,*) ncouch
      write(*,*)'Lecture du fichier intomodesVs'
      call lecintomo(lats,lons,nomtra,late,lone,err)

c     what kind of starting model
      write(*,*) 'Starting model: '
      write(*,*) '1. PREM'
      write(*,*) '2. 3SMAC'
      write(*,*) '3. User defined model from vsfin / anvsfin files'
      read(*,*)model_type
      if (model_type.eq.2) then
        write(*,*) 'Use 50km 3SMAC for all depths? [1=yes]'
        read(*,*) smac50
      endif

c     what kind of errors?
      write(*,*) 'Use data errors? [1=yes]'
      read(*,*) ierror
      if (ierror.ne.1) then
        write(*,*) 'Fixed error value (e.g. 0.01)'
        read(*,*) fixederr
      endif

c     apply errors to epicenter locations?
      write(*,*) 'Apply errors to path averaged models ? [1=yes]'
      read(*,*) irand_vel
      if(irand_vel.eq.1) then
        write(*,*) 'Half-width of the gaussian error to apply?'
        read(*,*) rand_vel
      endif

c     ask for required model step:
      write(*,*)'Model step in degrees: '
      read(*,*) gridstep
      intstep=gridstep/2.
      nlat=int(180/gridstep)
      nlon=int(360/gridstep)
      
c-----------------------------------------------
c     READ STARTING VS MODEL
c-----------------------------------------------
      if (model_type.eq.1) then
        call lecPREM(c,a1,a2,gridstep,ncouch,prof)
      else if (model_type.eq.2) then
        call lec3smac(c,a1,a2,gridstep,ncouch,prof,smac50)
      else if (model_type.eq.3) then
        call readfin(prof,c,a1,a2,ncouch,nlat,nlon)
      endif


c     ask for synthetic type
      write(*,*) 'Type of synthetic test desired:'
      write(*,*) '0: No perturbation'
      write(*,*) '1: Checkerboard'
      write(*,*) '2: Cylinders'
      write(*,*) '3: Temperature defined plume conduits'
      write(*,*) '4: Anisotropy perturbation only'
      write(*,*) '5: Lithospheric cooling model'
      write(*,*) '6: Slab model'
      write(*,*) '7: Plumes with Kaminski Ribe 2002 anisotropy'

      read(*,*) synt_type

c     -----------------------------------------------------------------
c     Read more paramters according to each synthetic type
c     -----------------------------------------------------------------

c     -----------------------------------------------------------------
c     Kaminski Ribe 2002
c     -----------------------------------------------------------------
      if (synt_type .eq. 7) then
        write(*,*) 'Scale size of plumes (km) : '
        read(*,*) scale_rad
        write(*,*) 'Plate motion azimuth (degrees from N) : '
        read(*,*) plate_azimuth
        plate_azimuth=plate_azimuth*PI/180.0
        write(*,*) 'Azimuthal anisotropy magnitude (e.g. 0.02) : '
        read(*,*) mag_anis

c       get the plume locations
        print *,'Number of cylindrical perturbations ',
     *         '(max ', NPERTMAX,')'
        read(*,*)n_cyl
        if (n_cyl.gt.NPERTMAX) then
         write(*,*) 'Too many perturbations (max = ',NPERTMAX,')'
         stop
        endif
        do  i=1,n_cyl
         print *,'Lat, lon, magnitude Vs perturbation (fraction) ', i
         write(*,*)'(e.g. 25, 115, 0.03)'
         read(*,*)cyl_lat(i),cyl_lon(i),cyl_mag(i)
         if (cyl_lon(i).lt.0) cyl_lon(i) = 360 + cyl_lon(i) 
         cyl_minprof(i)=prof(1)
         cyl_maxprof(i)=prof(ncouch)
         cyl_rad(i)=scale_rad
       enddo
c       turn scale_rad into degrees
        scale_rad=scale_rad/111.25

c       calculate rotation angle for rotated coordinate system
  		theta_rot = -plate_azimuth 
  


c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     slab model signature
c     -----------------------------------------------------------------
      else if (synt_type .eq. 6) then
        write(*,*) 'Size of constant slab perturbation (e.g. 0.03) : '
        read(*,*) slab_pert
        write(*,*) 'Number of .slb files to read : '
        read(*,*) n_slabs

        if (n_slabs .gt. NSLABSMAX) then
          write(*,*) 'Too many slabs.  Max is ',NSLABSMAX
          stop
        endif

c       get the slab file names
        do i = 1, n_slabs
          write(*,*) 'Full path of .slb file ', i
          read(*,'(a)') slab_file(i)
        enddo
c     -----------------------------------------------------------------

c     -----------------------------------------------------------------
c     cooling model signature
c     -----------------------------------------------------------------
      else if (synt_type .eq. 5 ) then
        write(*,*) 'Type of cooling signature:'
        write(*,*) '(1) Half space cooling'
        write(*,*) '(2) Plate model'
        read(*,*) cooling_type

c       standard parameters for all cooling models
        write(*,*) 'Mantle temperature Tm (K):'
        read(*,*) Tm
        write(*,*) 'Mantle diffusivity Kappa (km^2/Ma):'
        read(*,*) kappa
        write(*,*) 'dlovVsdT (%/100K):'
        read(*,*) dlogVsdT

c       parameters for plate models
        if (cooling_type .eq. 2) then
          write(*,*) 'Asymptotic plate thickness (km):'
          read(*,*) plate_th
        end if

c       read age file
        call read_age(age,gridstep)
c     -----------------------------------------------------------------

c     -----------------------------------------------------------------
c     simple anisotropy pattern
c     -----------------------------------------------------------------
      elseif (synt_type .eq. 4 ) then
        write(*,*) 'Magnitude of azimuthal anisotropy will be 2%.'
        write(*,*) 'Anisotropy type:'
        write(*,*) '(1) Uniform NS'
        write(*,*) '(2) Uniform EW'

        read(*,*)anis_type
c     -----------------------------------------------------------------

      else if (synt_type .eq. 3 ) then
c     -----------------------------------------------------------------
c     plume conduits
c     -----------------------------------------------------------------
        print *,'Number of plume perturbations ',
     *         '(max ', NPERTMAX,')'
        read(*,*)n_plume
        if (n_plume.gt.NPERTMAX) then
          write(*,*) 'Too many perturbations (max = ',NPERTMAX,')'
          stop
        endif
        do 222 i=1,n_plume
          print *,'Lat, lon, radius(km), mindepth, maxdepth, ',
     *            'magnitude (delta T in Kelvin) for plume', i
          write(*,*)'(e.g. 25, -15, 200, 100, 410, 300)'
          read(*,*)plume_lat(i),plume_lon(i),plume_rad(i),
     *             plume_minprof(i),plume_maxprof(i),plume_DT(i)
222     continue
c       hard code values for cammarano file
        n_adb=3
        n_plume_profs=27
        call setup_VsTgrad(plume_VsTgrad,plume_profs,plume_T,
     *                     n_plume_profs,n_adb)
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
     *           'magnitude (percentage) for cylinder', i
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
        write(*,*)'Magnitude (%) of perturbation? (e.g. 0.05)'
        read(*,*) mag
        write(*,*)'Anisotropic perturbation? (0=no, 1=yes)'
        read(*,*) check_anis
c     -----------------------------------------------------------------
      endif


c     -----------------------------------------------------------------
c     END USER INPUT
c     -----------------------------------------------------------------


c-----------------------------------------------
c     ADD PERTURBATIONS
c-----------------------------------------------
      write(*,*) 'Adding perturbations'

c------------------------------------------------------------------
c    start loop over depth layers
c------------------------------------------------------------------
      do 10 jp=1,ncouch
        write(*,*) 'Doing layer at ',prof(jp), 'km'

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
c           Kaminski Ribe 2002
c-----------------------------------------------------------------
          	if(synt_type .eq. 7) then
c				apply isotropic plume perturbation
                c(jp,i,j)=vs_cyl(c(jp,i,j),lat,lon,prof(jp),
     *                        cyl_lat, cyl_lon, scale_rad,
     *                        prof(1), prof(ncouch),
     *                        cyl_mag,n_cyl)
     
           		anis_angle = kr_anis_angle(lat,lon,scale_rad,
     *							cyl_lat,cyl_lon,
     *                          n_cyl,plate_angle,theta_rot)
     
     			anis_azimuth = PI/2.0 - anis_angle

c           	apply anisotropy
           		a1(jp,i,j) = (mag_anis*c(jp,i,j)/2.0)*
     *      		         cos(2*anis_azimuth)
           		a2(jp,i,j) = (mag_anis*c(jp,i,j)/2.0)*
     *      		         sin(2*anis_azimuth)

c-----------------------------------------------------------------
c           cooling models
c-----------------------------------------------------------------
          else if(synt_type .eq. 5) then
            if(age(i,j) .lt. 250) then
            if (cooling_type .eq. 1) then
              c(jp,i,j)=halfspcool_vs(age(i,j),c(jp,i,j),prof(jp),
     &                                dlogVsdT,Tm,kappa)
            else if (cooling_type .eq. 2) then
              c(jp,i,j)=platecool_vs(age(i,j),c(jp,i,j),prof(jp), 
     &                               dlogVsdT,Tm,kappa,plate_th)
            endif
            endif

c-----------------------------------------------------------------
c           anisotropy
c-----------------------------------------------------------------
          elseif(synt_type .eq. 4) then
            if(anis_type.eq.1) then
c             Uniform NS anisotropy 2%
c             a1(jp,i,j)= 0.045
c             a2(jp,i,j)= 0.0
              a1(jp,i,j) = 0.02*c(jp,i,j)/2.0
              a2(jp,i,j) = 0
            elseif(anis_type.eq.2) then
c             Uniform EW anisotropy 2%
c              a1(jp,i,j)= -0.045   
c              a2(jp,i,j)= 0.0
              a1(jp,i,j) = -0.02*c(jp,i,j)/2.0
              a2(jp,i,j) = 0
            endif
          elseif(synt_type .eq. 3) then
c-----------------------------------------------------------------
c           plumes
c-----------------------------------------------------------------
c           no point doing any of this outside our gradient bounds
c           i.e. 100-800 km depth
            if (prof(jp).ge.100 .and. prof(jp).le.800) then
              plume_grad=VsTgrad1300(prof(jp),plume_VsTgrad,
     *                   plume_profs,plume_T,n_plume_profs,n_adb)
              c(jp,i,j)=vs_plume(c(jp,i,j),lat,lon,prof(jp),plume_lat,
     *                  plume_lon,plume_rad,plume_maxprof,
     *                  plume_minprof,plume_DT,n_plume,plume_grad)
            endif

          elseif(synt_type .eq. 2) then
c-----------------------------------------------------------------
c           cylinders
c-----------------------------------------------------------------
            c(jp,i,j)=vs_cyl(c(jp,i,j),lat,lon,prof(jp),
     *                        cyl_lat, cyl_lon, cyl_rad,
     *                        cyl_maxprof, cyl_minprof,
     *                        cyl_mag,n_cyl)
          else
c-----------------------------------------------------------------
c           checkerboard
c-----------------------------------------------------------------
            c(jp,i,j)=vschecker(c(jp,i,j),lat,lon,prof(jp),pert_size,
     *                           profmin, profmax,mag)
            if (check_anis .eq. 1) then
              a1(jp,i,j)=0.0
              a2(jp,i,j)=a2checker(lat,lon,prof(jp),pert_size,
     *                           profmin, profmax)
            endif
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
25    vsom(jp)=0

      inum=-204673
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
        npat= npas+1

c  Boucle sur chaque point de l azimuth epicentre-station
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
            v(jp,jcle)=c(jp,i,j)
            v(jp,jcle)=v(jp,jcle)+a1(jp,i,j)*cos(2*degtorad*azb)+
     *                 a2(jp,i,j)*sin(2*degtorad*azb)
c           v(jp,jcle)=v(jp,jcle)+a1(jp,i,j)*cos(2*degtorad*azes)+
c    *                 a2(jp,i,j)*sin(2*degtorad*azes)
            if(v(jp,jcle).eq.0) then
              write(*,*)'V NUL TRA',itra
              stop
            endif
70        continue 
 
60      continue 

c  Il ne reste plus qu a integrer sur le trajet.
    
        do 80 jp=1,ncouch
          do 90 jcle=1,npat 
            len(jcle)=1./v(jp,jcle)
90        continue
c         
          call trapez(npat,len,dpas,sint)
          if(delta.eq.0) then 
            write(*,*)'DELTA NUL TRA',itra
            stop
          endif
c         vout(jp)=sint/(delta)
          vout(jp)=sint/(dpas*(npat-1))
          if(vout(jp).eq.0)then
            write(*,*)'VOUT NUL TRA',itra
            stop
          endif
          vout(jp)=1./vout(jp)

c         if we are adding gaussian errors, then here is where we do so
          if (irand_vel.eq.1) then
            vout(jp)=vout(jp)+rand_vel*gasdev(inum)
          endif
80      continue

        do 120 ii=1,ncouch
          vsom(ii)=vsom(ii)+vout(ii)
120     continue

        write(isor,1003) (vout(ii),ii=1,ncouch)
        if(ierror.eq.1) then
          write(isor,1003) (err(itra,ii),ii=1,ncouch)
        else
          write(isor,1003) (fixederr,ii=1,ncouch)
        endif

30    continue

c  Affichage de la moyenne des lenteurs a 100 km utiliser comme modele
c  a priori
     

      do 130 ii=1,ncouch
        write(*,*) prof(ii),vsom(ii)/nbtra,' .05 0. 0. 0.005'
130   continue

      close (isor)
1003  format(50f8.4)


      return

      end program



