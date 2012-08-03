! Test program for reading regiostack model

      program path_average

      implicit none
      include 'constants.h'
      
      integer, parameter :: MAXDEP = 5000 ! something large enough

      double precision mindepth, maxdepth, depth, dstep, r, theta, phi
      double precision dvs, drho, dvp
      double precision lon, lat, dist, moho
      double precision vs_slowness(MAXDEP), density(MAXDEP), depths(MAXDEP)
      integer i,k,npts, ndepths
      logical found_crust
      
      call read_crustal_model  
      call read_mantle_model  


      write(*,*) 'Number of points : '
      read(*,*) npts
      write(*,*) 'Input start depth for cut (km) : '
      read(*,*) mindepth
      write(*,*) 'Input end depth for cut (km) : '
      read(*,*) maxdepth
      write(*,*) 'Depth step (km) :'
      read(*,*) dstep

      open(unit=11,file='track.xyz',status='old')
  
!     initialise path averaging quantities to zero
      depths(:) = 0 ; vs_slowness(:) = 0; density(:) = 0
      do i=1,npts
        read(11,*) lon, lat, dist
        theta=(90-lat)*DEGREES_TO_RADIANS
        phi=lon*DEGREES_TO_RADIANS
        depth=mindepth
        k=0
        do while (depth .le. maxdepth)
          k=k+1
          depths(k) = depth
          r=1-depth*1000.0d0/R_EARTH
          call crustal_model(lat,lon,r,dvp,dvs,drho,moho,found_crust)
          if ( .not. found_crust ) then
            call mantle_model(r,theta,phi,dvs,dvp,drho)
          endif
          density(k) = density(k) + drho
          vs_slowness(k) = vs_slowness(k) + 1/dvs
          depth=depth+dstep
        enddo
        ndepths=k
      enddo

      vs_slowness = vs_slowness / dble(npts)
      density = density / dble(npts)

      close(11)

      open(unit=10,file='path-average.dat',status='unknown') 
      do k = 1, ndepths
        write(10,'(3(f8.2,1x))') depths(k), 1/vs_slowness(k), density(k)
      enddo
      close(10)
      

      end program path_average
