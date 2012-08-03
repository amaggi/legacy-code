! Test program for reading regiostack model

      program depth_cut

      implicit none
      include 'constants.h'

      double precision mindepth, maxdepth, depth, dstep, r, theta, phi
      double precision dvs, drho, dvp
      double precision lon, lat, dist, moho
      integer i,j,npts
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
      open(unit=10,file='cut.dat',status='unknown') 
  
      do i=1,npts
        read(11,*) lon, lat, dist
        theta=(90-lat)*DEGREES_TO_RADIANS
        phi=lon*DEGREES_TO_RADIANS
        depth=mindepth
        do while (depth .le. maxdepth)
          r=1-depth*1000.0d0/R_EARTH
          call crustal_model(lat,lon,r,dvp,dvs,drho,moho,found_crust)
          if ( .not. found_crust ) then
            call mantle_model(r,theta,phi,dvs,dvp,drho)
          endif
          write(10,'(7(f8.2,1x))') lon, lat, depth, dist, dvs, drho
          depth=depth+dstep
        enddo
      enddo

      close(10)
      close(11)
      

      end program depth_cut
