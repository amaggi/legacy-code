! Test program for reading regiostack model

      program depth_cut

      implicit none
      include 'constants.h'

      double precision mindepth, maxdepth, depth, dstep, r, theta, phi
      double precision dvs, drho, dvp
      double precision lon, lat, dist
      integer i,npts
      
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
          call mantle_model(r,theta,phi,dvs,dvp,drho)
          write(10,'(5(f8.2,1x))') lon, lat, depth, dist, dvs*100
          depth=depth+dstep
        enddo
      enddo

      close(10)
      close(11)
      

      end program depth_cut
