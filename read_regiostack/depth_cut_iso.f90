! Test program for reading regiostack model

      program depth_cut

      implicit none
      include 'constants.h'

      double precision depth, r, theta, phi
      double precision dvs, drho, dvp
      integer i,j
      
      call read_mantle_model  

      write(*,*) 'Input depth for cut (km) : '
      read(*,*) depth

      r=1-depth*1000.0d0/R_EARTH

      open(unit=10,file='cut.dat',status='unknown') 
  
      do j=1,360
        theta=j*DEGREES_TO_RADIANS/2.0d0
        do i=1,720
          phi=i*DEGREES_TO_RADIANS/2.0d0
          call mantle_model(r,theta,phi,dvs,dvp,drho)
          write(10,'(f8.2,1x,f8.2,1x,f8.2)') 90-j/2.0,i/2.0,dvs*100
        enddo
      enddo
      

      end program depth_cut
