! Test program for reading regiostack model

      program depth_cut

      implicit none
      include 'constants.h'

      integer, parameter :: NLAT = 180
      integer, parameter :: NLON = 360


      double precision depth, r, theta, phi
      double precision dvs, alpha, an_phi, ampG
      integer i,j
      
      call read_mantle_model  

      write(*,*) 'Input depth for cut (km) : '
      read(*,*) depth
      !depth=50

      r=1-depth*1000.0d0/R_EARTH

      open(unit=10,file='cut_aniso.dat',status='unknown') 
  
      do j=1,NLAT
        theta=j*DEGREES_TO_RADIANS
        do i=1,NLON
          phi=i*DEGREES_TO_RADIANS
!          print *, '(Th,Ph) = (',theta*RADIANS_TO_DEGREES, &
!                   phi*RADIANS_TO_DEGREES,')'
          call mantle_model(r,theta,phi,dvs,alpha,an_phi,ampG)
          write(10,'(f8.2,1x,f8.2,1x,f8.4,1x,f8.4,1x,f8.2,1x,f8.4)') &
           90-theta*RADIANS_TO_DEGREES,phi*RADIANS_TO_DEGREES,dvs,&
           alpha*100,90-(an_phi*RADIANS_TO_DEGREES),ampG
        enddo
      enddo

      open(unit=11,file='anis_xy.dat',status='unknown') 
      open(unit=10,file='cut_aniso_2.dat',status='unknown') 
      do i=1,633
        read(11,*) theta, phi
        !write(*,*) theta, phi
        theta=(90-theta)*DEGREES_TO_RADIANS
        phi=phi*DEGREES_TO_RADIANS
        call mantle_model(r,theta,phi,dvs,alpha,an_phi,ampG)
!        dvs=(dvs/orm-1)*100
        write(10,'(f8.2,1x,f8.2,1x,f8.4,1x,f8.4,1x,f8.2,1x,f8.4)') &
        90-theta*RADIANS_TO_DEGREES,phi*RADIANS_TO_DEGREES,dvs,&
        alpha*100,90-(an_phi*RADIANS_TO_DEGREES),ampG
      enddo
      close(10)
      close(11)


      end program depth_cut
