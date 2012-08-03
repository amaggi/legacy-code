! Test program for reading regiostack model

      program depth_cut

      include 'constants.h'

      integer, parameter :: NORMMASK=2575

      integer, parameter :: NLAT = 180
      integer, parameter :: NLON = 360

      double precision depth, r, theta, phi
      double precision dvs, alpha, an_phi, ampG
      integer i,j

      real olat(NORMMASK), olon(NORMMASK), bathy, age
      real orm

      
      call read_mantle_model  

      write(*,*) 'Input depth for cut (km) : '
      read(*,*) depth

      r=1-depth*1000.0d0/R_EARTH


!     calculate orm
      open(unit=10,file='/home/alessia/share/ORM/mask_30_70_pacific.xy',status='old')
      do i=1,NORMMASK
        read(10,*) olon(i), olat(i), bathy, age
        theta=(90-olat(i))*DEGREES_TO_RADIANS
        phi=olon(i)*DEGREES_TO_RADIANS
        call mantle_model(r,theta,phi,dvs,alpha,an_phi,ampG)
        orm=orm+dvs
      enddo
      close(10)
      orm=orm/NORMMASK
      write(*,*) 'ORM at ',depth,'km depth = ',orm


!      write(*,*) 'Lat, lon of point of interest :'
!      read(*,*) theta, phi
! 
!      theta=(90-theta)*DEGREES_TO_RADIANS
!      phi=phi*DEGREES_TO_RADIANS
!
!      call mantle_model(r,theta,phi,dvs,alpha,an_phi,ampG)
!      dvs=(dvs/orm-1)*100
!      write(*,'(f8.2,1x,f8.2,1x,f8.4,1x,f8.4,1x,f8.2,1x,f8.4)') &
!      90-theta*RADIANS_TO_DEGREES,phi*RADIANS_TO_DEGREES,dvs,&
!      alpha*100,90-an_phi*RADIANS_TO_DEGREES,ampG
  
      open(unit=10,file='cut_aniso.dat',status='unknown') 
      do j=1,NLAT
        theta=j*DEGREES_TO_RADIANS
        do i=1,NLON
          phi=i*DEGREES_TO_RADIANS
          call mantle_model(r,theta,phi,dvs,alpha,an_phi,ampG)
          dvs=(dvs/orm-1)*100
          write(10,'(f8.2,1x,f8.2,1x,f8.4,1x,f8.4,1x,f8.2,1x,f8.4)') &
           90-theta*RADIANS_TO_DEGREES,phi*RADIANS_TO_DEGREES,dvs,&
           alpha*100,90-(an_phi*RADIANS_TO_DEGREES),ampG
        enddo
      enddo

      open(unit=11,file='anis_xy.dat',status='unknown') 
      open(unit=10,file='cut_aniso_2.dat',status='unknown') 
      do i=1,633
!      do i=1,2856
        read(11,*) theta, phi
        !write(*,*) theta, phi
        theta=(90-theta)*DEGREES_TO_RADIANS
        phi=phi*DEGREES_TO_RADIANS
        call mantle_model(r,theta,phi,dvs,alpha,an_phi,ampG)
        dvs=(dvs/orm-1)*100
        write(10,'(f8.2,1x,f8.2,1x,f8.4,1x,f8.4,1x,f8.2,1x,f8.4)') &
        90-theta*RADIANS_TO_DEGREES,phi*RADIANS_TO_DEGREES,dvs,&
        alpha*100,90-(an_phi*RADIANS_TO_DEGREES),ampG
      enddo
      close(10)
      close(11)


      end program depth_cut
