! Test program for reading regiostack model

      program depth_cut

      include 'constants.h'

      parameter (NORMMASK=2575)

      integer, parameter :: NLAT = 180
      integer, parameter :: NLON = 360

      double precision depth, r, theta, phi
      double precision dvs, drho, dvp

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
        call mantle_model(r,theta,phi,dvs,dvp,drho)
        orm=orm+dvs
      enddo
      close(10)
      orm=orm/NORMMASK
      write(*,*) 'ORM at ',depth,'km depth = ',orm


      open(unit=10,file='cut.dat',status='unknown') 
  
      do j=1,NLAT
        theta=j*DEGREES_TO_RADIANS*1.0d0
        do i=1,NLON
          phi=i*DEGREES_TO_RADIANS*1.0d0
          call mantle_model(r,theta,phi,dvs,dvp,drho)
          dvs=(dvs/orm-1)*100
          write(10,'(f8.2,1x,f8.2,1x,f8.2)') 1.0d0*(90-j),1.0d0*i,dvs
        enddo
      enddo

      end program depth_cut
