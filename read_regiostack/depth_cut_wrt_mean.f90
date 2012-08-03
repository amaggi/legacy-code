! Test program for reading regiostack model

      program depth_cut

      implicit none
      include 'constants.h'

      integer, parameter :: NLAT = 180
      integer, parameter :: NLON = 360

      double precision depth, r, theta, phi, sumweight
      double precision dvs(NLAT,NLON), drho, dvp, dvsmean
      integer i,j
      
      call read_mantle_model  

      write(*,*) 'Input depth for cut (km) : '
      read(*,*) depth

      r=1-depth*1000.0d0/R_EARTH

      open(unit=10,file='cut.dat',status='unknown') 
  
      dvsmean=0
      sumweight=0
      do j=1,NLAT
        theta=j*DEGREES_TO_RADIANS*1.0d0
        do i=1,NLON
          phi=i*DEGREES_TO_RADIANS*1.0d0
          call mantle_model(r,theta,phi,dvs(j,i),dvp,drho)
          dvsmean=dvsmean+dvs(j,i)*sin(theta)
          sumweight=sumweight+sin(theta)
        enddo
      enddo

      dvsmean=dvsmean / sumweight
      
      do j=1,NLAT
        do i=1,NLON
          dvs(j,i)=dvs(j,i)-dvsmean
          write(10,'(f8.2,1x,f8.2,1x,f8.2)') 1.0d0*(90-j),1.0d0*i,dvs(j,i)*100
        enddo
      enddo

      end program depth_cut
