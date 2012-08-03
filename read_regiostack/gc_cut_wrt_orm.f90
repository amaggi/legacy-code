! Test program for reading regiostack model

      program depth_cut

      include 'constants.h'

      parameter (MAX_DEPTHS=1000,NORMMASK=2575)

      double precision mindepth, maxdepth, depth, dstep, r, theta, phi
      double precision dvs, drho, dvp
      double precision orm(MAX_DEPTHS)
      integer i,npts,ndepth

      real olat(NORMMASK), olon(NORMMASK), bathy, age
      real lat,lon,dist
      
      call read_mantle_model  


      write(*,*) 'Number of points : '
      read(*,*) npts
      write(*,*) 'Input start depth for cut (km) : '
      read(*,*) mindepth
      write(*,*) 'Input end depth for cut (km) : '
      read(*,*) maxdepth
      write(*,*) 'Depth step (km) :'
      read(*,*) dstep

      ndepth=(maxdepth-mindepth)/dstep

!     sanity check
      if (ndepth.gt.MAX_DEPTHS) then 
        stop 'Increase MAX_DEPTHS'
      endif

!     initialise orm
      do k=1,ndepth
        orm(k)=0
      enddo

!     calculate and output orm
      open(unit=10,file='/home/alessia/share/ORM/mask_30_70_pacific.xy',status='old')
      open(unit=12,file='orm.dat',status='unknown') 
      do i=1,NORMMASK
        read(10,*) olon(i), olat(i), bathy, age
        theta=(90-olat(i))*DEGREES_TO_RADIANS
        phi=olon(i)*DEGREES_TO_RADIANS
        do k = 1,ndepth
          depth=mindepth+(k-1)*dstep
          r=1-depth*1000.0d0/R_EARTH
          call mantle_model(r,theta,phi,dvs,dvp,drho)
          orm(k)=orm(k)+dvs
        enddo
      enddo
      do k=1,ndepth
        orm(k)=orm(k)/NORMMASK
        depth=mindepth+(k-1)*dstep
        write(12,'(f8.2,1x,f8.4)') depth, orm(k)
      enddo
      close(10)
      close(12)

!     get the absolute values of vs and turn them into
!     percentage variations wrt orm
!     also write out both vs profile and orm
      open(unit=10,file='cut.dat',status='unknown') 
      open(unit=11,file='track.xyz',status='old')
      do i=1,npts
        read(11,*) lon, lat, dist
        theta=(90-lat)*DEGREES_TO_RADIANS
        phi=lon*DEGREES_TO_RADIANS
        do k=1,ndepth
          depth=mindepth+(k-1)*dstep
          r=1-depth*1000.0d0/R_EARTH
          call mantle_model(r,theta,phi,dvs,dvp,drho)
!         turn absolute value of vs into perc wrt orm
          dvs=(dvs/orm(k)-1)*100
          write(10,'(5(f8.2,1x))') lon, lat, depth, dist, dvs
        enddo
      enddo
      close(10)
      close(11)

      end program depth_cut
