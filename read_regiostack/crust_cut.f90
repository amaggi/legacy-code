! Test program for reading 3smac

      program crust_cut

      implicit none
      include 'constants.h'

      double precision mindepth, maxdepth, depth, dstep, r, theta, phi
      double precision vs, rho, vp, moho
      double precision lon, lat, dist
      double precision scaleval
      integer i,j,npts
      logical found_crust

      scaleval=dsqrt(PI*GRAV*RHOAV)*R_EARTH
      
      call read_crustal_model  

      write(*,*) 'Number of points : '
      read(*,*) npts
      write(*,*) 'Input start depth for cut (km) : '
      read(*,*) mindepth
      write(*,*) 'Input end depth for cut (km) : '
      read(*,*) maxdepth
      write(*,*) 'Depth step (km) :'
      read(*,*) dstep

      open(unit=11,file='track.xyz',status='old')
      open(unit=10,file='crust.dat',status='unknown') 
  
      do i=1,npts
        read(11,*) lon, lat, dist
        theta=lat
        if (lon > 180.0) then
          phi = lon-360.0
        else
          phi = lon
        endif
        depth=mindepth
        do while (depth .le. maxdepth)
          r=1-depth*1000.0d0/R_EARTH
          call crustal_model(theta,phi,r,vp,vs,rho,moho,found_crust)
          write(10,'6(f8.2,1x)') lon, lat, depth, dist, vs*scaleval/1000.0d0, R_EARTH*(1-moho)/1000.0d0
          depth=depth+dstep
        enddo
      enddo

      close(10)
      close(11)
      

      end program crust_cut
