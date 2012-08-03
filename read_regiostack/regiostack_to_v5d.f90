! Test program for reading regiostack model and turning it into a v5d grid

      program depth_cut

      implicit none
      include 'constants.h'

      double precision depth_max, depth_min, depth_step 
      double precision depth, r, theta, phi
      double precision x,y,z
      double precision dvs, drho, dvp
      integer i,j,k,kcell(8),ik
      integer nplane,nrow
      
      call read_mantle_model  

      write(*,*) 'Input min depth (km) : '
      read(*,*) depth_min
      write(*,*) 'Input max depth (km) : '
      read(*,*) depth_max
      write(*,*) 'Input depth step (km) : '
      read(*,*) depth_step



      write (*,*) 'Generating points ...'

      open(unit=10,file='regiostack.xyz',status='unknown') 
  
      depth=depth_min
      do while (depth .le. depth_max)
        r=1-depth*1000.0d0/R_EARTH
        do j=1,180
          theta=j*DEGREES_TO_RADIANS*1.0d0
          z=r*cos(theta)
          do i=1,360
            phi=i*DEGREES_TO_RADIANS*1.0d0
            x=r*sin(theta)*cos(phi)
            y=r*sin(theta)*sin(phi)
            call mantle_model(r,theta,phi,dvs,dvp,drho)
            write(10,'(f8.6,1x,f8.6,1x,f8.6,1x,f8.3)') x,y,z,dvs*100
          enddo
        enddo
        depth=depth+depth_step
      enddo

      close(10)
      
      write (*,*) 'Generating cells ...'

      open(unit=10,file='regiostack.cell',status='unknown') 

      nplane=180*360
      nrow=360
      k=1
      depth=depth_min
      do while (depth .lt. depth_max)
        do j=1,179
          do i=1,360
            kcell(1) = (k-1)*nplane+(j-1)*nrow+i
            if (i.eq.360) then
              kcell(2) = (k-1)*nplane + (j-1)*nrow + 1
            else
              kcell(2) = kcell(1)+1
            endif
            kcell(3) = kcell(2)+nrow
            kcell(4) = kcell(1)+nrow
            kcell(5) = kcell(1)+nplane
            kcell(6) = kcell(2)+nplane
            kcell(7) = kcell(3)+nplane
            kcell(8) = kcell(4)+nplane
            write(10,'(8(i5,1x))') (kcell(ik), ik=1,8)
          enddo
        enddo
        depth=depth+depth_step
        k=k+1
      enddo

      close(10)

      end program 
