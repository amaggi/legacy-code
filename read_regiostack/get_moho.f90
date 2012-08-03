! Test program for reading regiostack model

      program get_moho

      implicit none
      include 'constants.h'

      double precision lat, lon, mylon
      double precision vs, rho, vp, moho
      double precision depth, r
      integer i,j
      logical found_crust
      
      call read_crustal_model  

      depth=2.0
      r=1-depth*1000.0d0/R_EARTH

      open(unit=10,file='moho.dat',status='unknown') 
  
      do j=1,720
        lon=j/2.0d0
        if (lon .gt. 180.0) then
          mylon=lon-360.0
        else 
          mylon=lon
        endif
        do i=1,360
          lat=90-(i/2.0d0)
          call crustal_model(lat,mylon,r,vp,vs,rho,moho,found_crust)
          write(10,'f8.2,1x,f8.2,1x,f8.2') lat, lon, R_EARTH*(1-moho)/1000.0d0
        enddo
      enddo
      

      end program get_moho
