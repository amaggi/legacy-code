

!     function to calculate the anistropy angle for a given point
      function kr_anis_angle(lat,lon,scale_rad,cyl_lat,cyl_lon,n_cyl,
     *                       plate_angle,theta_rot)
   
      parameter (epsilon=1e-5, c=4.0, d=7.4)
      integer n_cyl
      real lat,lon,cyl_lat,cyl_lon, scale_rad 
      real x,y,x_plume,y_plume,x_origin,y_origin,x_prime,y_prime
      
      dimension cyl_lat(*), cyl_lon(*)
      
!     scale coordinate system:
      x=lon/scale_rad
      y=lat/scale_rad


      do k = 1, n_cyl
      
        anis_angle = plate_angle
      
        x_plume = cyl_lon(k)/scale_rad
        y_plume = cyl_lat(k)/scale_rad
        
        x_origin = x_plume + c*sin(theta_rot)
        y_origin = y_plume + c*sin(theta_rot)
        
  	    call rotate(x_prime,y_prime, -1*theta_rot)
	   
	    ! if rotated point is between the two parablas then set the 
	    ! anisotropy angle to be
	    ! such that the fast direction points to the center of the plume
	    if (abs(x_prime) > epsilon .and. y_prime < d .and. 
     *    y_prime > ya(x_prime) .and. y_prime < yb(x_prime) ) 
     *    then
		  anis_angle = atan2(y - y_plume , x - x_plume)
	    end if 
 	  enddo

      end function kr_anis_angle


c     function describing the external parabola
      real function ya(x)
      parameter (a2=0.19)
      real x
	  ya = a2*x*x
      end function ya
      
c     function describing the internal parabola
      real function yb(x)
      parameter (b1=1.74, b2=0.35)
      real x
	  yb = b1+a2*x*x
      end function yb 
      

      subroutine rotate(x,y,th)

      real x,y 
      real th
  
      x = x*cos(th) - y*sin(th)
      y = x*sin(th) + y*cos(th)
  
      end subroutine  