c---------------------------------------------------------------------
c $Id:$
c---------------------------------------------------------------------
c	LITHOSPHERIC COOLING FUNCTIONS
c---------------------------------------------------------------------


      function halfspcool_vs(age,c,prof,dlogVsdT,Tm,kappa)
c     input parameters
      real age, c, prof, dlogVsdT, Tm, kappa
     
      real T, vs_pert

c     calculate temperature
      T=Tm*erf(prof/sqrt(4*kappa*age))

c     vs perturbation is calculated from delta T
c     note that dlogVsdT is in % / 100K
      vs_pert=dlogVsdT*(T-Tm)/100.0

c     calculate perturbed velocity
      halfspcool_vs=c*(1+vs_pert/100.0)
      end function halfspcool_vs

c---------------------------------------------------------------------

      function platecool_vs(age,c,prof,dlogVsdT,Tm,kappa,plate_th)
c     input parameters
      real age, c, prof, dlogVsdT, Tm, kappa, plate_th
     
      real T, vs_pert

c     calculate temperature
      T=plate_T(prof,age,plate_th,kappa)

c     vs perturbation is calculated from delta T
c     note that dlogVsdT is in % / 100K
      vs_pert=dlogVsdT*(T-Tm)/100.0

c     calculate perturbed velocity
      platecool_vs=c*(1+vs_pert/100.0)
      end function platecool_vs

c---------------------------------------------------------------------

      function plate_T(z,t,a,kappa)
      real z,t,a,kappa
      integer n, n_max
      real t_min, partial_t
      real pi

c     approximate formula used here is valid for ages > 10 Ma
c     for younger ages, the temperature converges to the value at 10Ma
      pi=3.1415926
      n_max=6
      t_min = 10.0

      if (z.gt.a) then
        plate_T = Tm
        return
      endif
 
      partial_t=0
      if (t .ge. t_min) then
        do n=1,n_max
          partial_t = partial_t + 2/(n*pi) * 
     &                exp (-1*(n*n*pi*pi*kappa*t)/(a*a)) *
     &                sin ((n*pi*z)/a)
        enddo
      else 
        do n=1,n_max
          partial_t = partial_t + 2/(n*pi) * 
     &                exp (-1*(n*n*pi*pi*kappa*t_min)/(a*a)) *
     &                sin ((n*pi*z)/a)
        enddo
      endif
      plate_T = Tm*(z/a + partial_t)

      end function plate_T


c---------------------------------------------------------------------

      subroutine read_age(age,gridstep)

      include 'param_V3.h'
      parameter(NLIN_TEMP=1800,NCOL_TEMP=3600)

      real age, gridstep
      real lat, lon, temp_gridstep, bid
      integer ilat,jlon, nlat, nlon
      integer ilat_temp, jlon_temp
      real temp_age(NLIN_TEMP,NCOL_TEMP)

      dimension age(NLIN,NCOL)

c     read the original Muller et al file
      temp_gridstep=0.1
      open(unit=11,file="/home/alessia/share/ocean-age/age_1.6.xyz")
c     latitude goes from +90 to -72 in age_1.6.xyz file
      do ilat=1,1621 
        do jlon=1,3600
          read(11,*) lat, lon, temp_age(ilat,jlon)
        enddo
          read(11,*) lat, lon, bid
      enddo
      close(11)

c     find the points corresponding to our 1 degree grid
      nlat=int(180/gridstep)
      nlon=int(360/gridstep)
      
      do ilat=1,nlat
        lat=getlat(ilat,gridstep)
        ilat_temp=getilin(lat,temp_gridstep)
        if (ilat_temp.gt.1621) then
          age(ilat,jlon) = 999
        else 
          do jlon=1,nlon
           lon=getlon(jlon,gridstep)
           jlon_temp=geticol(lon,temp_gridstep)
           age(ilat,jlon)=temp_age(ilat_temp,jlon_temp)
          enddo
        endif
      enddo

        
      end subroutine read_age

