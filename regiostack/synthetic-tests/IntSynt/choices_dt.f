c---------------------------------------------------------
c   set_uniform:
c     rets a uniform value to the c array
c---------------------------------------------------------
      SUBROUTINE set_uniform(c,cval)
      include 'param_V3.h'
      dimension c(NPROFMAX,NLIN,NCOL)
      real cval

      do i=1,nlon
        do j=1,nlat
           do k=1,ncouch
             c(k,j,i)=cval
           enddo
        enddo
      enddo

      end

c---------------------------------------------------------------------
c  dtchecker:
c    returns dt perturbed according to checkerboard pattern.  Calls
c    checkersign to determine sign of perturbation.
c---------------------------------------------------------------------

      function dtchecker(dt,lat,lon,prof,pert_size,pmin,pmax,mag)
      real dt,lat,lon, pert_size, mag, pmin, pmax, prof
      integer psign;

      if (prof.le.pmax.and.prof.ge.pmin) then
        psign=checkersign(lat,lon,pert_size)
        dtchecker=dt+ psign*mag
      else
        dtchecker=dt
      endif
      end function


c---------------------------------------------------------------------
c  dt_cyl:
c    returns dt perturbed according to a cylinder with parameters given
c    in the function call
c---------------------------------------------------------------------

      function dt_cyl(dt,lat,lon,prof,cyl_lat,cyl_lon,cyl_rad,
     *                cyl_maxprof,cyl_minprof,cyl_mag,n_cyl)
      real dt, lat, lon, prof, cyl_lat, cyl_lon, cyl_rad
      real cyl_maxprof, cyl_minprof, cyl_mag
      integer n_cyl, i
      
      dimension cyl_lat(*), cyl_lon(*), cyl_rad(*)
      dimension cyl_maxprof(*), cyl_minprof(*), cyl_mag(*)

      real az, baz, dkm, ddg
      real dt_tmp 

      dt_tmp=dt

      do 10, i=1, n_cyl
      if (prof.le.cyl_maxprof(i).and.prof.ge.cyl_minprof(i)) then
        call distaz(lat, lon, cyl_lat(i), cyl_lon(i), az, baz, ddg, dkm)
        if(dkm .le. cyl_rad(i)) then 
          dt_tmp=dt_tmp+cyl_mag(i)
        endif
      endif
10    continue

      dt_cyl = dt_tmp
      end function

