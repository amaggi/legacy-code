C---------------------------------------------------------------------
C $Id:$
C---------------------------------------------------------------------
C              CHECKERBOARD SPECIFIC FUNCTIONS
C---------------------------------------------------------------------

c---------------------------------------------------------------------
c checkersign: 
c   Given lat, lon and perturbation size in km, 
c   finds the closest spherical tesseral function on the sphere
c   and returns +/-1, the sign
c   of an alternating checkerboard pattern at (lat, lon).  Checks are
c   `squares' of equal area.  Pattern is symmetric about 0N and 210E
c   (centered for the pacific ocean).
c---------------------------------------------------------------------

      function checkersign(lat, lon, pert_size)

      real lat, lon, pert_size, checkersign
      real pi, deg2rad, radius
      real undeggc, lambdagc, nlambdagc, psign

      pi=3.1415926
      deg2rad=pi/180.
      radius=6371.

      undeggc=(2*pi*radius)/360. 
      lambdagc=2*pert_size/undeggc
      nlambdagc=int(360./lambdagc)

      psign=cos(nlambdagc*lat*deg2rad)*cos(nlambdagc*lon*deg2rad)
      if (psign.gt.0) psign=1
      if (psign.lt.0) psign=-1
      checkersign=psign

      end function

c---------------------------------------------------------------------
c  vschecker:
c    returns vs perturbed according to checkerboard pattern.  Calls
c    checkersign to determine sign of perturbation.
c---------------------------------------------------------------------

      function vschecker(c,lat,lon,prof,pert_size,pmin,pmax,mag)
      real c,lat,lon, pert_size, mag, pmin, pmax, prof
      integer psign;

      if (prof.le.pmax.and.prof.ge.pmin) then
        psign=checkersign(lat,lon,pert_size)
        vschecker=c+ psign*mag*c
      else
        vschecker=c
      endif
      end function

c---------------------------------------------------------------------
c  a2checker:
c    returns a2 perturbed according to checkerboard pattern. a2 varies
c    s.t. anisotropy either points 45N or 135N.  Calls checkersign to
c    obtain sign of perturbation.
c---------------------------------------------------------------------

      function a2checker(lat,lon,prof,pert_size,pmin,pmax)
      real lat,lon, pert_size, pmin, pmax, prof
      integer psign;

      if (prof.le.pmax.and.prof.ge.pmin) then
        psign=checkersign(lat,lon,pert_size)
        a2checker=psign*0.045
      else
        a2checker=0.0
      endif

      end function

c---------------------------------------------------------------------
c  vs_cyl:
c    returns vs perturbed according to a cylinder with parameters given
c    in the function call
c---------------------------------------------------------------------

      function vs_cyl(c,lat,lon,prof,cyl_lat,cyl_lon,cyl_rad,
     *                cyl_maxprof,cyl_minprof,cyl_mag,n_cyl)
      real c, lat, lon, prof, cyl_lat, cyl_lon, cyl_rad
      real cyl_maxprof, cyl_minprof, cyl_mag
      integer n_cyl, i
      
      dimension cyl_lat(*), cyl_lon(*), cyl_rad(*)
      dimension cyl_maxprof(*), cyl_minprof(*), cyl_mag(*)

      real az, baz, dkm, ddg
      real vs

      vs=c

      do 10, i=1, n_cyl
      if (prof.le.cyl_maxprof(i).and.prof.ge.cyl_minprof(i)) then
        call distaz(lat, lon, cyl_lat(i), cyl_lon(i), az, baz, ddg, dkm)
        if(dkm .le. cyl_rad(i)) then 
          vs=vs+cyl_mag(i)*vs
        endif
      endif
10    continue

      vs_cyl = vs
      end function
c---------------------------------------------------------------------
c  vs_plume:
c    returns vs perturbed according to a plume with parameters given
c    in the function call
c---------------------------------------------------------------------

      function vs_plume(c,lat,lon,prof,plume_lat,plume_lon,plume_rad,
     *                plume_maxprof,plume_minprof,plume_DT,n_plume,
     *                plume_grad)
      real c, lat, lon, prof, plume_lat, plume_lon, plume_rad
      real plume_maxprof, plume_minprof, plume_DT, plume_grad
      integer n_plume, i
      
      dimension plume_lat(*), plume_lon(*), plume_rad(*)
      dimension plume_maxprof(*), plume_minprof(*), plume_DT(*)

      real az, baz, dkm, ddg
      real vs

      vs=c

      do 10, i=1, n_plume
      if (prof.le.plume_maxprof(i).and.prof.ge.plume_minprof(i)) then
        call distaz(lat,lon,plume_lat(i),plume_lon(i),az,baz,ddg,dkm)
        if(dkm .le. plume_rad(i)) then 
          vs=vs+(plume_DT(i)/100.0)*plume_grad*vs/100.0
        endif
      endif
10    continue

      vs_plume = vs

      end function




