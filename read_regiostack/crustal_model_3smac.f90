!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 5
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology July 2004
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

!
! read 3SMAC model by Nataf & Ricard
! subroutine written by Alessia Maggi (Feb 2005)
!
module crustal_model_constants

  implicit none
  include 'constants.h'

!  crustal model parameters for 3-SMAC
   integer, parameter :: N_RADII = 11
   integer, parameter :: NLAT = 90
   integer, parameter :: NLON = 180

   integer, parameter :: TOP_INDEX_WITH_SEDIMENTS = 3
   integer, parameter :: TOP_INDEX_NO_SEDIMENTS = 5
   integer, parameter :: MOHO_INDEX = 8

   logical, parameter :: INCLUDE_SEDIMENTS_CRUST = .false.

   double precision, parameter :: XPAS = 2.0d0
   double precision, parameter :: YPAS = 2.0d0

   logical, parameter :: DO_CAP_SMOOTHING = .true.
   logical, parameter :: DO_RADIUS_SHIFTING = .true.
   

end module crustal_model_constants

module crustal_model_variables

  use crustal_model_constants
  implicit none

   double precision smac_vs(NLAT,NLON,N_RADII)
   double precision smac_vp(NLAT,NLON,N_RADII)
   double precision smac_rho(NLAT,NLON,N_RADII)

   double precision smac_radii(NLAT,NLON,N_RADII) 
   double precision smac_moho(NLAT,NLON) 
 
   double precision low_crust_thickness(NLAT,NLON)
   double precision upper_crust_thickness(NLAT,NLON)


   double precision grid_phi(1:NLON)
   double precision grid_theta(1:NLAT)

end module crustal_model_variables

!---------------------------
  subroutine crustal_model(xlat,xlon,r,vp,vs,rho,moho,found_crust)

  use crustal_model_variables
  implicit none

  integer, parameter :: NTHETA = 2
  integer, parameter :: NPHI = 10
  integer, parameter :: N_TEST_POINTS = NTHETA*NPHI
  ! cap size is the angle from due north to the border of the cap (half-angle)
  double precision, parameter :: CAP_SIZE_DEG = 2.0
  double precision, parameter :: CAP_SIZE = CAP_SIZE_DEG * DEGREES_TO_RADIANS

  ! argument variables
  double precision xlat,xlon,r,vp,vs,rho,moho
  logical found_crust

  ! local variables
  integer i,j,k, itheta, iphi
  double precision theta, phi, dtheta, dphi
  double precision cap_area
  double precision r_rot, theta_rot, phi_rot
  double precision test_theta(NTHETA), test_phi(NPHI) 
  double precision test_lat(N_TEST_POINTS), test_lon(N_TEST_POINTS) 
  double precision weight(N_TEST_POINTS), w, total_weight
  double precision sint, sinp, cost, cosp
  double precision rotation_matrix(3,3), x(3), xc(3)
  double precision test_vp, test_vs, test_rho, test_moho
  logical test_found_crust


  ! if smoothing flag is not set this is enough
  if (.not. DO_CAP_SMOOTHING) then
    call get_crustal_model(xlat,xlon,r,vp,vs,rho,moho,found_crust)
    return
  endif

  ! control reaches here only if smoothing flag is set 
    
  ! do some spherical geometry to find the xlat, xlon and integration weights 
  ! for all points in a spherical cap around the requested xlat, xlon
  theta = (90.0-xlat)*PI/180.0
  phi = xlon*PI/180.0

  sint = sin(theta)
  cost = cos(theta)
  sinp = sin(phi)
  cosp = cos(phi)

  ! set up rotation matrix to go from cap at North pole
  ! to cap around point of interest
  rotation_matrix(1,1) = cosp*cost
  rotation_matrix(1,2) = -sinp
  rotation_matrix(1,3) = cosp*sint
  rotation_matrix(2,1) = sinp*cost
  rotation_matrix(2,2) = cosp
  rotation_matrix(2,3) = sinp*sint
  rotation_matrix(3,1) = -sint
  rotation_matrix(3,2) = 0.0
  rotation_matrix(3,3) = cost

  dtheta = CAP_SIZE/dble(NTHETA)
  dphi = 2.0*PI/dble(NPHI)
  cap_area = 2.0*PI*(1.0-cos(CAP_SIZE))

  ! set up the locations of the points in this spherical cap that we will be
  ! using for the integration
  do itheta = 1, NTHETA
    test_theta(itheta) = 0.5*(2.0*itheta-1)*dtheta
  enddo
  do iphi = 1, NPHI
    test_phi(iphi) = (iphi-1)*dphi
  enddo

  ! set up the weights (do here for clarity - for speed would do later)
  k=0
  total_weight=0
  do itheta = 1, NTHETA
    w = sin(test_theta(itheta))*dtheta*dphi/cap_area
    do iphi = 1, NPHI
      k=k+1
      weight(k) = w
      total_weight = total_weight + weight(k)
    enddo
  enddo

  if(abs(total_weight-1.0) > 0.001) stop 'error in cap integration for 3smac'

  ! find the test_lat, test_lon corresponding to to the test_theta test_phi
  k=0
  do itheta = 1, NTHETA
    cost = cos(test_theta(itheta))
    sint = sin(test_theta(itheta))
    do iphi = 1, NPHI
      k=k+1
      cosp = cos(test_phi(iphi))
      sinp = sin(test_phi(iphi))
!     x,y,z coordinates of integration point in cap at North pole
      xc(1) = sint*cosp
      xc(2) = sint*sinp
      xc(3) = cost
!     get x,y,z coordinates in cap around point of interest
      do i=1,3
        x(i) = 0.0
        do j=1,3
          x(i) = x(i)+rotation_matrix(i,j)*xc(j)
        enddo
      enddo
!     get latitude and longitude (degrees) of integration point
      call xyz_2_rthetaphi_dble(x(1),x(2),x(3),r_rot,theta_rot,phi_rot)
      call reduce(theta_rot,phi_rot)
      test_lat(k) = (PI/2.0-theta_rot)*180.0/PI
      test_lon(k) = phi_rot*180.0/PI
      if(test_lon(k) > 180.0) test_lon(k) = test_lon(k)-360.0
    enddo
  enddo

  ! do a weighted integration over the cap
  vp=0
  vs=0
  rho=0
  moho=0
  do k=1,N_TEST_POINTS
    call get_crustal_model(test_lat(k), test_lon(k), r, test_vp, test_vs, test_rho, test_moho, test_found_crust)
    vp=vp+test_vp*weight(k)
    vs=vs+test_vs*weight(k)
    rho=rho+test_rho*weight(k)
    moho=moho+test_moho*weight(k)
  enddo

  found_crust = test_found_crust
!  write(*,*) 'DEBUG2, found_crust ', found_crust

  return

  end subroutine crustal_model


!---------------------------


  subroutine get_crustal_model(xlat,xlon,x,vp,vs,rho,moho,found_crust)

  use crustal_model_variables
  implicit none

  double precision xlat,xlon,x,vp,vs,rho,moho
  logical found_crust

  integer i,k, top_index
  integer jtheta, iphi, krad, kr1, kr2
  integer box_theta(1:8), box_phi(1:8)
  logical pole
  double precision radius, theta, phi, dtheta, dphi,dr
  double precision vs1, vs2, rho1, rho2, vp1, vp2
  double precision interp_radii(1:N_RADII), interp_vs(1:N_RADII)
  double precision interp_vp(1:N_RADII), interp_rho(1:N_RADII)
  double precision box(1:8) 
  double precision scaleval

  vp=0.0d0
  vs=0.0d0
  rho=0.0d0
  moho=0.0d0
  found_crust=.true.

  if (INCLUDE_SEDIMENTS_CRUST) then
    top_index = TOP_INDEX_WITH_SEDIMENTS
  else
    top_index = TOP_INDEX_NO_SEDIMENTS
  endif

  ! radius in SMAC is defined in km, not m
  radius = x*R_EARTH / 1000.0d0
  ! xlat is latitude in degrees, we want colatitude in radians
  theta = (90.0-xlat) * DEGREES_TO_RADIANS
  ! xlon is longitude in degrees, from -180 to 180 degrees 
  if (xlon < 0) then
    phi = (xlon + 360.0 ) * DEGREES_TO_RADIANS
  else
    phi = xlon * DEGREES_TO_RADIANS
  endif

  pole = .false.
  ! bracket the theta phi values
  call bracket_crust(grid_theta,NLAT,theta,jtheta)
  call bracket_crust(grid_phi,NLON,phi,iphi)
  if ((jtheta .eq. 0) .or. (jtheta .eq. NLAT)) pole = .true.

  if (.not. pole) then
    ! we are outside the polar regions, so we need to set up an interpolation box
    box_theta(1) = jtheta
    box_theta(2) = jtheta+1
    box_theta(3) = jtheta+1
    box_theta(4) = jtheta
    dtheta=(theta - grid_theta(jtheta))/(grid_theta(jtheta+1) - grid_theta(jtheta))
    if ((iphi.eq.0) .or. (iphi.eq.NLON)) then
      ! make interpolation box go over the discontinuity
      box_phi(1) = NLON
      box_phi(2) = NLON
      box_phi(3) = 1
      box_phi(4) = 1
      if (iphi.eq.0) then
        dphi=1+(phi - grid_phi(1))/(XPAS*DEGREES_TO_RADIANS)
      else
        dphi=(phi - grid_phi(NLON))/(XPAS*DEGREES_TO_RADIANS)
      endif
    else
      ! standard interpolation box
      box_phi(1) = iphi
      box_phi(2) = iphi
      box_phi(3) = iphi+1
      box_phi(4) = iphi+1
      dphi=(phi - grid_phi(iphi))/(grid_phi(iphi+1)-grid_phi(iphi))
    endif
  endif
  
  ! setup the interpolated arrays for my (theta,phi) point
  do k=1,N_RADII
    interp_radii(k) = 0
  enddo
  if (pole) then
  ! deal with polar regions
    do k=1,N_RADII
      do i=1,NLON
        if (jtheta.eq.0) then    
          interp_radii(k)=interp_radii(k)+smac_radii(1,i,k)
        else
          interp_radii(k)=interp_radii(k)+smac_radii(NLAT,i,k)
        endif
      enddo
      interp_radii(k)=interp_radii(k)/(1.0d0*NLON)
    enddo
  else
  ! non polar regions - use interpolation boxes
    do k=1,N_RADII
      box(1)=smac_radii(box_theta(1),box_phi(1),k)
      box(2)=smac_radii(box_theta(2),box_phi(2),k)
      box(3)=smac_radii(box_theta(3),box_phi(3),k)
      box(4)=smac_radii(box_theta(4),box_phi(4),k)
      interp_radii(k) = (1-dtheta)*(1-dphi)*box(1) + &
                        (dtheta)  *(1-dphi)*box(2) + &
                        (dtheta)  *(dphi)  *box(3) + &
                        (1-dtheta)*(dphi)  *box(4)
    enddo
  endif

  ! now have radii interpolated for my (theta,phi) position, so do some bracketing
  call bracket_crust(interp_radii,N_RADII,radius,krad)

  ! set up short cuts to the two interfaces, kr1=krad, kr2=krad+1
  ! and fix up kr1, kr2 to deal with the cases of krad=0 (above the model)
  ! and krad=N_RADII (below the last interface)
  kr1=krad
  kr2=krad+1
  if (krad .eq. N_RADII) then
    kr1=N_RADII
    kr2=N_RADII
  else if (krad < top_index ) then
    kr1=top_index
    kr2=top_index
  endif
  if( (kr1.ne.kr2) .and. ( interp_radii(kr1).eq.interp_radii(kr2) ) ) then
    ! we have somehow bracketed an infenitesimal space - this should be 
    ! impossible, but in case it does happen set kr1 = kr2
    kr1=kr2
  endif
 
  ! now know which interfaces we are in between, so find the values of vs, vp, rho
  ! on the two interfaces 1=krad, 2=krad+1
  if (pole) then
    vs1=0
    vs2=0
    vp1=0
    vp2=0
    rho1=0
    rho2=0
    moho=0
    do i=1,NLON
      if (jtheta.eq.0) then    
        vs1=vs1+smac_vs(1,i,kr1)
        vs2=vs2+smac_vs(1,i,kr2)
        vp1=vp1+smac_vp(1,i,kr1)
        vp2=vp2+smac_vp(1,i,kr2)
        rho1=rho1+smac_rho(1,i,kr1)
        rho2=rho2+smac_rho(1,i,kr2)
        moho=moho+smac_moho(1,i)
      else
        vs1=vs1+smac_vs(NLAT,i,kr1)
        vs2=vs2+smac_vs(NLAT,i,kr2)
        vp1=vp1+smac_vp(NLAT,i,kr1)
        vp2=vp2+smac_vp(NLAT,i,kr2)
        rho1=rho1+smac_rho(NLAT,i,kr1)
        rho2=rho2+smac_rho(NLAT,i,kr2)
        moho=moho+smac_moho(NLAT,i)
      endif
    enddo
    vs1=vs1/(1.0d0*NLON)
    vs2=vs2/(1.0d0*NLON)
    vp1=vp1/(1.0d0*NLON)
    vp2=vp2/(1.0d0*NLON)
    rho1=rho1/(1.0d0*NLON)
    rho2=rho2/(1.0d0*NLON)
    moho=moho/(1.0d0*NLON)
  else !(not on the pole)
    !vs1
    box(1)=smac_vs(box_theta(1),box_phi(1),kr1)
    box(2)=smac_vs(box_theta(2),box_phi(2),kr1)
    box(3)=smac_vs(box_theta(3),box_phi(3),kr1)
    box(4)=smac_vs(box_theta(4),box_phi(4),kr1)
    vs1 = (1-dtheta)*(1-dphi)*box(1) + &
          (dtheta)  *(1-dphi)*box(2) + &
          (dtheta)  *(dphi)  *box(3) + &
          (1-dtheta)*(dphi)  *box(4)
    !vs2
    box(1)=smac_vs(box_theta(1),box_phi(1),kr2)
    box(2)=smac_vs(box_theta(2),box_phi(2),kr2)
    box(3)=smac_vs(box_theta(3),box_phi(3),kr2)
    box(4)=smac_vs(box_theta(4),box_phi(4),kr2)
    vs2 = (1-dtheta)*(1-dphi)*box(1) + &
          (dtheta)  *(1-dphi)*box(2) + &
          (dtheta)  *(dphi)  *box(3) + &
          (1-dtheta)*(dphi)  *box(4)
    !vp1
    box(1)=smac_vp(box_theta(1),box_phi(1),kr1)
    box(2)=smac_vp(box_theta(2),box_phi(2),kr1)
    box(3)=smac_vp(box_theta(3),box_phi(3),kr1)
    box(4)=smac_vp(box_theta(4),box_phi(4),kr1)
    vp1 = (1-dtheta)*(1-dphi)*box(1) + &
          (dtheta)  *(1-dphi)*box(2) + &
          (dtheta)  *(dphi)  *box(3) + &
          (1-dtheta)*(dphi)  *box(4)
    !vp2
    box(1)=smac_vp(box_theta(1),box_phi(1),kr2)
    box(2)=smac_vp(box_theta(2),box_phi(2),kr2)
    box(3)=smac_vp(box_theta(3),box_phi(3),kr2)
    box(4)=smac_vp(box_theta(4),box_phi(4),kr2)
    vp2 = (1-dtheta)*(1-dphi)*box(1) + &
          (dtheta)  *(1-dphi)*box(2) + &
          (dtheta)  *(dphi)  *box(3) + &
          (1-dtheta)*(dphi)  *box(4)
    !rho1
    box(1)=smac_rho(box_theta(1),box_phi(1),kr1)
    box(2)=smac_rho(box_theta(2),box_phi(2),kr1)
    box(3)=smac_rho(box_theta(3),box_phi(3),kr1)
    box(4)=smac_rho(box_theta(4),box_phi(4),kr1)
    rho1 = (1-dtheta)*(1-dphi)*box(1) + &
          (dtheta)  *(1-dphi)*box(2) + &
          (dtheta)  *(dphi)  *box(3) + &
          (1-dtheta)*(dphi)  *box(4)
    !rho2
    box(1)=smac_rho(box_theta(1),box_phi(1),kr2)
    box(2)=smac_rho(box_theta(2),box_phi(2),kr2)
    box(3)=smac_rho(box_theta(3),box_phi(3),kr2)
    box(4)=smac_rho(box_theta(4),box_phi(4),kr2)
    rho2 = (1-dtheta)*(1-dphi)*box(1) + &
          (dtheta)  *(1-dphi)*box(2) + &
          (dtheta)  *(dphi)  *box(3) + &
          (1-dtheta)*(dphi)  *box(4)
    !moho
    box(1)=smac_moho(box_theta(1),box_phi(1))
    box(2)=smac_moho(box_theta(2),box_phi(2))
    box(3)=smac_moho(box_theta(3),box_phi(3))
    box(4)=smac_moho(box_theta(4),box_phi(4))
    moho = (1-dtheta)*(1-dphi)*box(1) + &
          (dtheta)  *(1-dphi)*box(2) + &
          (dtheta)  *(dphi)  *box(3) + &
          (1-dtheta)*(dphi)  *box(4)
  endif

  ! if below the moho then found_crust is false
  ! keep going to retrieve subcrustal vs, vp, rho
  
  if (radius < moho) then
    found_crust=.false.
  endif
!  write(*,*) 'DEBUG CRUST rad, moho', radius, moho, found_crust



  ! now have vs1,vs2,vp1,vp2,rho1,rho2 and the corresponding interface radii 
  ! interp_radii(kr1) and interp_radii(kr2) so can lineraly interpolate
  if (kr1.eq.kr2) then
    vs=vs1
    vp=vp1
    rho=rho1
  else
    dr=(radius-interp_radii(kr1)) / (interp_radii(kr2)-interp_radii(kr1))
    vs= (1-dr)*vs1  + dr*vs2
    vp= (1-dr)*vp1  + dr*vp2
    rho=(1-dr)*rho1 + dr*rho2
  endif

! non-dimensionalize
!  scaleval = dsqrt(PI*GRAV*RHOAV)
!  vp = vp*1000.0d0/(R_EARTH*scaleval)
!  vs = vs*1000.0d0/(R_EARTH*scaleval)
!  rho = rho*1000.0d0/RHOAV

! keep dimensions
  vp = vp
  vs = vs
  rho = rho

  moho = (R_EARTH-moho*1000.0d0)/R_EARTH

  return
  end subroutine get_crustal_model

!---------------------------

  subroutine read_crustal_model

  use crustal_model_variables
  implicit none
  character*80 radname,vpname,vsname,rhoname

! local variables
  integer i,j,k,m,n

  double precision get_theta_crust, get_phi_crust
  double precision r_corr

! set up the theta-phi coordinates for later searching
  do i=1,NLAT
    grid_theta(i)=get_theta_crust(i)
  enddo
  do i=1,NLON
    grid_phi(i)=get_phi_crust(i)
  enddo

! read the 3SMAC files
! start by reading the crustal thickness file
  open(unit=10,file="DATA/3-SMAC/brut/datacroute",status='old')
  read(10,*)
  read(10,*)
  read(10,*)
  do j=1,NLON
    ! crust thick in m
    read(10,*)(smac_moho(k,j),k=1,NLAT) 
    ! in km to be consistent with radii
  enddo
  close(10)
  do j=1,NLON
    do k=1,NLAT
      smac_moho(k,j)=(R_EARTH-smac_moho(k,j))/1000.0 
    enddo
  enddo
  
! now read the radii and parameter files
  do i=1,N_RADII
    ! setup filenames for each index
    if (i < 10 ) then
      write(radname,'("DATA/3-SMAC/radi/radi.",i1)') i
      write(vpname,'("DATA/3-SMAC/para/VP.",i1)') i
      write(vsname,'("DATA/3-SMAC/para/VS.",i1)') i
      write(rhoname,'("DATA/3-SMAC/para/RHO.",i1)') i
    else
      write(radname,'("DATA/3-SMAC/radi/radi.",i2)') i
      write(vpname,'("DATA/3-SMAC/para/VP.",i2)') i
      write(vsname,'("DATA/3-SMAC/para/VS.",i2)') i
      write(rhoname,'("DATA/3-SMAC/para/RHO.",i2)') i
    endif

    ! open all files
    open(unit=10,file=radname,status='old')
    open(unit=11,file=vpname,status='old')
    open(unit=12,file=vsname,status='old')
    open(unit=13,file=rhoname,status='old')

    ! read 3 blank lines from the top of each file
    do m=10,13
      read(m,*)
      read(m,*)
      read(m,*)
    enddo

    ! read the information from each file
    do j=1,NLON
      read(10,*)(smac_radii(k,j,i),k=1,NLAT)
      read(11,*)(smac_vp(k,j,i),k=1,NLAT)
      read(12,*)(smac_vs(k,j,i),k=1,NLAT)
      read(13,*)(smac_rho(k,j,i),k=1,NLAT)
    enddo

!   close all files 
    do m=10,13
      close(m)
    enddo

  enddo

  if(DO_RADIUS_SHIFTING) then
! correct the radii so that smac_radii(,,MOHO_INDEX) lies on the moho defined 
! by smac_moho(,)
    do j=1,NLON
      do k=1,NLAT
        r_corr=smac_moho(k,j)-smac_radii(k,j,MOHO_INDEX)
        do i=1,N_RADII
          smac_radii(k,j,i)=smac_radii(k,j,i)+r_corr
        enddo
      enddo
    enddo
  endif

  end subroutine read_crustal_model

!---------------------------

  double precision function get_theta_crust(iy)

  use crustal_model_variables
  implicit none

  integer iy

  get_theta_crust = ((iy-1) * YPAS+YPAS/2.0d0) * DEGREES_TO_RADIANS

  end function get_theta_crust

!----------------------------------

  double precision function get_phi_crust(ix)

  use crustal_model_variables
  implicit none

  integer ix

  get_phi_crust = ( (ix-1) * XPAS+XPAS/2.0d0) * DEGREES_TO_RADIANS

  end function get_phi_crust

!----------------------------------

  subroutine bracket_crust(xx,n,x,j)
  integer j,n
  double precision x, xx(1:n)
! Given an array xx(1:n) and a value x, returns j such that x is between
! xx(j) and xx(j+1).  xx(1:n) must be monotonic.  j=0 or j=n+1 are returned
! to indecate that x is out of range.  Modelled on numercal recipies locate.

  integer jl,jm,ju

  jl=0
  ju=n+1

! find values by bisection
  do while (ju-jl > 1) 
    jm=(ju+jl)/2
    if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm))) then
      jl=jm
    else
      ju=jm
    endif
  enddo

! setup return values
  if(x.eq.xx(1)) then
    j=1
  else if (x.eq.xx(n)) then
    j=n-1
  else
    j=jl
  endif

  return

  end subroutine bracket_crust


!----------------------------------

