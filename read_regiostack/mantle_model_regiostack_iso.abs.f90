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
!  Alessia Maggi - February 2005
!
!  subroutines for inserting regiostack isotropic vsfin models 
!  into specfem3D 
!
!=====================================================================

module three_d_mantle_model_constants

  implicit none
  include 'constants.h'

  ! set RMOHO to PREM's average moho depth
 ! double precision, parameter ::    RMOHO = 6346600.d0
  ! set RMOHO to 40km depth
  double precision, parameter ::    RMOHO = 6331000.d0
  ! set other PREM interfaces
  double precision, parameter ::    R80 = 6291000.d0
  double precision, parameter ::    R220 = 6151000.d0
  double precision, parameter ::    R400 = 5971000.d0
  double precision, parameter ::    R600 = 5771000.d0
  double precision, parameter ::    R670 = 5701000.d0
  double precision, parameter ::    R771 = 5600000.d0
  double precision, parameter ::    RTOPDDOUBLEPRIME = 3630000.d0
  double precision, parameter ::    RCMB = 3480000.d0
  double precision, parameter ::    RICB = 1221000.d0


  ! regiostack model at 1x1 degrees 
  integer, parameter :: NLAT = 180, NLON = 360, NLAY = 29
  double precision, parameter :: XPAS = 1.0d0
  double precision, parameter :: YPAS = 1.0d0
  double precision, parameter :: X0 = -180.0d0
  double precision, parameter :: Y0 = 0.0d0

  ! background model problem
  logical, parameter :: WRT_SMOOTHED_PREM = .true.
  logical, parameter :: OUTPUT_ABSOLUTE = .true.


end module three_d_mantle_model_constants

module three_d_mantle_model_variables

  use three_d_mantle_model_constants
  implicit none

  double precision grid_vs(1:NLAT,1:NLON,1:NLAY)
  double precision grid_r(1:NLAY)
  double precision grid_phi(1:NLON)
  double precision grid_theta(1:NLAT)

end module three_d_mantle_model_variables

!-----------------

  subroutine read_mantle_model

  use three_d_mantle_model_variables
  implicit none

  integer i,j,k
  double precision depth, rho
  double precision get_theta, get_phi
  double precision dr, bg_vs, smooth_prem(1:NLAY)

  double precision quick_prem_vsv

! read the smoothed prem data in case we want to use it
  open(unit=10,file='DATA/regiostack/des.prem.rovs.int',status='old')
  read(10,*)
  read(10,*)
  read(10,*)
  do k=1,NLAY
   read(10,*) depth, rho, smooth_prem(k)
  enddo
  close(10)

! read vsfin file
  open(unit=10,file='vsfin',status='old')
  do k=1,NLAY
    read(10,*) depth
    grid_r(k) = 1.0d0 - (depth*1000.d0)/R_EARTH
    if(WRT_SMOOTHED_PREM) then
      bg_vs = smooth_prem(k)
    else
      bg_vs = quick_prem_vsv(grid_r(k))
    endif
    ! read the absolute velocities
    do j=1,NLAT
      read(10,*) (grid_vs(j,i,k),i=1,NLON)  
      do i=1,NLON
        ! Turn absolute velocites into relative velocites
          grid_vs(j,i,k) = (grid_vs(j,i,k) / bg_vs) - 1.0d0
        ! If we want to output absolute velocities, then apply 
        ! the perturbation to PREM
        if (OUTPUT_ABSOLUTE) then
           grid_vs(j,i,k) = (grid_vs(j,i,k) + 1.0d0) * bg_vs
        endif
      enddo
    enddo
  enddo
  close(10)

! set up the theta / phi arrays as read from vsfin file
! phi : -PI -> PI
  do i=1,NLAT
    grid_theta(i)=get_theta(i)
  enddo
  do i=1,NLON
    grid_phi(i)=get_phi(i)
  enddo

  end subroutine read_mantle_model

!---------------------------


  subroutine mantle_model(radius,theta,phi,dvs,dvp,drho)

  use three_d_mantle_model_variables
  implicit none

! factor to convert perturbations in shear speed to perturbations in density
  double precision, parameter :: SCALE_RHO = 0.40d0

  double precision radius,theta,phi,dvs,dvp,drho

  integer i,j,k, kr1, kr2
  integer iphi,jtheta,krad
  double precision dvs1,dvs2
  double precision myphi,dr,dtheta,dphi
  double precision r_moho,r_model_bottom
  double precision vsbox(1:8)
  integer box_theta(1:4), box_phi(1:4)
  logical pole

  dvs = ZERO
  dvp = ZERO
  drho = ZERO
  pole = .false.

  r_moho = RMOHO / R_EARTH
  r_model_bottom = grid_r(NLAY) / R_EARTH

  ! if we are below the range of my model, default to d_stuff=0
  !if( (radius > r_moho) .or. (radius < r_model_bottom)) return
  if( (radius < r_model_bottom)) return

  ! find bracketing indexes for theta, phi, radius
  call bracket(grid_r,NLAY,radius,krad)
  call bracket(grid_theta,NLAT,theta,jtheta)
  if ((jtheta.eq.0) .or. (jtheta.eq.NLAT)) pole=.true.

  ! sort out the different periodicity problem
  ! do this on a copy myphi of original phi in case calling routine 
  ! wants to reuse the phi value
  if (phi>PI) then
    myphi=phi-TWO_PI
  else
    myphi=phi
  endif
  call bracket(grid_phi,NLON,myphi,iphi)
  
  ! recheck for radius limits
  if ((krad.eq.NLAY)) then
    return
  elseif (krad.eq.0) then
    kr1=1
    kr2=1
    dr=0
  else
    kr1=krad
    kr2=krad+1
    dr=(radius - grid_r(krad))/(grid_r(kr2) - grid_r(kr1))
  endif

  ! check for latitude limits
  if (pole) then 
    dvs1=0
    dvs2=0
    ! pole region : return average of points surrounding polar cap
    ! interpolating for depth
    do i=1,NLON
      if (jtheta.eq.0) then
        dvs1=dvs1+grid_vs(1,i,kr1)
        dvs2=dvs2+grid_vs(1,i,kr2)
      else
        dvs1=dvs1+grid_vs(1,i,kr1)
        dvs2=dvs2+grid_vs(1,i,kr2)
      endif
    enddo
    dvs1=dvs1/(1.0d0*NLON)
    dvs2=dvs2/(1.0d0*NLON)
    dvs=(1-dr)*dvs1 + dr*dvs2
    drho=dvs*SCALE_RHO
    return
  else
    box_theta(1)=jtheta
    box_theta(2)=jtheta+1
    box_theta(3)=jtheta+1
    box_theta(4)=jtheta
    dtheta=(theta - grid_theta(jtheta))/(grid_theta(jtheta+1) - grid_theta(jtheta))
  endif

  if ((iphi.eq.0) .or. (iphi.eq.NLON)) then
    box_phi(1)=NLON
    box_phi(2)=NLON
    box_phi(3)=1
    box_phi(4)=1
    if (iphi.eq.0) then
      dphi=1+(myphi - grid_phi(1))/(XPAS*DEGREES_TO_RADIANS)
    else
      dphi=(myphi - grid_phi(NLON))/(XPAS*DEGREES_TO_RADIANS)
    endif
  else
    box_phi(1)=iphi
    box_phi(2)=iphi
    box_phi(3)=iphi+1
    box_phi(4)=iphi+1
    dphi=(myphi - grid_phi(iphi))/(grid_phi(iphi+1) - grid_phi(iphi))
  endif

   ! make interpolation cube 
    vsbox(1)=grid_vs(box_theta(1),  box_phi(1),  kr1)
    vsbox(2)=grid_vs(box_theta(2),  box_phi(2),  kr1)
    vsbox(3)=grid_vs(box_theta(3),  box_phi(3),  kr1)
    vsbox(4)=grid_vs(box_theta(4),  box_phi(4),  kr1)
    vsbox(5)=grid_vs(box_theta(1),  box_phi(1),  kr2)
    vsbox(6)=grid_vs(box_theta(2),  box_phi(2),  kr2)
    vsbox(7)=grid_vs(box_theta(3),  box_phi(3),  kr2)
    vsbox(8)=grid_vs(box_theta(4),  box_phi(4),  kr2)

  ! interpolate over the 3D cube
 
  dvs=(1-dtheta)*(1-dphi)*(1-dr)*vsbox(1) + &
      (dtheta)  *(1-dphi)*(1-dr)*vsbox(2) + &
      (dtheta)  *(dphi)  *(1-dr)*vsbox(3) + &
      (1-dtheta)*(dphi)  *(1-dr)*vsbox(4) + &
      (1-dtheta)*(1-dphi)*  (dr)*vsbox(5) + &
      (dtheta)  *(1-dphi)*  (dr)*vsbox(6) + &
      (dtheta)  *(dphi)  *  (dr)*vsbox(7) + &
      (1-dtheta)*(dphi)  *  (dr)*vsbox(8) 

  drho=SCALE_RHO*dvs
    
  return

  end subroutine mantle_model

!----------------------------------

  subroutine bracket(xx,n,x,j)
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

  end subroutine bracket

!----------------------------------

  double precision function get_theta(iy)
! returns theta for vsfin coordinate iy (uses Y0)

!  include 'constants.h'
  use three_d_mantle_model_variables
  implicit none

  integer iy

  get_theta = (Y0 + (iy-1) * YPAS+YPAS/2.0d0) * DEGREES_TO_RADIANS

  end function get_theta


!----------------------------------

  double precision function get_phi(ix)
! returns phi for vsfin coordinate ix (uses X0)

!  include 'constants.h'
  use three_d_mantle_model_variables
  implicit none

  integer ix

  get_phi = (X0 + (ix-1) * XPAS+XPAS/2.0d0) * DEGREES_TO_RADIANS

  end function get_phi


!----------------------------------

!----------------------------------
!----------------------------------
  double precision function quick_prem_vsv(x)

! cut and pasted from prem_model.f90
  use three_d_mantle_model_variables
  implicit none

  double precision x

  double precision r
  double precision rho, vpv, vsv, vph, vsh, Qmu, Qkappa, eta_aniso


  r=x*R_EARTH

  if(r <= RTOPDDOUBLEPRIME ) then
    stop 'not a mantle depth in quick_prem_vsv'
!
!--- mantle: below d670
!
  else if(r > RTOPDDOUBLEPRIME .and. r <= R771) then
    rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
    vpv=24.9520d0-40.4673d0*x+51.4832d0*x*x-26.6419d0*x*x*x
    vsv=11.1671d0-13.7818d0*x+17.4575d0*x*x-9.2777d0*x*x*x
    vph=vpv
    vsh=vsv
    Qmu=312.0d0
    Qkappa=57827.0d0
  else if(r > R771 .and. r <= R670) then
    rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
    vpv=29.2766d0-23.6027d0*x+5.5242d0*x*x-2.5514d0*x*x*x
    vsv=22.3459d0-17.2473d0*x-2.0834d0*x*x+0.9783d0*x*x*x
    vph=vpv
    vsh=vsv
    Qmu=312.0d0
    Qkappa=57827.0d0
!
!--- mantle: above d670
!
  else if(r > R670 .and. r <= R600) then
    rho=5.3197d0-1.4836d0*x
    vpv=19.0957d0-9.8672d0*x
    vsv=9.9839d0-4.9324d0*x
    vph=vpv
    vsh=vsv
    Qmu=143.0d0
    Qkappa=57827.0d0
  else if(r > R600 .and. r <= R400) then
    rho=11.2494d0-8.0298d0*x
    vpv=39.7027d0-32.6166d0*x
    vsv=22.3512d0-18.5856d0*x
    vph=vpv
    vsh=vsv
    Qmu=143.0d0
    Qkappa=57827.0d0
  else if(r > R400 .and. r <= R220) then
    rho=7.1089d0-3.8045d0*x
    vpv=20.3926d0-12.2569d0*x
    vsv=8.9496d0-4.4597d0*x
    vph=vpv
    vsh=vsv
    Qmu=143.0d0
    Qkappa=57827.0d0
  else if(r > R220 .and. r <= R80) then
!
! anisotropy in PREM only above 220 km
!
    rho=2.6910d0+0.6924d0*x
    vpv=0.8317d0+7.2180d0*x
    vph=3.5908d0+4.6172d0*x
    vsv=5.8582d0-1.4678d0*x
    vsh=-1.0839d0+5.7176d0*x
    eta_aniso=3.3687d0-2.4778d0*x
    Qmu=80.0d0
    Qkappa=57827.0d0

  else if(r > R80) then
    rho=2.6910d0+0.6924d0*x
    vpv=0.8317d0+7.2180d0*x
    vph=3.5908d0+4.6172d0*x
    vsv=5.8582d0-1.4678d0*x
    vsh=-1.0839d0+5.7176d0*x
    eta_aniso=3.3687d0-2.4778d0*x
    Qmu=600.0d0
    Qkappa=57827.0d0
  endif

  quick_prem_vsv = vsv
 

  end function quick_prem_vsv

