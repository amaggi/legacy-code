!=====================================================================
!
!  Alessia Maggi - February 2005
!
!  subroutines for reading regiostack isotropic vsfin models and
!  azimuthally anisotropic anvsfin models
! 
!  amp(G)=2*ro*betaV*[(2*sqrt((da2**2)+(da1**2))/c)]
!
!  Single depth version.
!
!=====================================================================

module three_d_mantle_model_constants

  implicit none
  include 'constants.h'

  ! set PREM interfaces
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
  integer, parameter :: NLAT = 180, NLON = 360
  double precision, parameter :: XPAS = 1.0d0
  double precision, parameter :: YPAS = 1.0d0
  double precision, parameter :: X0 = -180.0d0
  double precision, parameter :: Y0 = 0.0d0

  ! crustal model parameters for 3-SMAC
   integer, parameter :: NLAT_CRUST = 90
   integer, parameter :: NLON_CRUST = 180
   double precision, parameter :: XPAS_CRUST = 2.0d0
   double precision, parameter :: YPAS_CRUST = 2.0d0

  ! absolute or relative output
  logical, parameter :: OUTPUT_ABSOLUTE = .true.


end module three_d_mantle_model_constants

module three_d_mantle_model_variables

  use three_d_mantle_model_constants
  implicit none
  
  double precision grid_vs(1:NLAT,1:NLON)
  double precision grid_a1(1:NLAT,1:NLON)
  double precision grid_a2(1:NLAT,1:NLON)
  double precision grid_r(1:NLAT,1:NLON)
  ! note, z=0 refers to radius of moho

  double precision grid_phi(1:NLON)
  double precision grid_theta(1:NLAT)

  ! 3smac moho depth
  double precision smac_phi(1:NLON_CRUST)
  double precision smac_theta(1:NLAT_CRUST)
  double precision smac_moho(NLAT_CRUST,NLON_CRUST) 

  ! radius of information
  double precision radius

end module three_d_mantle_model_variables

!-----------------

  subroutine read_mantle_model

  use three_d_mantle_model_variables
  implicit none

  integer i,j,k,krad
  double precision depth 
  double precision get_theta, get_phi
  double precision dr, bg_rho, bg_vs 

  double precision vp,vs,rho

  ! set up the theta / phi arrays as read from vsfin file
  ! phi : -PI -> PI
  do i=1,NLAT
    grid_theta(i)=get_theta(i)
  enddo
  do i=1,NLON
    grid_phi(i)=get_phi(i)
  enddo


  ! read vsfin and anvsfin file
  open(unit=10,file='./vsfin',status='old')
  open(unit=11,file='./anvsfin',status='old')
  ! read only one depth slice from the files (the top one)
  read(10,*) depth
  read(11,*) depth
  radius = 1.0d0 - (depth*1000.d0)/R_EARTH

  ! read the absolute velocities and set the radius for this depth layer
  do j=1,NLAT
    read(10,*) (grid_vs(j,i),i=1,NLON) 
  enddo
  do j=1,NLAT
    read(11,*) (grid_a1(j,i),i=1,NLON) 
  enddo
  do j=1,NLAT
    read(11,*) (grid_a2(j,i),i=1,NLON) 
  enddo
  close(10)
  close(11)

  end subroutine read_mantle_model

!---------------------------


  subroutine mantle_model(theta,phi,dvs,alpha,am_phi,ampG)

  use three_d_mantle_model_variables
  implicit none

  double precision theta,phi,dvs,alpha,am_phi,ampG

  integer i,j,k
  integer iphi,jtheta
  double precision a1, a2
  double precision bg_vs, bg_rho
  double precision myphi,dr,dtheta,dphi
  double precision vsbox(1:4), a1box(1:4), a2box(1:4)
  double precision interpolate_box
  integer box_theta(1:4), box_phi(1:4)
  logical pole


  dvs = ZERO
  a1 = ZERO
  a2 = ZERO
  alpha = ZERO
  am_phi = ZERO
  ampG = ZERO
  pole = .false.

  ! find bracketing indexes for theta (i.e. the two indexes of the theta array
  ! that bracket the theta value requested)
  call bracket(grid_theta,NLAT,theta,jtheta)
  
  ! are we at the poles?
  if ((jtheta.eq.0) .or. (jtheta.eq.NLAT)) pole=.true.

  ! sort out the periodicity problem for phi
  ! do this on a copy myphi of original phi in case calling routine 
  ! wants to reuse the phi value
  if (phi>PI) then
    myphi=phi-TWO_PI
  else
    myphi=phi
  endif

  ! find bracketing indexes for phi
  call bracket(grid_phi,NLON,myphi,iphi)

  ! if we are in a polar region, then return the average of the
  ! points surrounding the polar cap, interpolating with depth
  if (pole) then 
    do i=1,NLON
      if (jtheta.eq.0) then
        ! north pole
        dvs=dvs+grid_vs(1,i)
        a1=a1+grid_a1(1,i)
        a2=a2+grid_a2(1,i)
      else
        ! south pole	  
        dvs=dvs+grid_vs(NLAT,i)
        a1=a1+grid_a1(NLAT,i)
        a2=a2+grid_a2(NLAT,i)
      endif
    enddo

    dvs=dvs/(1.0d0*NLON)
    a1=a1/(1.0d0*NLON)
    a2=a2/(1.0d0*NLON)

    ! get reference model values of beta and rho
    call quick_prem_rho_vsv(radius,bg_rho,bg_vs)
    ! alpha is peak to peak anisotropy value as a fraction
    alpha = 2*sqrt(a1**2+a2**2)/dvs
    ! calculate G in GPa (* 1000 to transform bg_vs from km/s to m/s)
    ampG=2*bg_rho*bg_vs*alpha*1000/10**9
    ! calculate angle (in radians) of fast SV direction 
    am_phi=atan2(a2,a1)/2.0
    ! if require relative velocities, then calculate them
    if(.not.OUTPUT_ABSOLUTE) then
      dvs = (dvs-bg_vs)/bg_vs
    endif

    return

  else
    ! we are not at a pole, so set up the latitude values of the 
    ! four corners of the interpolation box, and find interpolation
    ! distance dtheta
    box_theta(1)=jtheta
    box_theta(2)=jtheta+1
    box_theta(3)=jtheta+1
    box_theta(4)=jtheta
    dtheta=(theta - grid_theta(jtheta))/(grid_theta(jtheta+1) - grid_theta(jtheta))
  endif

  ! if we are at the longitude discontinuity, set up the longitude 
  ! values of the corners of the interpolation box to straddle the 
  ! discontinuity
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
    ! set up the corners of the interpolation box in the normal manner
    box_phi(1)=iphi
    box_phi(2)=iphi
    box_phi(3)=iphi+1
    box_phi(4)=iphi+1
    dphi=(myphi - grid_phi(iphi))/(grid_phi(iphi+1) - grid_phi(iphi))
  endif

   ! set up the interpolation box : vs
    vsbox(1)=grid_vs(box_theta(1),  box_phi(1))
    vsbox(2)=grid_vs(box_theta(2),  box_phi(2))
    vsbox(3)=grid_vs(box_theta(3),  box_phi(3))
    vsbox(4)=grid_vs(box_theta(4),  box_phi(4))

   ! set up the interpolation box : a1
    a1box(1)=grid_a1(box_theta(1),  box_phi(1))
    a1box(2)=grid_a1(box_theta(2),  box_phi(2))
    a1box(3)=grid_a1(box_theta(3),  box_phi(3))
    a1box(4)=grid_a1(box_theta(4),  box_phi(4))

   ! set up the interpolation box : a2
    a2box(1)=grid_a2(box_theta(1),  box_phi(1))
    a2box(2)=grid_a2(box_theta(2),  box_phi(2))
    a2box(3)=grid_a2(box_theta(3),  box_phi(3))
    a2box(4)=grid_a2(box_theta(4),  box_phi(4))



  ! interpolate over the 3D boxes

  dvs=interpolate_box(vsbox,dtheta,dphi)
  a1 =interpolate_box(a1box,dtheta,dphi) 
  a2 =interpolate_box(a2box,dtheta,dphi)

  ! get reference model values of beta and rho
  call quick_prem_rho_vsv(radius,bg_rho,bg_vs)
  ! alpha is peak to peak anisotropy value as a fraction
   alpha = 2*sqrt(a1**2+a2**2)/dvs
  ! calculate G in GPa (* 1000 to transform bg_vs from km/s to m/s)
  ampG=2*bg_rho*bg_vs*alpha*1000/10**9
  ! calculate angle (in radians) of fast SV direction 
  am_phi=atan2(a2,a1)/2.0
  ! if require relative velocities, then calculate them
  if(.not.OUTPUT_ABSOLUTE) then
    dvs = (dvs-bg_vs)/bg_vs
  endif

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

 double precision function interpolate_box(box,dtheta,dphi)
 
 double precision box(1:4)
 double precision dtheta, dphi
 
 interpolate_box =  (1-dtheta)*(1-dphi)*box(1) + &
                    (dtheta)  *(1-dphi)*box(2) + &
                    (dtheta)  *(dphi)  *box(3) + &
                    (1-dtheta)*(dphi)  *box(4) 

 
 end function interpolate_box

!----------------------------------


!----------------------------------
  subroutine quick_prem_rho_vsv(x,prem_rho,prem_vsv)

! Routine to calculate the prem Vsv values at a given depth from the
! original coefficients of the PREM model

  use three_d_mantle_model_variables
  implicit none

  double precision x

  double precision r,prem_rho,prem_vsv
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

  prem_rho = rho
  prem_vsv = vsv
  
  return
 
  end subroutine quick_prem_rho_vsv

! -------------------------------

  subroutine get_moho(theta, phi, moho)
    
  use three_d_mantle_model_variables
  double precision theta, phi, moho

  double precision myphi, dtheta, dphi, mohobox(4)
  integer i,jtheta, iphi, box_theta(4), box_phi(4)
  logical pole

  pole = .false.

  ! 3SMAC is defined between 0 and 360 degrees, so fix the
  ! longitude accordingly
  if (phi < 0) then
    myphi = phi + TWO_PI
  else
	myphi = phi
  endif

  ! bracket the latitude
  call bracket(smac_theta,NLAT_CRUST,theta,jtheta)
	
  ! are we at the pole?
  if ((jtheta.eq.0) .or. (jtheta.eq.NLAT_CRUST)) pole=.true.

  ! bracket the longitude
  call bracket(smac_phi,NLON_CRUST,myphi,iphi)

  ! if we are in a polar region, then return the average of the
  ! points surrounding the polar cap
  if (pole) then 
	moho=0
    ! pole region : return average of points surrounding polar cap
    ! interpolating for depth
    do i=1,NLON_CRUST
      if (jtheta.eq.0) then
        moho=moho+smac_moho(1,i)
      else
        moho=moho+smac_moho(NLAT_CRUST,i)
      endif
    enddo
    moho=moho/(1.0d0*NLON_CRUST)
    return
  else
    ! we are not at the pole, so set up the interpolation box coordinates
    box_theta(1)=jtheta
    box_theta(2)=jtheta+1
    box_theta(3)=jtheta+1
    box_theta(4)=jtheta
    dtheta=(theta - smac_theta(jtheta))/(smac_theta(jtheta+1) - smac_theta(jtheta))
  endif

  ! if we are at the longitude discontinuity, then set up the interpolation
  ! box coordinates accordingly
  if ((iphi.eq.0) .or. (iphi.eq.NLON_CRUST)) then
    box_phi(1)=NLON_CRUST
    box_phi(2)=NLON_CRUST
    box_phi(3)=1
    box_phi(4)=1
    if (iphi.eq.0) then
      dphi=1+(myphi - smac_phi(1))/(XPAS_CRUST*DEGREES_TO_RADIANS)
    else
      dphi=(myphi - smac_phi(NLON_CRUST))/(XPAS_CRUST*DEGREES_TO_RADIANS)
    endif
  else
    ! we can set up the longitude interpolation box coordinates in the normal manner
    box_phi(1)=iphi
    box_phi(2)=iphi
    box_phi(3)=iphi+1
    box_phi(4)=iphi+1
    dphi=(myphi - smac_phi(iphi))/(smac_phi(iphi+1) - smac_phi(iphi))
  endif

  ! set up the interpolation box
  mohobox(1)=smac_moho(box_theta(1),  box_phi(1))
  mohobox(2)=smac_moho(box_theta(2),  box_phi(2))
  mohobox(3)=smac_moho(box_theta(3),  box_phi(3))
  mohobox(4)=smac_moho(box_theta(4),  box_phi(4))

  ! perform the moho interpolation
  moho=(1-dtheta)*(1-dphi)*mohobox(1) + &
      (dtheta)  *(1-dphi) *mohobox(2) + &
      (dtheta)  *(dphi)   *mohobox(3) + &
      (1-dtheta)*(dphi)   *mohobox(4) 
 

  end subroutine get_moho

! -------------------------------
  subroutine read_moho

  use three_d_mantle_model_variables

  integer j,k

  open(unit=10,file="datacroute",status='old')
  read(10,*)
  read(10,*)
  read(10,*)
  do j=1,NLON_CRUST
    ! crust thick in m
    read(10,*)(smac_moho(k,j),k=1,NLAT_CRUST) 
  enddo
  close(10)


  do j=1,NLON_CRUST
    do k=1,NLAT_CRUST
      smac_moho(k,j)=(R_EARTH-smac_moho(k,j))/R_EARTH
    enddo
  enddo


  do j=1,NLON_CRUST
    smac_phi(j) = ((j-1) * XPAS_CRUST+XPAS_CRUST/2.0d0) * DEGREES_TO_RADIANS
  enddo
  do j=1,NLAT_CRUST
    smac_theta(j) = ((j-1) * YPAS_CRUST+YPAS_CRUST/2.0d0) * DEGREES_TO_RADIANS
  enddo
  
  end subroutine read_moho
  
  
