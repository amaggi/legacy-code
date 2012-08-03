!=====================================================================
!
!  Alessia Maggi - February 2005
!
!  subroutines for reading regiostack isotropic vsfin models and
!  azimuthally anisotropic anvsfin models
! 
!  amp(G)=2*ro*betaV*[(2*sqrt((da2**2)+(da1**2))/c)]
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

  ! value of smoothed PREM at 40km depth - for fixing up shallow mantle model
  double precision, parameter ::    PREM_40 = 4.400

  ! regiostack model at 1x1 degrees 
  integer, parameter :: NLAT = 180, NLON = 360, NLAY = 20
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
  
  double precision grid_vs(1:NLAT,1:NLON,0:NLAY)
  double precision grid_a1(1:NLAT,1:NLON,0:NLAY)
  double precision grid_a2(1:NLAT,1:NLON,0:NLAY)
  double precision grid_r(1:NLAT,1:NLON,0:NLAY)
  ! note, z=0 refers to radius of moho

  double precision grid_phi(1:NLON)
  double precision grid_theta(1:NLAT)

  ! 3smac moho depth
  double precision smac_phi(1:NLON_CRUST)
  double precision smac_theta(1:NLAT_CRUST)
  double precision smac_moho(NLAT_CRUST,NLON_CRUST) 

  ! actual number of layers
  integer n_layers

end module three_d_mantle_model_variables

!-----------------

  subroutine read_mantle_model

  use three_d_mantle_model_variables
  implicit none

  integer i,j,k,krad
  double precision depth, radius
  double precision get_theta, get_phi
  double precision dr, bg_rho, bg_vs, smooth_prem(0:NLAY)

  double precision vp,vs,rho,moho

  ! read the moho model, because we will need the moho depths to
  ! shift the vs-depths as required
  call read_moho

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
  do k=1,NLAY
    read(10,*,end=666) depth
    read(11,*,end=666) depth
    radius = 1.0d0 - (depth*1000.d0)/R_EARTH
    ! read the absolute velocities and set the radius for this depth layer
    do j=1,NLAT
      read(10,*) (grid_vs(j,i,k),i=1,NLON) 
      do i=1,NLON
        grid_r(j,i,k) = radius
     enddo
    enddo
    do j=1,NLAT
      read(11,*) (grid_a1(j,i,k),i=1,NLON) 
    enddo
    do j=1,NLAT
      read(11,*) (grid_a2(j,i,k),i=1,NLON) 
    enddo
    !print *, 'tt = ', depth, ' a1(159,33): ', grid_a1(159,33,k)
    !print *, 'tt = ', depth, ' a2(159,33): ', grid_a2(159,33,k)

  enddo
666 continue
  n_layers=k
  print *, 'Read ', n_layers, ' layers in model.'
  close(10)
  close(11)


  ! check for moho depth, and if necessary correct grid_r
  do j=1,NLAT
    do i=1,NLON

      ! set layer 0 values to smoothed PREM at 40 km depth 
      grid_vs(j,i,0) = PREM_40
      grid_a1(j,i,0) = 0.0
      grid_a2(j,i,0) = 0.0

      ! set layer 0 radius to moho depth from 3SMAC
      call get_moho(grid_theta(j), grid_phi(i),moho)
      grid_r(j,i,0) = moho

      ! moho could be deeper than the layer1 (40km) or layer2 (50km)
      ! if so, we need to shift those layers deeper by splitting the
      ! difference between the moho depth and the first layer deeper than
      ! the moho depth

      ! find the first unaffected layer index krad
      krad=0
      do while (grid_r(j,i,krad) >= grid_r(j,i,0))
        krad=krad+1
      enddo

      ! fix up all layers above the unaffected layer
      do k=1,krad-1
        ! if k=1, k-1=0, so moho depth is used to correct 40km layer
        ! if k=2, k-1=1, so corrected 40km is used to correct 50km 
        ! layer and so on
        grid_r(j,i,k) = (grid_r(j,i,k-1) + grid_r(j,i,krad))/2.0
      enddo

    enddo
  enddo

  end subroutine read_mantle_model

!---------------------------


  subroutine mantle_model(radius,theta,phi,dvs,alpha,am_phi,ampG)

  use three_d_mantle_model_variables
  implicit none

  double precision radius,theta,phi,dvs,alpha,am_phi,ampG

  integer i,j,k, kr1, kr2
  integer iphi,jtheta,krad
  double precision dvs1,dvs2
  double precision a1, a2, a1_1,a1_2,a2_1,a2_2
  double precision bg_vs, bg_rho
  double precision myphi,dr,dtheta,dphi
  double precision r_moho,r_model_bottom, r_40
  double precision vsbox(1:8), a1box(1:8), a2box(1:8)
  double precision r_tmp(0:NLAY)
  double precision interpolate_cube
  integer box_theta(1:4), box_phi(1:4)
  logical pole


  dvs = ZERO
  alpha = ZERO
  am_phi = ZERO
  ampG = ZERO
  pole = .false.

  ! find the bottom of the model (the top of the model will be found later)
  r_model_bottom = grid_r(1,1,n_layers) 

  ! if we are below the range of the model, return
  if( (radius < r_model_bottom)) return

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

  ! set the moho depth for this lat/lon point
  r_moho = grid_r(jtheta,iphi,0) 
  r_40 = grid_r(jtheta,iphi,1)
  
  ! if we are above the mantle model at this lat/lon value, return
  if (radius > r_40) return

  ! write a temporary array of layer depths for use in later interpolation
  ! use n_layers+1 because n_layers refers to the number of layers in the 
  ! regiostack model, however we have added a 0th layer to simulate the 
  ! moho for correct interpolation
  r_tmp(0:n_layers)=grid_r(jtheta,iphi,0:n_layers)
  
  ! bracket the depth index : this will return the two indexes of the depth
  ! array that bracket the requested depth given the distribution of layer depths
  ! at jtheta, iphi.  Here and in the following we assume that the lateral variation
  ! of radii is locally smooth, and much smaller than the separation between the layer
  ! depths, so that we can ignore differences in layer depth between neighbouring 
  ! lat/lon points 
  call bracket(r_tmp,n_layers+1,radius,krad)
  
  ! recheck for radius limits
  ! note that krad was returned assuming array was 1 based, but
  ! depth array is now 0 based, so we need to correct krad
  if ((krad.eq.n_layers+1) .or. krad.eq.0) then
    ! we are below the maximum depth or above the moho - so return 
    ! we checked for this earlier, so we should never have to get here anyway ...
    return
  else
    ! set the bracketing indexes and find interpolation distance dr
    ! krad = 1 means the first value of the 0 based array, kr=0
    kr1=krad-1
    kr2=kr1+1
    dr=(radius - grid_r(jtheta,iphi,kr1))/ &
       (grid_r(jtheta,iphi,kr2) - grid_r(jtheta,iphi,kr1))
  endif

  ! if we are in a polar region, then return the average of the
  ! points surrounding the polar cap, interpolating with depth
  if (pole) then 
    dvs1=0
    dvs2=0
    a1_1=0
    a1_2=0
    a2_1=0
    a2_2=0
    
    do i=1,NLON
      if (jtheta.eq.0) then
        ! north pole
        dvs1=dvs1+grid_vs(1,i,kr1)
        dvs2=dvs2+grid_vs(1,i,kr2)
        a1_1=a1_1+grid_a1(1,i,kr1)
        a1_2=a1_2+grid_a1(1,i,kr2)
        a2_1=a2_1+grid_a2(1,i,kr1)
        a2_2=a2_2+grid_a2(1,i,kr2)
      else
        ! south pole	  
        dvs1=dvs1+grid_vs(NLAT,i,kr1)
        dvs2=dvs2+grid_vs(NLAT,i,kr2)
        a1_1=a1_1+grid_a1(NLAT,i,kr1)
        a1_2=a1_2+grid_a1(NLAT,i,kr2)
        a2_1=a2_1+grid_a2(NLAT,i,kr1)
        a2_2=a2_2+grid_a2(NLAT,i,kr2)
      endif
    enddo

    dvs1=dvs1/(1.0d0*NLON)
    dvs2=dvs2/(1.0d0*NLON)
    a1_1=a1_1/(1.0d0*NLON)
    a1_2=a1_2/(1.0d0*NLON)
    a2_1=a2_1/(1.0d0*NLON)
    a2_2=a2_2/(1.0d0*NLON)

    dvs=(1-dr)*dvs1 + dr*dvs2
    a1=(1-dr)*a1_1 + dr*a1_2
    a2=(1-dr)*a2_1 + dr*a2_2

    ! get reference model values of beta and rho
    call quick_prem_rho_vsv(radius,bg_rho,bg_vs)
    ! alpha is peak to peak anisotropy value as a fraction
    alpha = 2*sqrt(a1**2+a2**2)/dvs
    ! calculate G in GPa (* 1000 to transform bg_vs from km/s to m/s)
    ampG=2*bg_rho*bg_vs*alpha*1000/10**9
    ! calculate angle (in radians) of fast SV direction 
    am_phi=atan2(2*a2/(alpha*dvs),2*a1/(alpha*dvs))/2.0
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

   ! set up the interpolation cubes : vs
    vsbox(1)=grid_vs(box_theta(1),  box_phi(1),  kr1)
    vsbox(2)=grid_vs(box_theta(2),  box_phi(2),  kr1)
    vsbox(3)=grid_vs(box_theta(3),  box_phi(3),  kr1)
    vsbox(4)=grid_vs(box_theta(4),  box_phi(4),  kr1)
    vsbox(5)=grid_vs(box_theta(1),  box_phi(1),  kr2)
    vsbox(6)=grid_vs(box_theta(2),  box_phi(2),  kr2)
    vsbox(7)=grid_vs(box_theta(3),  box_phi(3),  kr2)
    vsbox(8)=grid_vs(box_theta(4),  box_phi(4),  kr2)

   ! set up the interpolation cubes : a1
    a1box(1)=grid_a1(box_theta(1),  box_phi(1),  kr1)
    a1box(2)=grid_a1(box_theta(2),  box_phi(2),  kr1)
    a1box(3)=grid_a1(box_theta(3),  box_phi(3),  kr1)
    a1box(4)=grid_a1(box_theta(4),  box_phi(4),  kr1)
    a1box(5)=grid_a1(box_theta(1),  box_phi(1),  kr2)
    a1box(6)=grid_a1(box_theta(2),  box_phi(2),  kr2)
    a1box(7)=grid_a1(box_theta(3),  box_phi(3),  kr2)
    a1box(8)=grid_a1(box_theta(4),  box_phi(4),  kr2)

   ! set up the interpolation cubes : a2
    a2box(1)=grid_a2(box_theta(1),  box_phi(1),  kr1)
    a2box(2)=grid_a2(box_theta(2),  box_phi(2),  kr1)
    a2box(3)=grid_a2(box_theta(3),  box_phi(3),  kr1)
    a2box(4)=grid_a2(box_theta(4),  box_phi(4),  kr1)
    a2box(5)=grid_a2(box_theta(1),  box_phi(1),  kr2)
    a2box(6)=grid_a2(box_theta(2),  box_phi(2),  kr2)
    a2box(7)=grid_a2(box_theta(3),  box_phi(3),  kr2)
    a2box(8)=grid_a2(box_theta(4),  box_phi(4),  kr2)



  ! interpolate over the 3D cubes

  dvs=interpolate_cube(vsbox,dr,dtheta,dphi)
  a1 =interpolate_cube(a1box,dr,dtheta,dphi) 
  a2 =interpolate_cube(a2box,dr,dtheta,dphi)

  ! get reference model values of beta and rho
  call quick_prem_rho_vsv(radius,bg_rho,bg_vs)
  ! alpha is peak to peak anisotropy value as a fraction
   alpha = 2*sqrt(a1**2+a2**2)/dvs
  ! calculate G in GPa (* 1000 to transform bg_vs from km/s to m/s)
  ampG=2*bg_rho*bg_vs*alpha*1000/10**9
  ! calculate angle (in radians) of fast SV direction 
  am_phi=atan2(2*a2/(alpha*dvs),2*a1/(alpha*dvs))/2.0
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

 double precision function interpolate_cube(cube,dr,dtheta,dphi)
 
 double precision cube(1:8)
 double precision dr, dtheta, dphi
 
 interpolate_cube = (1-dtheta)*(1-dphi)*(1-dr)*cube(1) + &
                    (dtheta)  *(1-dphi)*(1-dr)*cube(2) + &
                    (dtheta)  *(dphi)  *(1-dr)*cube(3) + &
                    (1-dtheta)*(dphi)  *(1-dr)*cube(4) + &
                    (1-dtheta)*(1-dphi)*  (dr)*cube(5) + &
                    (dtheta)  *(1-dphi)*  (dr)*cube(6) + &
                    (dtheta)  *(dphi)  *  (dr)*cube(7) + &
                    (1-dtheta)*(dphi)  *  (dr)*cube(8) 

 
 end function interpolate_cube

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
  
  
