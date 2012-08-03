!----------------------------------
  double precision function quick_prem_vsv(x)
  implicit none

  double precision, parameter ::    R_EARTH = 6371000.d0
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


! cut and pasted from prem_model.f90

  double precision x

  double precision r
  double precision rho, vpv, vsv, vph, vsh, Qmu, Qkappa, eta_aniso


  r=x*R_EARTH
  print *, x,r

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
  print *, 'DEBUG quick_prem_vsv=',quick_prem_vsv
 

  end function quick_prem_vsv



       program test_prem
       double precision, parameter ::    R_EARTH = 6371000.d0
       double precision :: rad_fraction
       double precision :: vsv


       rad_fraction=0.99
       print *, quick_prem_vsv(rad_fraction)
!      write(*,*) 'Rad=',rad_fraction, &
!      'depth=',R_EARTH*(1-rad_fraction)/1000.0, &
!      ' vsv=',vsv



       end program


