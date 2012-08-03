  program generate_LPO
  
  integer, parameter :: CO_PLANAR = 1
  integer, parameter :: CO_LINEAR = 2
  real, parameter :: BOX_SIZE = 8
  real, parameter :: X0 = 1
  real, parameter :: Y0 = 1
  real, parameter :: Z0 = 1
  real, parameter :: PI = 3.1415926

  real :: alignment_fraction, random
  real, dimension (:), allocatable :: x1,y1,z1,x2,y2,z2
  real :: theta_aligned, phi_aligned, theta_random, phi_random
  real :: theta, phi,x,y,z,x_norm,y_norm,z_norm
  
  logical :: aligned

  integer :: a_type
  integer :: i

  character*120 :: filename

  ! USER INTERFACE

  ! obtain number of crystals
  write(*,*) 'Give number of anisotropic crystals :'
  read (*,*) n_crystals

  ! obtain alignment type 
  write(*,*) 'Alignment type : (1) co-planar, (2) co-linear. '
  read(*,*) a_type

  ! obtain alignment parameters 
  if(a_type.eq.CO_PLANAR) then
    write(*,*) 'Give theta, phi for normal to plane (degrees)'
    read(*,*) theta_aligned, phi_aligned
  else if(a_type.eq.CO_LINEAR) then
    write(*,*) 'Give theta, phi for line direction (degrees)'
    read(*,*) theta_aligned, phi_aligned
  else
    stop 'Incorrect value for alignment type.'
  end if
  theta_aligned = theta_aligned*PI/180.0
  phi_aligned = phi_aligned*PI/180.0

  ! obtain alignment fraction
  write(*,*) 'Fraction of crystals for which alignement applies: [0-1] :'
  read(*,*)   alignment_fraction 

  if (alignment_fraction<0 .or. alignment_fraction>1) then
    stop 'Incorrect value for alignment fraction.'
  endif

  ! obtain output filename
  write(*,*) 'Output filename'
  read(*,'(a)') filename

  ! END USER INTERFACE

  ! allocate memory
  allocate(x1(n_crystals))
  allocate(x2(n_crystals))
  allocate(y1(n_crystals))
  allocate(y2(n_crystals))
  allocate(z1(n_crystals))
  allocate(z2(n_crystals))

  ! loop over crystals
  do i = 1, N_CRYSTALS

    ! is this crystal going to be aligned ?
    call random_number(random)
    if (random < alignment_fraction) then
      aligned = .true.
    else
      aligned = .false.
    endif

    ! position end point of crystal
    call random_number(random)
    x1(i) = X0+random*BOX_SIZE
    call random_number(random)
    y1(i) = Y0+random*BOX_SIZE
    call random_number(random)
    z1(i) = Z0+random*BOX_SIZE

    ! if not aligned, then create a random orientation
    if (.not. aligned) then
      call random_number(random)
      theta = random*PI/2.0
      call random_number(random)
      phi = random*2*PI
      call thetaphi2xyz(theta,phi,x2(i),y2(i),z2(i))
    else
      if(a_type.eq.CO_PLANAR) then
        ! create normal vector
        call thetaphi2xyz(theta_aligned,phi_aligned,x_norm,y_norm,z_norm)
        ! create a random vector
        call random_number(random)
        theta_random = random*PI
        call random_number(random)
        phi_random = random*2*PI
        call thetaphi2xyz(theta_random,phi_random,x,y,z)
        ! take cross product
        call crossproduct(x_norm,y_norm,z_norm,x,y,z,x2(i),y2(i),z2(i))
        ! normalise resulting vector to 1
        call normalise(x2(i),y2(i),z2(i))

      else if (a_type.eq.CO_LINEAR) then
        theta=theta_aligned
        phi=phi_aligned
        call thetaphi2xyz(theta,phi,x2(i),y2(i),z2(i))
      else
        stop 'Incorrect value for alignment type.'
      end if
    endif
  enddo

  ! translate x2 etc by x1 etc
  x2=x2+x1
  y2=y2+y1
  z2=z2+z1

  ! write output to file
  open(unit=11,file=filename,status='unknown')
  if(a_type.eq.CO_PLANAR) then
    write(11,*) '# CO_PLANAR (theta,phi) = ',theta_aligned*180.0/PI, phi_aligned*180/PI
  else if(a_type.eq.CO_LINEAR) then
    write(11,*) '# CO_LINEAR (theta,phi) = ',theta_aligned*180.0/PI, phi_aligned*180/PI
  else
    stop 'Incorrect value for alignment type.'
  end if
  write(11,*) '# Percentage aligned = ', alignment_fraction*100
  write(11,*) '# '
  write(11,*) '# x1   y1   z1 '
  write(11,*) '# x2   y2   z2 '
  do i = 1, N_CRYSTALS
    write(11,'(3(f6.3,1x))') x1(i), y1(i), z1(i) 
    write(11,'(3(f6.3,1x))') x2(i), y2(i), z2(i) 
    write(11,'(a)') '>'
  enddo

  close(11)
  
  ! deallocate memeory
  deallocate(x1,x2,y1,y2,z1,z2)
  
  
  end program generate_LPO

  subroutine thetaphi2xyz(theta,phi,x,y,z)
  real, intent(in) :: theta, phi
  real, intent(out) :: x,y,z
  x=sin(theta)*cos(phi)
  y=sin(theta)*sin(phi)
  z=cos(theta)
  end subroutine

  subroutine crossproduct(ax,ay,az,bx,by,bz,cx,cy,cz)
  real, intent(in) :: ax,ay,az,bx,by,bz
  real, intent(out) :: cx,cy,cz
  cx=ay*bz-az*by
  cy=az*bx-ax*bz  
  cz=ax*by-ay*bx
  end subroutine

  subroutine normalise(x,y,z)
  real, intent(inout) :: x,y,z
  real :: amplitude

  amplitude=sqrt(x*x+y*y+z*z)
  x=x/amplitude
  y=y/amplitude
  z=z/amplitude
  end subroutine
