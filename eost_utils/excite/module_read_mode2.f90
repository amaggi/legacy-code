!********************************************************************
module read_mode
!********************************************************************
  implicit none
  integer, parameter :: DDP = kind(1.0d0)
  
  type model
    real :: r, rho
    type (model), pointer :: suivant
  end type model

  type mode
    real(DDP) :: fmhz
    character :: type
    integer :: nv, lv
    type (mode), pointer :: suivant
  end type mode

contains
!********************************************************************
  subroutine read_table (table, nlay, nmod, deb, deb2)
    implicit none
    character (len=50), intent(in) :: table
    integer, intent(out) :: nlay, nmod
    type (model), pointer :: deb, courant, courant1
    character (len=50) :: nom_mod
    integer, parameter :: fin_fich = -1
    integer :: err
    type (mode), pointer :: deb2, courant2, courant3
    real :: r, rho
    real(DDP) :: fmhz
    integer :: nv, lv
    character :: type
  
    nullify(deb)
    nullify(deb2)
    print *, 'DEBUG : table = ', table, '***'
    open(7,file=table,status='old',form='formatted')
    read(7,'(4x, a50/)') nom_mod
    read(7,'(//)')
    nlay = 1
    nmod = 1
    allocate(courant)
    read(7,'(9x, f9.1, 3x, f9.3)') r, rho
    courant%r = r
    courant%rho = rho
    deb => courant
    do 
      allocate(courant1)
      read(7,'(9x, f9.1, 3x, f9.3)') r, rho
      if (rho == 0) exit
      courant1%r = r
      courant1%rho = rho
      nullify(courant1%suivant)
      courant%suivant => courant1
      courant => courant1
      nlay = nlay+1
    end do 
    read(7,'(////)')  
    allocate(courant2)
    read(7,'(1x, i4, 1x, a1, 1x, i4, 17x, f13.7)',iostat=err) nv, type, lv, fmhz
    courant2%nv = nv
    courant2%type = type
    courant2%lv = lv
    courant2%fmhz = fmhz
    deb2 => courant2
    do
      allocate(courant3) 
      read(7,'(1x, i4, 1x, a1, 1x, i4, 17x, f13.7)',iostat=err) nv, type, lv, fmhz
      if (err == fin_fich) exit
      courant3%nv = nv
      courant3%type = type
      courant3%lv = lv
      courant3%fmhz = fmhz
      nullify(courant3%suivant)
      courant2%suivant => courant3
      courant2 => courant3
      nmod = nmod+1
    end do
    close (7)
  end subroutine read_table

!********************************************************************
  subroutine read_R(mode, nlay, u, up, ww, q) 
    implicit none
    integer, parameter :: DDP = kind(1.0d0)
    character (len=50), intent(in) :: mode
    integer, intent(in) :: nlay
    real(DDP), intent(out), dimension(nlay) :: u, up
    real, intent(out) :: ww, q
    integer :: n, l
    integer :: i, err
    integer, parameter :: fin_fich = -1
    real :: c, u1, up1
  
    open(9, file = mode, form = 'formatted', status = 'old')
    read(9, '(1x, i4, 1x, i4, 1x, e14.6, 1x, e14.6, 1x, e14.6)') n, l, ww, q, c
    do i= 1, nlay
      read(9, '(e12.5, 2x, e12.5)', iostat=err) u1, up1
      if(err == fin_fich) exit
      u(i)=u1
      up(i)=up1
    end do
    close(9)
  end subroutine read_R

!********************************************************************  
  subroutine read_S(mode, nlay, u, up, v, vp, ww, q) 
    implicit none
    integer, parameter :: DDP = kind(1.0d0)
    character (len=50), intent(in) :: mode
    integer, intent(in) :: nlay
    real(DDP), intent(out), dimension(nlay) :: u, up, v, vp
    real, intent(out) :: ww, q
    integer :: n, l
    integer :: i, err
    integer, parameter :: fin_fich = -1
    real :: c, u1, up1, v1, vp1
  
    open(10, file = mode, form = 'formatted', status = 'old')
    read(10, '(1x, i4, 1x, i4, 1x, e14.6, 1x, e14.6, 1x, e14.6)') n, l, ww, q, c
    do i= 1, nlay
      read(10, '(e12.5, 2x, e12.5, 2x, e12.5, 2x, e12.5)', iostat=err) u1, up1, v1, vp1
      if(err == fin_fich) exit
      u(i)=u1
      up(i)=up1
      v(i)=v1
      vp(i)=vp1
    end do
    close(10)
  end subroutine read_S

!*********************************************************************
  subroutine read_T(mode, nlay, w, wp, ww, q) 
    implicit none
    integer, parameter :: DDP = kind(1.0d0)
    character (len=50), intent(in) :: mode
    integer, intent(in) :: nlay
    real(DDP), intent(out), dimension(nlay) :: w, wp
    real, intent(out) :: ww, q
    integer :: n, l
    integer :: i, err
    integer, parameter :: fin_fich = -1
    real :: c, w1, wp1
  
    open(10, file = mode, form = 'formatted', status = 'old')
    read(10, '(1x, i4, 1x, i4, 1x, e14.6, 1x, e14.6, 1x, e14.6)') n, l, ww, q, c
    do i= 1, nlay
      read(10, '(56x, e12.5, 2x, e12.5)', iostat=err) w1, wp1
      if(err == fin_fich) exit
      w(i)=w1
      wp(i)=wp1
    end do
    close(10)
  end subroutine read_T

!********************************************************************  
  subroutine nom_mode(n,l,nom,c)
    implicit none
    integer, intent(in) :: n, l
    character, intent(in) :: c
    character (len = 7), intent(out) :: nom
    integer :: j, i, e
    integer, dimension(3) :: d1, d2
    integer, dimension(2,10) :: tab
  
    nom = '000' // c // '000'
    tab(1,:) = (/ (i,i=0,9) /)
    tab(2,:) = (/ (i,i=48,57) /)
  
    d1(1) = n/100
    e = n-d1(1)*100
    d1(2) = e/10
    d1(3) = e - d1(2)*10
  
    d2(1) = l/100
    e = l-d2(1)*100
    d2(2) = e/10
    d2(3) = e - d2(2)*10
  
    do j = 1,3
      do i = 0,9
        if(d1(j) == i) then
           nom(j:j) = achar(tab(2,i+1))
        end if
        if (d2(j) == i) then
           nom(j+4:j+4) = achar(tab(2,i+1))
        end if
      end do
    end do
  end subroutine nom_mode
  
end module read_mode
!*************************************************************