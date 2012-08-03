subroutine modes_R(filename, u, up ,nmod, nlay, wv, qv)
  implicit none
  character (len=50), intent(in) :: filename
  integer, intent(in) :: nlay, nmod
  integer :: jmode, a, j
  integer :: nv, lv
  real :: gv
  real, dimension(nmod), intent(out) :: wv, qv
  real, dimension(nmod,nlay), intent(out) :: u, up

  jmode = 1
    
  open(8, file=filename, form='unformatted', status='old')
  do jmode = 1, nmod
    read(8) nv, lv, wv(jmode), qv(jmode), gv, (u(jmode,j),j=1,nlay),&
             (up(jmode,j), j=1,nlay)
!    print *, nv, lv, wv, qv, gv
!    print *, (u(jmode,j),j=1,808), (up(jmode,j), j=1,808)
  end do
  close (8)
end