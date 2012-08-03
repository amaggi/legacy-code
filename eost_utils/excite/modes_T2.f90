subroutine modes_T(filename, w, wp, nmod, nlay, wv, qv)
  implicit none
  character (len=50), intent(in) :: filename
  integer, intent(in) :: nmod, nlay
  integer :: jmode, a, j
  integer :: nv, lv
  integer, parameter :: DDP = kind(1.0d0)
  real :: gv
  real, dimension(nmod), intent(out) :: wv, qv
  real, dimension(nmod,nlay) :: w, wp

  jmode = 1
    
  open(8, file=filename, form='unformatted', status='old')
  do jmode = 1, nmod
    read(8) nv, lv, wv(jmode), qv(jmode), gv, (w(jmode,j),j=1,nlay), (wp(jmode,j), j=1,nlay)
!    print *, nv, lv, wv, qv, gv
!    print *, (w(jmode,j),j=1,808), (wp(jmode,j), j=1,808)
  end do
  close (8)
end