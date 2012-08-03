subroutine modes_S(filename, u, up, v, vp, nmod, nlay, wv, qv)
  implicit none
  character (len=50), intent(in) :: filename
  integer, intent(in) :: nmod, nlay
  integer :: jmode, a, j
  integer :: nv, lv
  integer, parameter :: DDP = kind(1.0d0)
  real, dimension(nmod), intent(out) :: wv, qv
  real, dimension(nmod,nlay), intent(out) :: u, up, v, vp
  real, dimension(nlay) :: p, pp
  real :: gv
 
  jmode = 1
    
  open(9, file=filename, form='unformatted', status='old')
  do jmode = 1, nmod
    read(9) nv, lv, wv(jmode), qv(jmode), gv,&
          (u(jmode,j),j=1,nlay), (up(jmode,j), j=1,nlay),&
          (v(jmode,j), j=1,nlay), (vp(jmode,j), j=1,nlay),&
          (p(j), j=1,nlay), (pp(j), j=1,nlay)
!    print *, nv, lv, wv, qv, gv
!    print *, (u(jmode,j),j=1,808), (up(jmode,j), j=1,808), &
!              (v(jmode,j), j=1,nlay), (vp(jmode,j), j=1,nlay), (p(j), j=1,nlay), (pp(j), j=1,nlay)
  end do
  close (9)
end