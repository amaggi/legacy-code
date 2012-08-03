program testinf90

  implicit none

  integer,parameter :: NMAX = 5000
  real*8 :: t(NMAX),corr(NMAX),hil(NMAX),env(NMAX)
  integer :: ind(NMAX)

  character(len=100) :: filename

  integer :: read_sac_double,xcorr,hilbert,envelope
  integer :: npts,npts_corr,nhil,nenv,i

  print *,'Input SAC file name'
  !read(*,'(a)') filename
  filename = 'sacfile'
  
  npts = read_sac_double(trim(filename)//char(0),t);
  if (npts == 0) stop 'Incorrect number of points'
  print *,'Total number of points in time series is ',npts

  nhil = hilbert(npts,t,hil)
  print *,'NHIL : ',nhil
  if (nhil /= 0) stop 'Error in hilbert transform'
  
  nenv = envelope(npts,t,env)
  print *,'NENV : ',nenv
  if (nenv /= 0) stop 'Error taking envelope'

  print *,'Writing result to file'

  open(12,file = 'testf90.env')
  do i = 1,npts
     write(12,'(i7,g15.5,g15.5,g15.5)') i,t(i),hil(i),env(i)
  enddo
  close(12)

end program testinf90

  
