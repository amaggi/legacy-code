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

  npts_corr = xcorr(t,npts,t,npts,corr,ind)
  print *,'Total Number of points for auto-correlation',npts_corr
  if (npts_corr == 0) stop 'Error in x correlation'
  
  nhil = hilbert(npts,t,hil)
  print *,'NHIL : ',nhil
  if (nhil /= 0) stop 'Error in hilbert transform'
  
  nenv = envelope(npts,t,env)
  print *,'NENV : ',nenv
  if (nenv /= 0) stop 'Error taking envelope'

  print *,'Writing result to file'

  open(11,file = 'testf90.corr')
  do i = 1,npts_corr
     write(11,'(i7,g15.5)') ind(i),corr(i)
  enddo
  close(11)
  
  open(12,file = 'testf90.env')
  do i = 1,npts
     write(12,'(i7,g15.5,g15.5,g15.5)') i,t(i),hil(i),env(i)
  enddo
  close(12)

end program testinf90

  
