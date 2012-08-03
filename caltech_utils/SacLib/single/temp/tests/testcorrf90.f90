program testinf90

  implicit none

  integer,parameter :: NMAX = 10000
  real*8 :: t1(NMAX),t2(NMAX),corr(NMAX),corr2(NMAX),conv(NMAX)
  integer :: ind(NMAX),ndata,nsyn,npts_corr,zero_index,ios,i
  integer :: xcorr,crosscorr,autocorr,convolve
  real*8 :: t0,dt,t0_s,dt_s
 
  call get_ascfile_value("data.ascii",t0,dt,ndata,t1)
  call get_ascfile_value("syn.ascii",t0_s,dt_s,nsyn,t2);

  print *, 'Number of data : ',ndata
  print *, 'Number of syn : ', nsyn
  
  npts_corr = xcorr(t1,ndata,t2,nsyn,corr,ind)
  print *,'Total number of points for cross-correlation',npts_corr
  if (npts_corr == 0) stop 'Error in x correlation'
  
  open(11,file = 'testf90.corr1')
  do i = 1,npts_corr
     write(11,'(i7,g15.5)') ind(i),corr(i)
  enddo
  close(11)

  zero_index = crosscorr(t1,ndata,t2,nsyn,corr2)
  print *, 'Zero index of cross-correlation ', zero_index
  open(11,file = 'testf90.corr2')
  do i = 1,ndata+nsyn-1
     write(11,'(i7,g15.5)') i-zero_index,corr2(i)
  enddo
  close(11)

  npts_corr = xcorr(t1,ndata,t1,ndata,corr,ind)
  print *,'Total number of points for auto-correlation',npts_corr
  if (npts_corr == 0) stop 'Error in x correlation'
  
  open(11,file = 'testf90.autocorr1')
  do i = 1,npts_corr
     write(11,'(i7,g15.5)') ind(i),corr(i)
  enddo
  close(11)

  zero_index = autocorr(t1,ndata,corr2)
  print *, 'Zero index for auto-correlation', zero_index
  open(11,file = 'testf90.autocorr2')
  do i = 1,ndata+ndata-1
     write(11,'(i7,g15.5)') i-zero_index,corr2(i)
  enddo
  close(11)

  npts_corr = convolve(t1,ndata,t2,nsyn,conv)
  if (npts_corr == 0) stop 'Error in convolution'
  if (npts_corr > NMAX) stop 'Too large convolution series'
  print *, 'Total number of points for convolution', npts_corr
  open(11,file = 'testf90.conv')
  do i = 1,npts_corr
     write(11,'(g15.5 ,g15.5)') -15.95+(i-1)*0.05 ,conv(i)
  enddo
  close(11)

end program testinf90

  
