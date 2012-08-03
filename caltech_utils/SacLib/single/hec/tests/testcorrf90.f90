program testinf90

  implicit none

  integer,parameter :: NMAX = 10000
  real*4 :: data(NMAX),syn(NMAX),corr(NMAX),corr2(NMAX),conv(NMAX)
  real*4 :: autocorr(NMAX),autocorr2(NMAX),data_hil(NMAX),data_env(NMAX)
  integer :: ind(NMAX),ndata,nsyn,npts_corr,zero_index,ios,i
  integer :: xcorr_c,crosscorr_c,autocorr_c,convolve_c
  integer :: hilbert_c,envelope_c
  real*4 :: t0,dt,t0_s,dt_s
 
! read sac file
  call read_sacfile_c("data.sac"//char(0),t0,dt,ndata,data)
  call read_sacfile_c("syn.sac"//char(0),t0_s,dt_s,nsyn,syn);
  if (abs(dt-dt_s) > 1.0e-4) stop 'Not same sampling inteval'

  print *, 'Number of data : ',ndata
  print *, 'Number of syn : ', nsyn

! xcorr
  npts_corr = xcorr_c(data,ndata,syn,nsyn,corr,ind)
  if (npts_corr == 0) stop 'Error in x correlation'
  print *,'Total number of points for cross-correlation',npts_corr

! correlate
  zero_index = crosscorr_c(data,ndata,syn,nsyn,corr2)
  print *, 'Zero index of cross-correlation ', zero_index

  call write_ascfile_c("corr1.ascii"//char(0), 0, dt,npts_corr,corr)
  call write_ascfile_c("corr2.ascii"//char(0),0, dt,ndata+nsyn-1,corr2)


  npts_corr = xcorr_c(data,ndata,data,ndata,autocorr,ind)
  print *,'Total number of points for auto-correlation',npts_corr
  if (npts_corr == 0) stop 'Error in x correlation'
 
  zero_index = autocorr_c(data,ndata,autocorr2)
  print *, 'Zero index for auto-correlation', zero_index

  call write_ascfile_c("auto1.ascii"//char(0),0, dt,npts_corr,autocorr)
  call write_ascfile_c("auto2.ascii"//char(0),0, dt,ndata+ndata-1,autocorr2)


! convolve
  npts_corr = convolve_c(data,ndata,syn,nsyn,conv)
  if (npts_corr == 0) stop 'Error in convolution'
  if (npts_corr > NMAX) stop 'Too large convolution series'
  print *, 'Total number of points for convolution', npts_corr
  
  call write_ascfile_c("conv.ascii"//char(0), t0 + t0_s, dt, npts_corr,conv)

! hilbert and envelope
  print *, 'Hilbert and envelope transform of data array'
  if ( hilbert_c(ndata,data,data_hil) /= 0) stop 'Error hilbert tf'
  if ( envelope_c(ndata,data,data_env) /= 0) stop 'Error env tr'
  print *, 'ndata = ', ndata
  call write_ascfile_c("hil.ascii"//char(0),t0,dt,ndata,data_hil)
  print *, 'ndata = ', ndata
  call write_ascfile_c("env.ascii"//char(0),t0,dt,ndata,data_env)

  
end program testinf90

  
