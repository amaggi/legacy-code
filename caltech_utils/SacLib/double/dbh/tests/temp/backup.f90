! test the linkage of fortran 90 programs with c subroutines
program testsacf90

implicit none
 

integer :: saclst_iheader_c,saclst_fheader_c,saclst_kheader_c,read_sac_float_c
character(len=100) filename,sta
integer npts,klen
real*8 t0,delta,sacdata(10000),sacdata_bp(10000)
real*8 trbdndw,flo,fhi,a
integer i,iord,passes


! read in header info
print *, '1'
filename = 'PAS.CI.BHZ'
call dread_sacfile_c(trim(filename)//char(0),t0,delta,npts,sacdata)
!print *,' npts = ',npts, ' dt = ', delta, ' t0 = ', t0

! backup data
sacdata_bp(1:npts) = sacdata(1:npts)

! apply filter
trbdndw = 0.3
a = 30.
iord = 4
flo = 1/100.
fhi = 1/10.
passes = 2
call xapiir(sacdata,npts,'BU',trbdndw,a,iord,'BP',flo,fhi,delta,passes)

call dwrite_sacfile_c(trim(filename)//char(0),'1.tmp'//char(0),t0,npts,sacdata)


filename = 'PAS.CI.BHZ'
print *, '2'
call dread_sacfile_c(trim(filename)//char(0),t0,delta,npts,sacdata)
flo = 1./35
fhi = 1./6
call rtrend(sacdata,npts)
call rmean(sacdata,npts)
!print *, 'input to program',npts,trbdndw,a,iord,flo,fhi,delta,passes
call xapiir(sacdata,npts,'BU',trbdndw,a,iord,'BP',flo,fhi,delta,passes)
call dwrite_ascfile_c('2.tmp'//char(0),t0,delta,npts,sacdata)

end program testsacf90


