! test the linkage of fortran 90 programs with c subroutines
program testsacf90

implicit none
 

integer :: saclst_iheader_c,saclst_fheader_c,saclst_kheader_c,read_sac_float_c
character(len=100) filename,sta
integer npts,klen
real t0,delta,sacdata(10000),sacdata_bp(10000)
real trbdndw,flo,fhi,a
integer i,iord,passes


! read in header info
filename = 'gni.sac'
call read_sacfile_c(trim(filename)//char(0),t0,delta,npts,sacdata)
print *,' npts = ',npts, ' dt = ', delta, ' t0 = ', t0

! backup data
sacdata_bp(1:npts) = sacdata(1:npts)

! apply filter
trbdndw = 0.3
a = 30.
iord = 4
flo = 1/100.
fhi = 1/40.
passes = 2
call xapiir(sacdata,npts,'BU',trbdndw,a,iord,'BP',flo,fhi,delta,passes)

call write_sacfile_c(trim(filename)//char(0),'sac.tmp'//char(0),t0,npts,sacdata)

end program testsacf90
