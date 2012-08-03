
program testsacf90

implicit none
 

integer :: saclst_iheader_c,saclst_fheader_c,saclst_kheader_c,read_sac_float_c
character(len=100) :: filename,pzfile,sta,kpfrom(2),kpto(2)
integer npts,klen
real t0,delta,sacdata(10000),sacdata_bp(10000),f(4)
real trbdndw,flo,fhi,a
integer i,iord,passes,nerr


! read in header info
filename = 'majo.sac'
pzfile = 'majo.pz'
call read_sacfile_c(trim(filename)//char(0),t0,delta,npts,sacdata)
print *,' npts = ',npts, ' dt = ', delta, ' t0 = ', t0


! backup data
sacdata_bp(1:npts) = sacdata(1:npts)

kpfrom(1) = 'POLEZERO'
kpfrom(2) = pzfile

kpto(1) = 'NONE'
kpto(2) = ''

f(1) = 1/120.0
f(2) = 1/100.0
f(3) = 1/40.0
f(4) = 1/35.0

call rtrend(sacdata,npts)
call rmean(sacdata,npts)
call transfer(sacdata,npts,delta,kpfrom,kpto,f,nerr)
if (nerr /= 0) stop 'Error transfering'

call write_sacfile_c(trim(filename)//char(0),'tran.sac'//char(0),t0,npts,sacdata)

end program testsacf90
