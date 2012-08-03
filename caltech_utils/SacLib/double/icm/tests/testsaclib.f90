
program testsacf90

implicit none
 

integer :: saclst_iheader_c,dsaclst_fheader_c,saclst_kheader_c,dread_sacdata_c
character(len=100) :: filename,pzfile,sta,kpfrom(2),kpto(2)
integer npts,klen
real*8 delta,sacdata(10000),sacdata_bp(10000),f(4)
real*8 trbdndw,flo,fhi,a
integer i,iord,passes,nerr


! read in header info
filename = 'majo.sac'
pzfile = 'majo.pz'
if (saclst_iheader_c(trim(filename)//char(0),'npts'//char(0), npts) /=0) &
     stop 'Error npts'
if (dsaclst_fheader_c(trim(filename)//char(0),'delta'//char(0),delta) /=0) &
     stop 'Error delta'
if (saclst_kheader_c(trim(filename)//char(0),'kstnm'//char(0),sta,klen) &
     /=0) stop 'Error kstnm'
print *,' npts = ',npts, ' dt = ', delta, ' sta = ', sta

! read in data
if (dread_sacdata_c(trim(filename)//char(0),sacdata) /= npts) &
     stop 'Error reading sac file'

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

! write out result
open(11,file='sac.tmp')
do i = 1,npts
   write(11,*) i*delta-1181.2,sacdata(i)
enddo
close(11)

end program testsacf90
