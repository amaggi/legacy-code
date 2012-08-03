! test the linkage of fortran 90 programs with c subroutines
program testsacf90

implicit none
 
character(len=100) filename,sta
integer npts,klen,nerr
real*8 t0,dt,sacdata(10000)
real*8 trbdndw,flo,fhi,a,tcut_start,tcut_end,f(4)
integer i,iord,passes
character(len=100):: pz_file,kpfrom(2),kpto(2)

! apply filter
trbdndw = 0.3
a = 30.
iord = 4
flo = 1/35.
fhi = 1/6.
passes = 2




! read in header info
print *, 'Data file'
filename = '9703873.CI.PAS.BHZ.sac'
tcut_start = -20
tcut_end = 170
call dread_sacfile_cut_c(trim(filename)//char(0),tcut_start,tcut_end,t0,dt,npts,sacdata)
print *,' npts = ',npts, ' dt = ', dt, ' t0 = ', t0
pz_file = 'CI_PAS_BHZ_19960201_30000101'
kpfrom(1) = 'POLEZERO'
kpfrom(2) = pz_file
kpto(1) = 'NONE'
kpto(2) = ' '
f(2) = 1./35
f(3) = 1./6
f(1) = f(2) * 0.8
f(4) = f(3) * 1.2
call transfer(sacdata,npts,dt,kpfrom,kpto,f,nerr)
if (nerr /= 0) stop 'Error transfering'

call xapiir(sacdata,npts,'BU',trbdndw,a,iord,'BP',flo,fhi,dt,passes)
call dwrite_sacfile_c(trim(filename)//char(0),'data.sac'//char(0),t0,npts,sacdata)

print *, '2'
filename = 'PAS.CI.BHZ'
call dread_sacfile_c(trim(filename)//char(0),t0,dt,npts,sacdata)
print *,' npts = ',npts, ' dt = ', dt, ' t0 = ', t0
call rtrend(sacdata,npts)
call rmean(sacdata,npts)
call xapiir(sacdata,npts,'BU',trbdndw,a,iord,'BP',flo,fhi,dt,passes)
call dwrite_sacfile_c(trim(filename)//char(0),'syn.sac'//char(0),t0,npts,sacdata)

end program testsacf90


