! test the linkage of fortran 90 programs with c subroutines
program testsacf90

implicit none

character(len=100) :: data_file,syn_file,pz_file,ext
real*8 :: tcut_start, tcut_end,f(4),trbdndw,a
character(len=30):: kpfrom(2), kpto(2)
integer :: passes,iord
 
real*8 :: data(10000), dt, syn(10000), dt_s, t0, t0_s
integer :: npts, npts_s, nerr

data_file = '9703873.CI.PAS.BHZ.sac'
syn_file = 'PAS.CI.BHZ'
pz_file = 'CI_PAS_BHZ_19960201_30000101'

tcut_start = -20
tcut_end = 170
ext = 'ext'
kpfrom(1) = 'POLEZERO'
kpfrom(2) = pz_file
kpto(1) = 'NONE'
kpto(2) = ' '
   
f(2) = 1./35
f(3) = 1./6
f(1) = f(2) * 0.8
f(4) = f(3) * 1.2
passes = 2
trbdndw = 0.3
a = 30.
iord = 4


!    read cutted data array
call dread_sacfile_cut_c(trim(data_file)//char(0),tcut_start,tcut_end,t0,dt,npts,data)
!call xapiir(data,npts,'BU',trbdndw,a,iord,'BP',f(2),f(3),dt,passes)

! write sac file is the line that causes the problem!
call dwrite_sacfile_c(trim(data_file)//char(0),"data.tmp"//char(0),t0,npts,data)

! read all syn file
call dread_sacfile_c(trim(syn_file)//char(0),t0_s,dt_s,npts_s,syn)
call dwrite_ascfile_c('syn.before'//char(0),t0_s,dt_s,npts_s,syn)
call xapiir(syn,npts_s,'BU',trbdndw,a,iord,'BP',f(2),f(3),dt_s,passes)
call dwrite_ascfile_c('syn.after'//char(0),t0_s,dt_s,npts_s,syn)

end program testsacf90



