c
c $Id: plt_rtns.f,v 1.1.1.1 2002/07/12 11:15:20 maggi Exp $
c $Log: plt_rtns.f,v $
c Revision 1.1.1.1  2002/07/12 11:15:20  maggi
c
c
c Revision 1.1  2002/05/23 10:28:35  maggi
c Initial revision
c
c
      subroutine dspsms(string,ndata,idsp,b,a)
c **** display seismograms
      dimension b(*),a(*)
      character*(*) string 
      
      call init_plot

      if(idsp.gt.0) then
c  display data and synthetics
        call minmax(b,ndata,ymin,ymax)
        call minmax(a,ndata,zmin,zmax)
        ymin=amin1(ymin,zmin)
        ymax=amax1(ymax,zmax)
        if(ymin.eq.0..and.ymax.eq.0.) then
          ymin=-1.
          ymax=1.
        endif
        call tkplt(0,b,ndata,ymin,ymax,string)
        call tkplt(1,a,ndata,ymin,ymax,string)
      else
c  display only synthetics
        if(idsp.eq.0)call minmax(b,ndata,ymin,ymax)
        if(ymin.eq.0..and.ymax.eq.0.) then
          ymin=-1.
          ymax=1.
        endif
        call tkplt(0,b,ndata,ymin,ymax,string)
      endif
      call finitt
      return
      end

      subroutine window(dsdm,mrow,nlen,ncol,string)
      include '../include/sizes.inc'
      include '../include/commons.inc'
      character*(*) string 
      character*1 chr
      real*8 dsdm
      common/array/dat(lgrm),syn(lgrm)
      common/work/a(maxsect),b(lgrm)
      common/dat/ndata(maxrecs),isb(maxcol),add(maxcol),idat	 
      dimension dsdm(mrow,*)
      call minmax(dat,nlen,ymin,ymax)
      call minmax(syn,nlen,zmin,zmax)
      ymin=amin1(ymin,zmin)
      ymax=amax1(ymax,zmax)
      if(ymin.eq.0..and.ymax.eq.0.) then
        ymin=-1.
        ymax=1.
      endif
    6 ntot=nlen
      call tkplt(0,dat,ntot,ymin,ymax,string)
      call tkplt(1,syn,ntot,ymin,ymax,string)
      call finitt
	  
  902 write(*,'(/,a)') ' <s> to choose start, <e> to choose end: '
       call cursor(xstrt,yy,chr)
       if(chr.ne.'s'.and.chr.ne.'S') goto 902
	  
  905 call cursor(xstop,yy,chr)
       if(chr.ne.'e'.and.chr.ne.'E') then
         write(*,'(a)') ' choose end with <e>!'
         go to 905
       endif

       ntot=(xstop-xstrt)/dt+1
       ip=xstrt/dt
       if(ip.lt.0) ip=0
       ntot=min(ntot,nlen-ip)
       do 5 i=1,ntot
       ii=ip+i
       a(i)=dat(ii)
    5  b(i)=syn(ii)
       call tkplt(0,a,ntot,ymin,ymax,string)
       call tkplt(1,b,ntot,ymin,ymax,string)
       call finitt

       write(*,'(/,a)') ' okay? '
       read(*,'(a1)') chr 
       if(chr.eq.'n'.or.chr.eq.'N') go to 6

C      get the start_time	  
       call datetime(iy,id,ih,im,ss,xstrt)
       do 7 i=1,ntot
       dat(i)=a(i)
   7   syn(i)=b(i)
       nlen=ntot
       if(ip.ne.1) then
        do 8 i=1,ncol
        do  8 j=1,ntot
        jk=j+idat
        jj=ip+jk
   8   dsdm(jk,i)=dsdm(jj,i)
       endif
       return
       end

      subroutine tkplt(iplt,b,npt,ymin,ymax,string)
c  if iplt.gt.0 then don't clear screen
c  if iplt.eq.0 clear screen for immediate plot
      dimension b(*)
      include '../include/sizes.inc'
      include '../include/commons.inc'
      character*(*) string
      if(iplt.le.0) then
        call init_plot
        xlong=(npt-1)*dt
        dxtic=xlong*.01
        dxticb=dxtic*10.
        dytic=abs(ymax-ymin)*.01
        dyticb=dytic*10.
        call setfor(240.,.5,1.)
        call axis(8.,5.5,.5,1.,1.5,1.,xlong,0.,ymax,ymin,
     & dxtic,dxticb,dytic,dyticb,'(*)','(e9.2)','seconds',
     & 'counts',string,1)
c       xctr=xlong*.5
c       yctr=ymax+.2
        
c        call text(xctr,yctr,0.,3,string,0)
      endif
      if(iplt.eq.1) call setfor(0.,.5,1.)
      xo=0.
      yo=b(1)
      do i=2,npt
      x=(i-1)*dt
      call line(xo,yo,x,b(i),.05,0,0)
      xo=x
      yo=b(i)
      enddo
      return
      end
	  
       subroutine init_plot
c  initialize plot dimensions, etc.
c      character*40 psname
        ssize=.6
        xwin=0.
        ywin=0.

c       write(*,'(A)') 'Input the postscript name:'
c       read(*,*) psname
        call initt(1,'','','fitting',ssize,xwin,ywin)
c  clear the screen
        call clear
c  set foreground color to pale blue
        call setfor(240.,.3,1.)
c  set background color to red-green
        call setbac(60.,.5,1.)
c  set horizontal and vertical dimensions of plot,coords of lower left corner
        call setdim(10.,7.5,0.,0.)
c  set to no scaling of data values
        call setscl(0.,1.,0.,1.)
c  select a font 
        call setaxf(132)
        return
        end
