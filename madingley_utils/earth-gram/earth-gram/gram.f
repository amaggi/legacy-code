c 
c $Id: gram.f,v 1.1.1.1 2002/07/12 11:15:19 maggi Exp $
c $Log: gram.f,v $
c Revision 1.1.1.1  2002/07/12 11:15:19  maggi
c
c
c Revision 1.2  2002/05/23 11:50:27  maggi
c Added debug statements.
c
c Revision 1.1  2002/05/23 10:28:31  maggi
c Initial revision
c
c
      program grams
c   surface wave greens function program
c   to be used with a moment tensor as specified for normal modes.
c   greens functions are computed for a step function source with
c   scalar moment equal to 10x27 dyne-cms.
c   greens functions will be output
      implicit integer*4(i-n)
      include '../include/sizes.inc'
      include '../include/commons.inc'
      include '../include/units.inc'
	  
      character outfile*27,string*8,chr*1
      real*8 c1,c2,s1,s2,caz,saz

      common/strain/strn(lgrm,9)
      dimension seis(lgrm)
      common/pitsa/selev,xphy,xdist
      character ko*1,sfx*3,stacode*4
      character*3 cmp(9)
      character*2 chan(9)
      data rad/57.29578d0/
      data cmp(1),cmp(2),cmp(3)/'err','epp','ezz'/
      data cmp(4),cmp(5),cmp(6)/'epz','erz','erp'/
      data cmp(7),cmp(8),cmp(9)/'vrt','rad','trn'/
      data chan(1),chan(2),chan(3)/'rr','pp','zz'/
      data chan(4),chan(5),chan(6)/'pz','rz','rp'/
      data chan(7),chan(8),chan(9)/'v ','r ','t '/

c initial variables:
c  mdmin & mdmx are min. & max. modes to include in the synthetic
      mdmin=0
      mdmx=1
7     init=1
c
      write(*,'(a$)') 
     &   ' enter 1 to model data,2 for synthetics only: '
      read(*,*)idat
	  	   
      if(idat.eq.1) then
c *** section to model data on disk ***
c get information from input file
        write(*,'(a$)') ' enter data file name (preceding suffix): '
        read(*,'(a27)') outfile
        write(*,'(a$)') ' enter file suffix: '
	read(*,'(a3)') sfx
c read in data and header (includes receiver location and source info)
        call in_sac(iinf1,seis,outfile,sfx)
		 
	call read_source(iinf5)
        stacode=sta
        call get_station(stacode,iinf5)

c compute source & receiver dependent parameters
        call distaz(th,ph,t1,p1,dist,phy,c1,c2,s1,s2,caz,saz)
c start time relative to origin time
        to=((((iy-jy)*365.+(id-jd))*24.+(ih-jh))*60.+(im-jm))*60.+
     +   ss-sec
	init=2
        write(*,'(a)') ' Pitsa-file Header'
        call prnt_hdr
       endif

c *** section to create synthetics when data are not on disk
c     or modify data file parameters ***
      call sins(c1,s1,c2,s2,phy,caz,saz)
      if(nscan.eq.0) go to 9999

      write(*,103) sta,chn,dist,phy
  103 format(/,'*** sta. ',a4,1x,a4,': dist.(km)=',f12.5,
     & ' source-sta. az=',f10.5)
 
       call strnmk(iinf2,c1,s1,c2,s2)
       if(npts.eq.0)go to 9999

c write synthetic to output file
      write(*,'(/,a$)') ' enter output file name: '
      read(*,'(a27)') outfile
      write(*,'(/,a$)') ' omit straingram output? '
      read(*,'(a1)') ko
      mstart=1
      if(ko.ne.'n'.and.ko.ne.'N') mstart=7
	  
      string(1:4)=sta
      string(5:5)='/'
      write(*,'(/,a)')
     &   ' Enter an <n> if you do not want to save the synthetic'
      do j=mstart,9
      write (*,*) 'Must start writing files now'
      call minmax(strn(1,j),nscan,ymin,ymax)
c      Overflow
       if(ymax-ymin.le.2e-15) then
       ymax=(ymax-ymin)/2+1e-15
       ymin=(ymax-ymin)/2-1e-15
       endif
      string(6:8)=cmp(j)
c     call tkplt(0,strn(1,j),nscan,ymin,ymax,string) 
c     call cursor(xx,yy,chr)
c     call finitt
c     if(chr.ne.'n'.and.chr.ne.'N') then
        chn(3:4)=chan(j)
c	some debugging output:
        write (*,*) 'Calling out_sac'
        call out_sac(iouf1,strn(1,j),outfile,cmp(j))
c     endif
      enddo

c **** get next data file on disk or change synthetics  ****
      init=init+1
      write(*,'(/,a$)') 
     & ' want to procede with the next set of straingrams? '
      read(*,'(a1)') ko
      if(ko.eq.'y'.or.ko.eq.'Y') go to 7
 9999 call finitt
      call hdkild
      stop
      end
 

      subroutine sins(c1,s1,c2,s2,phy,caz,saz)
c  initialize or alter variables
      implicit integer*4(i-n)
      include '../include/sizes.inc'
      include '../include/commons.inc'
      include '../include/units.inc'
	  
      real*8 c1,c2,s1,s2,caz,saz
	  
      character ans*1,stacode*4
	  
      data rad/.0174532/
c set default values
      if(init.eq.1)then
        ivel=1
	jy=1995
	jd=180
	jh=11
	jm=0
	sec=0.
	call read_source(iinf5)

       dt=.5
       sta='TPNV'
       stacode='TPNV'
       typ='nsn '
       chn='bbv'
       to=40.
       iy=1995
       id=180
       ih=0
       im=0
       ss=40.
       call get_station(stacode,iinf5)

       call distaz(th,ph,t1,p1,dist,phy,
     &           c1,c2,s1,s2,caz,saz)
	   
       nscan=512
      endif

      index=-1
  910 write(*,'(a)') ' '
      write(*,104) phy,dist
  104 format(' 1:    source-sta. azimuth= ',f7.3,' dist.= ',f10.5,'km')
      write(*,105) sta
  105 format(' 2:    instrument at station ',a4)
      write(*,106) nscan
  106 format(' 3:    npts in straingram is ',i5)
      write(*,107) dt
  107 format(' 4:    sample interval is ',f4.1)
      write(*,108) d0
  108 format(' 5:    source depth is ',f6.2)
      write(*,109) iy,id,ih,im,ss
  109 format(' 6:    start (yr,day,hrs,mins,secs)=',i4,':',i3,2(':',
     + i2),':',f4.1)
      write(*,110) to
  110 format(' 7:    start time:',f7.2)
      write(*,112)   sig,del,gam
112   format(' 8:    strike,dip,slip=',3(f8.3,1x))
      write(*,113)   (sol(j),j=1,6)
113   format('         moment tensor=',6(f8.3,1x))
      smoment=fmom*1e27
      write(*,114) smoment
114   format(' 9:    scalar moment=',g15.7)
      write(*,115) tconst
115   format(' 10:    source time constant=',f7.3)
      if(ivel.gt.0) then
        write(*,116) 
116     format(' 11:    output seismograms will be velocity (km/sec)')
      else
        write(*,117) 
117     format(' 11:    output seismograms will be displacement (cm)')
      endif

      write(*,'(a)') ' '
	  
      if(index.eq.0) return
      write(*,'(a$)') ' okay? '
      read(*,'(a1)') ans
      if(ans.ne.'n'.and.ans.ne.'N') return
	  
  900 write(*,'(a$)') ' enter index of parameter to change(0 quits): '
      read(*,*)index
      if(index.eq.0) go to 910
	  
      if(index.eq.1)then
      write(*,'(a$)') ' want to enter azimuth,dist.: '
      read(*,'(a1)') ans

      if(ans.eq.'n'.or.ans.eq.'N')then
        write(*,'(a$)')
     +  ' enter source lat.,lon. (0,0 to use disk file): '
        read(*,*) th,ph
	if(th.eq.0..and.ph.eq.0.) then
          write(*,'(a$)') 
     +     ' enter  yr dy hr mn sec near origin time ' 
          read(*,*) jy,jd,jh,jm,sec
	  call read_source(iinf5)
	else
          th=90.-th
          if(ph.lt.0.) ph=360.+ph
	endif
        write(*,'(a$)') '  enter new station code: '
        read(*,'(a4)') stacode
	call get_station(stacode,iinf5)
		
        call distaz(th,ph,t1,p1,dist,phy,c1,c2,s1,s2,caz,saz)
       else
         write(*,'(a$)') ' enter source-station azimuth,dist.: '
         read(*,*) phy,dist
         phl=(180.-phy)*rad
         c1=cos(phl)
         s1=sin(phl)
         c2=c1*c1-s1*s1
         s2=2.d0*c1*s1
       endif
	   
      elseif(index.eq.2)then
        write(*,'(a$)') '  enter station code: '
        read(*,'(a4)') stacode
	call get_station(stacode,iinf5)
		
        call distaz(th,ph,t1,p1,dist,phy,c1,c2,s1,s2,caz,saz)

      elseif(index.eq.3)then
       write(*,'(a$)') ' enter new straingram npts: '
       read(*,*)nscan

      elseif(index.eq.4)then
       write(*,'(a$)') ' enter sample interval: '
       read(*,*)dt

      elseif(index.eq.5)then
       write(*,'(a$)') ' enter new source depth '
       read(*,*)d0

      elseif(index.eq.6)then
       write(*,'(a$)') ' enter new <yr:day:hrs:mins:secs>: '
       read(*,*) iy,id,ih,im,ss
       jdays=0
       if(iy.ne.jy) call juldat(jy,12,31,jdays)	  
       to= ((float(id+jdays-jd)*24.+
     1     float(ih-jh))*60.+float(im-jm))*60.+(ss-sec)

      elseif(index.eq.7) then
       write(*,'(a$)') ' enter increment to start time: '
       read(*,*) tnew
       call update(tnew)
      elseif(index.eq.8) then
        write(*,'(a$)') ' do you want to use strike,dip and slip ? '
        read(*,'(a1)') ans
        if(ans.eq.'y'.or.ans.eq.'Y') then
          call fault
        else
          write(*,'(a$)') ' enter moment tensor elements f(x6) : '
          read(*,*) (sol(j),j=1,6)
	endif
	  
      elseif(index.eq.9) then
	   write(*,'(a$)') ' enter scalar moment: '
	   read(*,*) xmom
	   fmom=xmom*1.d-27

      elseif(index.eq.10) then
	   write(*,'(a$)') ' enter source time constant: '
	   read(*,*) tconst

      elseif(index.eq.11) then
           ivel=-ivel
      endif

      go to 900

      return
      end

      subroutine update(tnew)
c subroutine corrects start time and data series for new start time
      implicit integer*4(i-n)
      include '../include/sizes.inc'
      include '../include/commons.inc'

      to=to+tnew
      ss=ss+tnew
      if(ss.ge.60.)then
       imn=ss/60.
       ss=ss-float(imn)*60.
       im=im+imn
      endif
      if(im.ge.60)then
       jhr=im/60
       im=im-jhr*60
       ih=ih+jhr
      endif
      if(ih.ge.24)then
       idy=ih/24
       ih=ih-idy*24
       id=id+idy
      endif
      if(id.ge.366)then
       iyr=id/366
       id=id-iyr*366
       iy=iy+iyr
      endif
      return
      end


       subroutine strnmk(idfl,c1,s1,c2,s2)
      implicit integer*4(i-n)
      include '../include/sizes.inc'
      include '../include/commons.inc'
      real*8 c1,c2,s1,s2
      character*1 vary
	  
      lsou=idfl
      call startup(vary,dist,distr,distv,
     &   lsou,lrec,lref,lvar) 
      if(vary.eq.'e') then
        call mkhomog(lsou,c1,s1,c2,s2)
      else
        call mkhetero(vary,lsou,lrec,lref,lvar,
     &      distr,distv,c1,s1,c2,s2)
      endif
      return
      end


       subroutine startup(vary,dist,distr,distv,
     &                  lsou,lrec,lref,lvar)


       character filename*64, vary*1
	   
 1       write(*,'(/,a)') ' vary structure at/along:'
	 write(*,'(a)') '   <s> source end,'
	 write(*,'(a)') '   <r> receiver end,'	   
	 write(*,'(a)') '   <m> middle of path,'
	 write(*,'(a$)') '   <e> entire path: '
	 read(*,'(a1)') vary
	 if(vary.eq.'E') then
	   vary='e'
	 elseif(vary.eq.'M') then
  	   vary='m'
	 elseif(vary.eq.'R') then
	   vary='r'
	 elseif(vary.eq.'S') then
	   vary='s'
	 endif
	 if(vary.ne.'s'.and.vary.ne.'r'.
     &      and.vary.ne.'m'.and.vary.ne.'e') goto 1

         if(vary.ne.'e') then
2          write(*,'(/,a$)') ' enter length of path to vary: '
	   read(*,*) distv
	   if(distv.ge.dist) then
	     write(*,'(a)') ' length exceeds total path length!'
	     goto 2
	   else
	     distr=dist-distv
	   endif
	 else
	   distr=0.
	   distv=dist
	 endif
	   
       if(vary.eq.'s') then
         lrec=lsou+1
         lref=lrec
         lvar=lsou
         write(*,'(/,a)') ' Enter Rayleigh file name(s): '
	 write(*,'(a$)') '  name of source end greens function file: '
	 read(*,'(a64)') filename
	 open(lsou,file=filename,form='unformatted',status='old')
	 write(*,'(a$)') '  name of receiver end greens function file: '
	 read(*,'(a64)') filename
	 open(lrec,file=filename,form='unformatted',status='old')
         write(*,'(/,a)') ' Enter Love file name(s): '
	 write(*,'(a$)') '  name of source end greens function file: '
	 read(*,'(a64)') filename
	 open(lsou+2,file=filename,form='unformatted',status='old')
	 write(*,'(a$)') '  name of receiver end greens function file: '
	 read(*,'(a64)') filename
	 open(lrec+2,file=filename,form='unformatted',status='old')
		 
       elseif(vary.eq.'r') then
         lrec=lsou+1
         lref=lsou
         lvar=lrec
         write(*,'(/,a)') ' Enter Rayleigh file name(s): '
	 write(*,'(a$)') '  name of receiver end greens function file: '
	 read(*,'(a64)') filename
	 open(lrec,file=filename,form='unformatted',status='old')
	 write(*,'(a$)') '  name of source end greens function file: '
         read(*,'(a64)') filename
	 open(lsou,file=filename,form='unformatted',status='old')
         write(*,'(/,a)') ' Enter Love file name(s): '
	 write(*,'(a$)') '  name of receiver end greens function file: '
	 read(*,'(a64)') filename
	 open(lrec+2,file=filename,form='unformatted',status='old')
	 write(*,'(a$)') '  name of source end greens function file: '
         read(*,'(a64)') filename
	 open(lsou+2,file=filename,form='unformatted',status='old')
		 
       elseif(vary.eq.'m') then
         lrec=lsou
         lref=lsou
         lvar=lsou+1
         write(*,'(/,a)') ' Enter Rayleigh file name(s): '
	 write(*,'(a$)') '  name of middle section greens function file: '
	 read(*,'(a64)') filename
	 open(lvar,file=filename,form='unformatted',status='old')
	 write(*,'(a$)') '  name of receiver/source ends greens function file: '
	 read(*,'(a64)') filename
	 open(lref,file=filename,form='unformatted',status='old')
         write(*,'(/,a)') ' Enter Love file name(s): '
	 write(*,'(a$)') '  name of middle section greens function file: '
	 read(*,'(a64)') filename
	 open(lvar+2,file=filename,form='unformatted',status='old')
	 write(*,'(a$)') '  name of receiver/source ends greens function file: '
	 read(*,'(a64)') filename
	 open(lref+2,file=filename,form='unformatted',status='old')

       elseif(vary.eq.'e') then
         write(*,'(/,a$)')
     &  ' Enter Rayleigh greens function file name(s): '
	 read(*,'(a64)') filename
	 open(lsou,file=filename,form='unformatted',status='old')
         write(*,'(/,a$)') 
     &   ' Enter Love greens function file name(s): '
	 read(*,'(a64)') filename
         llove=lsou+2
	 open(llove,file=filename,form='unformatted',status='old')
       endif
       
       write(*,'(a)') ' '
       return
       end
