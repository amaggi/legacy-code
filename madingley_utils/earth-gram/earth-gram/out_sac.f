c
c $Id: out_sac.f,v 1.1.1.1 2002/07/12 11:15:20 maggi Exp $
c
c $Log: out_sac.f,v $
c Revision 1.1.1.1  2002/07/12 11:15:20  maggi
c
c
c Revision 1.2  2002/05/23 11:50:45  maggi
c Removed all calculation of dates and times.  The sac files now start at o=b=0.  Event and station locations are put into the header, and distances and azimuths are calculated from them by SAC.
c
c Revision 1.1  2002/05/23 10:28:33  maggi
c Initial revision
c
c
      subroutine out_sac(ioutp,grn,sacnam,sfx)
C     Output a SAC binary format file. Same function as out_pitsa.

      character*(*) sacnam,sfx
      character*40 sacname
      dimension grn(*)
      include '../include/sizes.inc'
      include '../include/commons.inc'

      common/pitsa/selev,xphy,xdist
	  
      data izero,zero/0,0.0/
	  
      nchars=nchar(sacnam)
      sacname=sacnam(1:nchars)
      if(sfx.ne.' ') sacname=sacnam(1:nchars)//'.'//sfx
c      if(sfx.ne.' ') then
c	  open(ioutp,file=pitsnam(1:nchars)//'.'//sfx,status='unknown',
c     &  form='formatted',iostat=ierr)
c        if(ierr.ne.0) then
c          write(*,'(a,a)') 
c     1		' error opening pitsa file ',pitsnam(1:nchars)//'.'//sfx
c          return
c        endif
c      else
c	  open(ioutp,file=pitsnam,status='unknown',
c     &  form='formatted',iostat=ierr)
c        if(ierr.ne.0) then
c          write(*,'(a,a)') 
c     1		' error opening pitsa file ',pitsnam
c          return
c        endif
c      endif
c	  
c      if(ierr.ne.0) then
c        write(iflog,'(a,a)') 
c     1		' error opening pitsa file ',pitsnam
c	    return
c	  endif

c      srate=1./dt

c     Some debuggin statements
c      do i = 1 , nscan
c        write (*,*) i, grn(i)
c      enddo

       CALL NEWHDR
       CALL SETNHV('NPTS',nscan,NERR)
       CALL SETIHV('IFTYPE','ITIME',NERR)
       CALL SETLHV('LEVEN',.TRUE.,NERR)
c      write(ioutp,'(a)')  '#samp_freq'
c      write(ioutp,'(f10.6)') srate 
       CALL SETFHV('DELTA',dt,NERR)
C
C     set the event origin time as reference time, calculate the begin time
C
       CALL SETIHV('IZTYPE','IO',NERR)
       CALL SETFHV('B',0.,NERR)

c      call gredat(jy,jmon,jday,jd)	  
c      write(ioutp,'(a)')  '#event_time'
c      write(ioutp,'(i4,1x,4(i2,1x),f9.6)')
c     & jy,jmon,jday,jh,jm,sec
c      CALL SETNHV('NZYEAR',jy,NERR)
c      CALL SETNHV('NZJDAY',jd,NERR)
c      CALL SETNHV('NZHOUR',jh,NERR)
c      CALL SETNHV('NZMIN',jm,NERR)
c      jsec=int(sec)
c      CALL SETNHV('NZSEC',jsec,NERR)
c      jmsec=int((sec-jsec)*1000)
c      CALL SETNHV('NZMSEC',jmsec,NERR)
c      write(ioutp,'(a)')  '#start_time'
c      call gredat(iy,mon,iday,id)
c      write(ioutp,'(i4,1x,4(i2,1x),f9.6)')
c     & iy,mon,iday,ih,im,ss
c      istart_time=( ( ((iy-jy)*365+id-jd)*24 +ih-jh)*60+im-jm)*60
c      start_time=float(istart_time)+ss-sec
c      CALL SETFHV('B',start_time,NERR)
 
       slat=90.-t1
       slon=p1
       if(p1.gt.180) slon=p1-360.
c      write(ioutp,'(a)')  '#station_coor_1'
c      write(ioutp,'(f11.4)') slat
       CALL SETFHV('STLA',slat,NERR)
c      write(ioutp,'(a)')  '#station_coor_2'
c      write(ioutp,'(f11.4)') slon
       CALL SETFHV('STLO',slon,NERR)
c      write(ioutp,'(a)')  '#station_coor_3'
c      write(ioutp,'(f11.4)') selev
       CALL SETFHV('STEL',selev,NERR)
c      write(ioutp,'(a)')  '#station_coor_type'
c      write(ioutp,'(i1)')  izero
c      write(ioutp,'(a)')  '#station_code'
c      write(ioutp,'(a4)') sta   
       CALL SETKHV('KSTNM',sta,NERR)
c      write(ioutp,'(a)')  '#station_channel'
c      write(ioutp,'(a4)') chn
       CALL SETKHV('KCMPNM',chn,NERR)

       elat=90.-th
       elon=ph
       if(ph.gt.180) elon=ph-360.
c      write(ioutp,'(a)')  '#event_coor_1'
c      write(ioutp,'(f11.4)') elat
       CALL SETFHV('EVLA',elat,NERR) 
c      write(ioutp,'(a)')  '#event_coor_2'
c      write(ioutp,'(f11.4)') elon
       CALL SETFHV('EVLO',elon,NERR)
c      write(ioutp,'(a)')  '#event_coor_3'
c      write(ioutp,'(f11.4)') d0
       CALL SETFHV('EVDP',d0*1000,NERR)
c      write(ioutp,'(a)')  '#event_coor_type'
c      write(ioutp,'(i1)')  izero
c
c       CALL SETLHV('LCALDA','false',NERR)
c      write(ioutp,'(a)')  '#event_az'
c      write(ioutp,'(f11.4)') phy
c       CALL SETFHV('AZ',phy,NERR)
c      write(ioutp,'(a)')  '#event_back_az'
c      write(ioutp,'(f3.1)') zero
c      write(ioutp,'(a)')  '#event_hypo_dist'
c      write(ioutp,'(f11.4)') dist
c       CALL SETFHV('DIST',dist,NERR)
c      write(ioutp,'(a)')  '#event_epi_dist'
c      write(ioutp,'(f11.4)') dist
c      write(ioutp,'(a)')  '#event_local_mag'
c      write(ioutp,'(f3.1)') zero
c      write(ioutp,'(a)')  '#event_body_wave_mag'
c      write(ioutp,'(f3.1)') zero
c      write(ioutp,'(a)')  '#event_surface_wave_mag'
c      write(ioutp,'(f3.1)') zero
c	  
c      do j=1,nscan
c        write(ioutp,'(g15.7)',err=99) grn(j)
c      enddo
c
c      close(ioutp)
      CALL WSAC0(sacname,XDUMMY,grn(1),NERR)
      if(NERR.ne.0) goto 99
	  return

99    write(*,'(a)') ' error writing greens functions'
      return
	  end
