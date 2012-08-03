c
c $Id: in_sac.f,v 1.1.1.1 2002/07/12 11:15:19 maggi Exp $
c $Log: in_sac.f,v $
c Revision 1.1.1.1  2002/07/12 11:15:19  maggi
c
c
c Revision 1.1  2002/05/23 10:28:31  maggi
c Initial revision
c
c
c
      subroutine in_sac(inp,grn,sacnam,sfx)
C     Input a SAC binary format file. Same function as in_pitsa

      parameter(MAXP=2050)
      character*(*) sacnam,sfx
      character*40 sacname,kztype
      dimension grn(*)
      include '../include/sizes.inc'
      include '../include/commons.inc'
	  
      common/pitsa/selev,xphy,xdist

      nchars=nchar(sacnam)
      if(sfx.ne.' ') sacname=sacnam(1:nchars)//'.'//sfx
c      if(sfx.ne.' ') then
c	  open(inp,file=pitsnam(1:nchars)//'.'//sfx,status='old',
c     &  form='formatted',iostat=ierr)
c        if(ierr.ne.0) then
c          write(*,'(a,a)') 
c     1		' error opening pitsa file ',pitsnam(1:nchars)//'.'//sfx
c          nscan=0
c          return
c        endif
c      else
c	  open(inp,file=pitsnam,status='old',
c     &  form='formatted',iostat=ierr)
c        if(ierr.ne.0) then
c          write(*,'(a,a)') 
c     1		' error opening pitsa file ',pitsnam
c          nscan=0
c          return
c        endif
c      endif

      CALL RSAC1(sacname,grn,nscan,begin,dt,MAXP,NERRR)

c      read(inp,'(a1)')  dum
c      read(inp,*)  jy,jmon,jday,jh,jm,sec
c      call juldat(jy,jmon,jday,jd)	
c      read(inp,'(a1)')  dum
c      read(inp,*) srate 
c      dt=1./srate
c      read(inp,'(a1)')  dum
c      read(inp,*) iy,mon,iday,ih,im,ss
c      call juldat(iy,mon,iday,id)

C     read in event_time(origin time), then calculate start_time(begin time)

      CALL GETIHV('IZTYPE',kztype,NERR)
      if(kztype.ne.'IO') write(*,*) 'Wrong reference time!'
      CALL GETNHV('NZYEAR',jy,NERR)
      CALL GETNHV('NZJDAY',jd,NERR)
      CALL GETNHV('NZHOUR',jh,NERR)
      CALL GETNHV('NZMIN', jm,NERR)
      CALL GETNHV('NZSEC',jsec,NERR)
      CALL GETNHV('NZMSEC',jmsec,NERR)
      sec=jsec+float(jmsec)/1000
      iy=jy
      id=jd
      ih=jh
      im=jm
      ss=sec
      call datetime(iy,id,ih,im,ss,begin)

c      read(inp,'(a1)')  dum
c      read(inp,*) slat
       CALL GETFHV('STLA',slat,NERR)
	  t1=90.-slat
	  
c      read(inp,'(a1)')  dum
c      read(inp,*) slon
       CALL GETFHV('STLO',slon,NERR)
          p1=slon
          if(slon.lt.0.) p1=p1+360.

c      read(inp,'(a1)')  dum
c      read(inp,*) selev
       CALL GETFHV('STEL',selev,NERR)
c      read(inp,'(a1)')  dum
c      read(inp,*)  izero
c      read(inp,'(a1)')  dum
c      read(inp,'(a4)') sta 
       CALL GETKHV('KSTNM',sta,NERR)  
c      read(inp,'(a1)')  dum
c      read(inp,'(a4)') chn
       CALL GETKHV('KCMPNM',chn,NERR)

c      read(inp,'(a1)')  dum
c      read(inp,*)  elat 
       CALL GETFHV('EVLA', elat, NERR) 
	  th=90.-elat

c      read(inp,'(a1)')  dum
c      read(inp,*)  elon
       CALL GETFHV('EVLO', elon, NERR)
          ph=elon
          if(ph.le.180) ph=elon+360.

c      read(inp,'(a1)')  dum
c      read(inp,*)  d0
       CALL GETFHV('EVDP',d0,NERR)
       d0=d0/1000
c      read(inp,'(a1)')  dum
c      read(inp,*)  izero

  

c      read(inp,'(a1)')  dum
c      read(inp,*) phy
       CALL GETFHV('AZ',xphy,NERR)
c      read(inp,'(a1)')  dum
c      read(inp,*)  zero
c      read(inp,'(a1)')  dum
c      read(inp,*) dist
c      read(inp,'(a1)')  dum
c      read(inp,*) xdist
c      read(inp,'(a1)')  dum
       CALL GETFHV('DIST', xdist, NERR)
c      read(inp,*) zero
c      read(inp,'(a1)')  dum
c      read(inp,*) zero
c      read(inp,'(a1)')  dum
c      read(inp,*) zero
c	  
c	  j=1
c 2    read(inp,*,end=88,err=99) grn(j)
c 	  j=j+1
c	  goto 2
c88    nscan=j-1
c
c      close(inp)

       if(NERRR.ne.0) goto 99
       return

99     write(*,'(a)') ' error reading input'
c      nscan=j-1
	  write(*,*) nscan,sta,chn,typ,ky,kd,kh,km,xss,dt
       write(*,*) d0,th,ph,jy,jd,jh,jm,sec
       write(*,*) dist,t1,p1
	  
       write(*,*) selev,phy
c	  
c      close(inp)
       return
	  end
