c
c $Id: peripheral.f,v 1.1.1.1 2002/07/12 11:15:20 maggi Exp $
c $Log: peripheral.f,v $
c Revision 1.1.1.1  2002/07/12 11:15:20  maggi
c
c
c Revision 1.1  2002/05/23 10:28:34  maggi
c Initial revision
c
c
      subroutine minmax(buf,n,rmin,rmax)
      dimension buf(*)
      rmin=0.
      rmax=0.
      do 5 i=1,n
      rmin=amin1(rmin,buf(i))
      rmax=amax1(rmax,buf(i))
    5 continue
      return
      end subroutine minmax

      subroutine get_station(stacode,iunit)
c **  subroutine returns colatitude, longitude, and station type
c     for station stacode.
      include '../include/sizes.inc'
      include '../include/commons.inc'

      character*4 stacode
      logical yeah
      character*52 filename,sfpath
  
c     filename='/home/solaris/xuelin/sfdata/stations'
      call getenv('SFDIR',sfpath)
      filename=sfpath(:lnblnk(sfpath))//'/stations'
1     inquire(file=filename,exist=yeah)
      if(yeah) then
        open(iunit,file=filename,status='old',iostat=ierr)
        if(ierr.ne.0) then
          write(*,'(a,a52)') ' error opening file ',filename
c         goto 7
        endif 
      else
        write(*,'(a)') ' station file ',filename,' does not exist!'
7       write(*,'(a)') 
     &   ' enter new file name or blank for interactive entry: '
        read(*,'(a52)') filename
        if(filename.eq.' ') goto 4
        goto 1
      endif

2     read(iunit,*,end=3) sta,t1,p1,elev,typ
      if(sta.eq.stacode) goto 5
      goto 2
   
5     t1=90.-t1
      if(p1.lt.0.) p1=360.+p1
      close(iunit)
      return
       
3     write(*,'(a)') ' station ',stacode,' not found'
      close(iunit)
   
4     write(*,'(a$)') 
     &        '  enter station type<sro ,asro,dwws,rstn,ida,nsn >: '
      read(*,'(a3)') typ
      write(*,'(a$)') ' enter station lat.,lon.: '
      read(*,*) t1,p1
      t1=90.-t1
      if(p1.lt.0.) p1=360.+p1

      return
      end subroutine get_station

      subroutine read_source(iunit)
c  routine reads source location information from file
c  or interactively
      implicit integer*4(i-n)
      include '../include/sizes.inc'
      include '../include/commons.inc'
      character*52 filename,sfpath
      logical yeah
      dimension xsol(6)
  
      xfound=99999.
  
c     filename='/home/solaris/xuelin/sfdata/sources'
      call getenv('SFDIR',sfpath)
      filename=sfpath(:lnblnk(sfpath))//'/sources'
1     inquire(file=filename,exist=yeah)
      if(yeah) then
        open(iunit,file=filename,status='old',iostat=ierr)
        if(ierr.ne.0) then
          write(*,'(a,a52)') ' error opening file ',filename
          goto 2
        endif
      else
        write(*,'(a,a52,a)') ' source file ',filename,' does not exist!'
        write(*,'(a$)') 
     &   ' enter new file name or blank for interactive entry: '
        read(*,'(a52)') filename
        if(filename.eq.' ') goto 2
        goto 1
      endif

      read(iunit,*)
  
   3  read(iunit,*,end=4) iyy,idd,ihh,imm,siec,xlat,xlon,xd0,
     & xstk,xdip,xslip,(xsol(j),j=1,6),xmom,xconst,qlength

      idif=( ( ((iyy-jy)*365+idd-jd)*24 +ihh-jh)*60+imm-jm)*60
      dif= abs(float(idif)+siec-sec)

      if(dif.lt.xfound) then
        xfound=dif
	ky=iyy
	kd=idd
	kh=ihh
	km=imm
	skec=siec
	th=90.-xlat
	ph=xlon
        if(xlon.lt.0.) ph=360.+xlon
	d0=xd0
	sig=xstk
	del=xdip
	gam=xslip
	fmom=xmom*1.d-27
        tconst=xconst
	do j=1,6
	  sol(j)=xsol(j)
	enddo
	xlength=qlength
      endif
		
      goto 3

	  
2     write(*,'(a)') ' enter source parameters,e.g.'
      write(*,'(a)') '  yr  dy hr mn sec     lat.   long.   depth  ' 
      write(*,'(a)') ' 1978 66  2 48 47.6   32.00  -137.61    12 '
      read(*,*) jy,jd,jh,jm,sec,xlat,xlon,d0
	  th=90.-xlat
	  ph=xlon
      if(xlon.lt.0.) ph=360.+xlon
      write(*,'(a,a)') 
     & '  strike dip slip or  ',
     & 'mzz mxx myy mxz myz mxy, mo   rise_time  length' 
      write(*,'(a,a)') 
     & '    30    50  90       ',
     & '0   0   0   0   0   0   1e27   1         10.'
      read(*,*)  sig,del,gam,(sol(j),j=1,6),fmom,tconst,xlength
      fmom=fmom*1.d-27
      goto 5	  
	  
   4  close(15)
          jy=ky
          jd=kd
          jh=kh
          jm=km
          sec=skec
   5  if(sol(1).eq.0.d0.and.sol(2).eq.0.d0.and.sol(3).eq.0.d0.
     & and.sol(4).eq.0.d0.and.sol(5).eq.0.d0.and.sol(6).eq.0.d0)    
     & call udc(sig,del,gam,sol)

      return
      end


      subroutine distaz(th,ph,tr,pr,dist,phy,c1,c2,s1,s2,caz,saz)
c subroutine computes source-receiver distance and azimuth
      implicit integer*4(i-n),real*8(a-h,o-z)
      real*4 th,ph,t1,p1,dist,phy,tr,pr
      data rad/57.29578d0/,ckm/111.195d0/
	  t1=tr/rad
	  p1=pr/rad
      t0=th/rad
      p0=ph/rad
      t2=t1
      c0=dcos(t0)
      s0=dsin(t0)
      c1=dcos(t2)
      s1=dsin(t2)
      dp=p1-p0
      co=c0*c1+s0*s1*dcos(dp)
      si=dsqrt(1.d0-co*co)
c   caz,saz refer to reciever-source azimuth measured cw from north
      saz=-s0*dsin(dp)/si
      caz=(c0-co*c1)/(si*s1)
      dist=(datan(si/co))*rad
      if(dist.le.0.d0)dist=dist+180.d0
      dist=dist*ckm
      s1=s1*dsin(dp)/si
      c1=(c0*co-c1)/(s0*si)
c   c1,s1 refer to source-reciever azimuth measured cw from north
      phy=datan(s1/c1)*rad
      if(c1.lt.0.d0) phy=phy+180.d0
      phy=180.d0-phy
      c2=c1**2-s1**2
      s2=2.d0*c1*s1
      return
      end

      subroutine udc(sig,del,gam,f)
c
c  udc computes the unit normal, unit slip and unit moment tensor
c  given strike(sig), dip(del) and slip(gam) in degrees.
c
      implicit real*8(a-h,o-z)
      dimension un(6),us(6),f(*)
      data rad/57.29578/
      s=sig/rad
      d=del/rad
      g=gam/rad
      cs=dcos(s)
      ss=dsin(s)
      cd=dcos(d)
      sd=dsin(d)
      cg=dcos(g)
      sg=dsin(g)
      un(1)=cd
      un(2)=sd*ss
      un(3)=sd*cs
      us(1)=-sg*sd
      us(2)=sg*cd*ss-cg*cs
      us(3)=sg*cd*cs+cg*ss
      f(1)=2.*un(1)*us(1)
      f(2)=2.*un(2)*us(2)
      f(3)=2.*un(3)*us(3)
      f(4)=un(1)*us(2)+us(1)*un(2)
      f(5)=un(1)*us(3)+us(1)*un(3)
      f(6)=un(2)*us(3)+us(2)*un(3)
      return
      end
	  
	  subroutine find_sfx(filename,sfx)
      character*(*) filename,sfx
      nlen=len(filename)
c  be sure blanks are removed from the end of the filename
	  do i=1,nlen
	  if(filename(i:i).eq.' ') goto 2
	  enddo
	  goto 3
2     nlen=i-1

c  check for the last '.' in the file name
3	  do i=1,nlen
	  ii=nlen+1-i
	  if(filename(ii:ii).eq.'.') goto 1
	  enddo
c  file name has no suffix
      sfx='   '
	  return
c  file name has a suffix
1	  nlen=min(nlen,ii+4)
      sfx=filename(ii+1:nlen)
	  filename(ii:nlen)=' '
	  return
	  end

      integer function nchar(pitsnam)
	  character*(*) pitsnam
	  nlen=len(pitsnam)
	  do i=1,nlen 
	  nchar=i
	  if(pitsnam(i:i).eq.' ') goto 2
	  enddo
      nchar=nlen
	  return
2     nchar=nchar-1
      return
	  end

      subroutine inresp(nt)
c  compute instrument response
      implicit integer*4(i-n)
       include '../include/sizes.inc'
       include '../include/commons.inc'
      real*8 df,pi,sh
      complex evres,slresp,irisresp
      data pi/3.14159265358979d0/
c      call pzall(nsta,nchn,ntyp,iy,id)
      nh=nt/2
      nha=nh+1
      sh=dt
      df=2.d0*pi/(sh*nt)
      do 5 i=2,nha
      wr=(i-1)*df
c  excitations are for displacements in cm, gdsn gains
c  are counts/meter and responses are for acceleration
c    5 resp(i)=.01*wr*wr*evres(wr)
      resp(i)=cmplx(1.,0.)

c  excitations are for displacements in cm, some IRIS(eg. SLR)
c  gains are counts/meter and responses are for displacement.
c  but some other IRIS(eg. SUR LBTB BOSA) gains are counts(meter/s)
c  and responses are for velocity. in irisresp, the responses have been
c  converted to displacement responses with unit counts/meter,
c  which is equivalent to resp(i)=.01*i*wr*velresp.
      if(typ.eq.'iris') resp(i)=.01*irisresp(wr)

c  excitations are for displacements in cm, SLR gains
c  are counts/meter and responses are for displacements
   5  if(typ.eq.'dwws') resp(i)=.01*slresp(wr)
      return
      end

      subroutine disply
c  subroutine to display options for synthetics
      implicit integer*4(i-n)
       include '../include/sizes.inc'
       include '../include/commons.inc'
	  
   99 write(*,'(a)') ' '
      write(*,100) sdep(nsdp)
  100 format(' 1:    source depth is ',f6.2)
      write(*,101) (i,sdep(i),i=1,nsrce)
  101 format(8x,'available depths are:',/,
     &  4(9x,3(i2,1x,f6.2,2x),/))
      write(*,102) mdmin,mdmx
  102 format(' 2:    mode',i3,1x,'to mode',i3,' included')
      write(*,103) nbran
  103 format(8x,i3,' modes may be included')

      write(*,'(/,a$)') ' enter index to change (0 quits): '
      read(*,*) index
      if(index.eq.0) return
      if(index.eq.1) then
       write(*,'(a$)') ' enter index of source depth desired: '
       read(*,*) nsdp
      elseif(index.eq.2) then
       write(*,'(a$)') ' enter min. and max. modes to include: '
       read(*,*)mdmin,mdmx
      endif
      go to 99
      end

      subroutine fault
      implicit real*8(a-h,o-z)
       include '../include/sizes.inc'
       include '../include/commons.inc'
      write(*,'(a$)') ' enter strike, dip and slip in degrees : '
      read(*,*) sig,del,gam
      call udc(sig,del,gam,sol)
      return
      end


      function nfac(n0)
c
c $$$$$$ calls only library routines $$$$$$
c
c  nfac returns the largest integer le n0 which will be acceptable to
c  singleton's fft routine (1966 version).  see the comments at the
c  start of fft for details
c   programmed on 28 aug 1979 by ray buland
c
      dimension npr(9),ipr(9)
      data npr/2,3,5,7,11,13,17,19,23/
      n = iabs(n0)
      if(n.le.1) go to 13
c  initialize the prime factor search
 4    k = n
      nn = 0
      do 5 i=1,9
 5    ipr(i) = 0
c  extract all prime factors up to 23
      do 1 i=1,9
 3    l = k/npr(i)
      if(l*npr(i)-k.ne.0) go to 1
      k = l
      nn = nn + 1
      ipr(i) = mod(ipr(i)+1,2)
      if(k-1) 2,2,3
 1    continue
c  fall through if the largest prime factor is greater than 23
 6    n = n - 1
      go to 4
c  try again if more than 208 prime factors have been found
 2    if(nn.ge.209) go to 6
      k = 1
c  check the square free portion
      do 7 i=1,9
      if(ipr(i).gt.0) k = k*npr(i)
 7    continue
c  the product of the square free factors must by less than 211
      if(k.gt.210) go to 6
 13   nfac = n
      return
      end


      subroutine srcget
c subroutine finds source depth in input file nearest to desired depth
       include '../include/sizes.inc'
       include '../include/commons.inc'
      dnear=9999.
      do 9 i=1,nsrce
      dif=abs(d0-sdep(i))
      if(dif.gt.dnear) go to 9
      nsdp=i
      dnear=dif
    9 continue
      return
      end
