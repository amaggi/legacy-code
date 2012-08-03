c
c $Id: mkhetero.f,v 1.1.1.1 2002/07/12 11:15:20 maggi Exp $
c $Log: mkhetero.f,v $
c Revision 1.1.1.1  2002/07/12 11:15:20  maggi
c
c
c Revision 1.1  2002/05/23 10:28:32  maggi
c Initial revision
c
c
       subroutine mkhetero(vary,lsou,lrec,lref,lvar,
     &          distr,distv,c1,s1,c2,s2)
c  subroutine reads in output from the excitation program.
c  equations are taken from mendiguren (1977) but the signs of
c  imaginary parts are opposite to those given.  the radial component
c  also is multiplied by -1.  greens functions have been tested
c  against those computed from a normal mode sum for azimuths in all
c  four quadrants.  differences with mendiguren possibly stem from the
c  different coordinate system used in the computation of the excitations
c  and ambiguities in mendiguren's paper.  this program does agree with aki
c  and richards derivation (1981).
      implicit integer*4(i-n)
      include '../include/sizes.inc'
      include '../include/commons.inc'
      include '../include/units.inc'
      real*8 sddep,w,ccc,gmmm,py1,py2,py3,ddt,ddtv
      real*8 c1,c2,s1,s2,b1,b2,gb1,gb2
      character*1 vary
	  
      common/strain/strn(lgrm,9)
      complex u,wu,dudr,dudz,dwdr,dwdz,erz,gr
      complex v,dvdr,dvdz,erp,gl
      dimension sddep(lsd),py(3)
      data degrad,pi/.017453293,3.141592653535897933/

      c2p1=0.5d0*(1.d0+c2)
      c2m1=0.5d0*(1.d0-c2)
      s22=0.5d0*s2
	   
      read(lref) nsrce,npts,ddt,jcom,nbran
      read(lvar) nsrcv,nptv,ddtv,jcomv,nbranv
      if(nsrce.ne.nsrcv) then
        write(*,'(/,a,2i5)')
     &  ' nsrce differs between files, nsrce=',nsrce,nsrcv
        lref=lref+2
        lvar=lvar+2
        npts=0
        goto 1111
      elseif(nptv.ne.npts) then
        write(*,'(/,a,2i5)')
     &  ' npts differs between files, npts=',npts,nptv
        lref=lref+2
        lvar=lvar+2
        npts=0
        goto 1111
      elseif(ddt.ne.ddtv) then
        write(*,'(/,a,f8.4,1x,f8.4)')
     &  ' dt differs between files, dt=',ddt,ddtv
        lref=lref+2
        lvar=lvar+2
        npts=0
        goto 1111
      elseif(jcom.ne.jcomv) then
        write(*,'(/,a,2i5)')
     &  ' jcom differs between files, jcom=',jcom,jcomv
        lref=lref+2
        lvar=lvar+2
        npts=0
        goto 1111
      elseif(nbran.ne.nbranv) then
        write(*,'(/,a,2i5)')
     &  ' nbran differs between files, nbran=',nbran,nbranv
        lref=lref+2
        lvar=lvar+2
        npts=0
        goto 1111
      endif
	  
c  check that sample intervals are consistent
      if(dt.ne.sngl(ddt)) then
        dt=sngl(ddt)
        write(*,103) dt
  103   format(/,' sample interval changed to ',f5.2)
      endif
c  check that component and input files are consistent
      if((jcom.ne.1)) then
       write(*,'(a)') ' ***wrong file type***'
       npts=0
       lref=lref+2
       lvar=lvar+2
       goto 1111
      endif
	  
c  choose source depth
      if(vary.ne.'m') then
        do  i=1,nsrce
        read(lrec) sddep(i)
        read(lsou) sddep(i)
        sdep(i)=sngl(sddep(i))
        enddo
      else
        do  i=1,nsrce
        read(lvar) sddep(i)
        read(lref) sddep(i)
        sdep(i)=sngl(sddep(i))
        enddo
      endif
c  find nearest source depth in file to event depth
      call srcget(d0,nsdp)
      sourcd=sdep(nsdp)
c initialize arrays and constants

      nper=npts/2
      nh1=nper+1
      npts2=npts+2
      trec=npts*dt
      df=2.*pi/trec
      nlen=npts
 
      do j=1,9
      do 5 i=1,npts2
    5 strn(i,j)=0.
      enddo
	  
c calculate instrument response and set number of modes
      call inresp(npts)
      call disply

c   fact accounts for distance, sign conventions, and fft normalization
c   fft sign conventions are as in mendiguren
      fact=-fmom*sqrt(1.e3/dist)*2./trec

c  read in vertical and radial excitation functions
       do 8888 kk=1,nbran
       read(lvar,end=10) nb
       read(lref,end=10) nbr
       if(nb.ne.nbr) then
         write(*,'(/,a,2i5)') 
     &  ' incompatible mode number in files! ',nb,nbr
         npts=0
         lref=lref+2
         lvar=lvar+2
         goto 1111
       endif
        do 888 j=1,nh1
c read in rayleigh excitations
        iki=2*(nper+2-j)
	ikr=iki-1
        if(vary.eq.'s') then
          read(lsou,err=223) w,ccc,gmmm,b1,b2,gb1,gb2
          ccv=sngl(ccc)
          gmv=sngl(gmmm)
          do k=1,nsrce
          read(lsou) py1,py2,py3
          if(k.eq.nsdp)then
            py(1)=sngl(py1)
            py(2)=sngl(py2)
            py(3)=sngl(py3)
          endif
          enddo
 223      read(lrec,err=222) w,ccc,gmmm,b1,b2,gb1,gb2
          ccr=sngl(ccc)
          gmr=sngl(gmmm)
          do k=1,nsrce
          read(lrec) py1,py2,py3
          enddo
        elseif(vary.eq.'r') then
          read(lsou,err=323) w,ccc,gmmm,b1,b2,gb1,gb2
          ccr=sngl(ccc)
          gmr=sngl(gmmm)
          do k=1,nsrce
          read(lsou) py1,py2,py3
          if(k.eq.nsdp)then
            py(1)=sngl(py1)
            py(2)=sngl(py2)
            py(3)=sngl(py3)
          endif
          enddo
 323      read(lrec,err=222) w,ccc,gmmm,b1,b2,gb1,gb2
          ccv=sngl(ccc)
          gmv=sngl(gmmm)
          do k=1,nsrce
          read(lrec) py1,py2,py3
          enddo
       else
          read(lvar,err=224) w,ccc,gmmm,b1,b2,gb1,gb2
          ccv=sngl(ccc)
          gmv=sngl(gmmm)
          do k=1,nsrce
          read(lvar) py1,py2,py3
          enddo
224       read(lref,err=222) w,ccc,gmmm,b1,b2,gb1,gb2
          ccr=sngl(ccc)
          gmr=sngl(gmmm)
          do k=1,nsrce
          read(lref) py1,py2,py3
          if(k.eq.nsdp)then
            py(1)=sngl(py1)
            py(2)=sngl(py2)
            py(3)=sngl(py3)
          endif
          enddo
       endif

       if(nb.ge.mdmin.and.nb.le.mdmx)then
	xkv=distv/ccv
	xkr=distr/ccr
        wvd=(xkv+xkr)*w
        cwd=cos(wvd)
        swd=sin(wvd)
		
	gr=cmplx( py(1)*(s2*sol(6)+c2m1*sol(3)+
     & c2p1*sol(2))+py(3)*sol(1),py(2)*(c1*sol(4)+s1*sol(5)))
	 
        atn=exp(-(gmr*distr+gmv*distv))
		
	u=gr*cmplx(atn*b2,0.)*cmplx(swd,cwd)*1.e-5
	wu=gr*cmplx(atn*b1,0.)*cmplx(cwd,-swd)*1.e-5
		
	dudr=-cmplx(0.,xk)*u
	dudz=-gb2*u/b2
		
	dwdr=-cmplx(0.,xk)*wu
	dwdz=-gb1*wu/b1
		
	erz=.5*(dwdr+dudz)
		
        if(ivel.gt.0) then
          strn(ikr,7)=strn(ikr,7)+ dimag(wu*w)
          strn(iki,7)=strn(iki,7)- dreal(wu*w)
          strn(ikr,8)=strn(ikr,8)+ dimag(u*w)
          strn(iki,8)=strn(iki,8)- dreal(u*w)
        else
          strn(ikr,7)=strn(ikr,7)+ real(wu)*1.e5
          strn(iki,7)=strn(iki,7)+ aimag(wu)*1.e5
          strn(ikr,8)=strn(ikr,8)+ real(u)*1.e5
          strn(iki,8)=strn(iki,8)+ aimag(u)*1.e5
        endif
		
        strn(ikr,1)=strn(ikr,1)+ real(dudr)
        strn(iki,1)=strn(iki,1)+ aimag(dudr)

        strn(ikr,3)=strn(ikr,3)+ real(dwdz)
        strn(iki,3)=strn(iki,3)+ aimag(dwdz)

        strn(ikr,5)=strn(ikr,5)+ real(erz)
        strn(iki,5)=strn(iki,5)+ aimag(erz)

       endif
  888  continue
  222  write(*,333) nb
       if(nb.eq.mdmx) go to 11
 8888  continue

c  read in love excitations
11    lref=lref+2
      lvar=lvar+2
      lsou=lsou+2
      lrec=lrec+2

      read(lref) nsrce,npts,ddt,jcom,nbran
      read(lvar) nsrcv,nptv,ddtv,jcomv,nbranv
      if(nsrce.ne.nsrcv) then
        write(*,'(/,a,2i5)')
     &  ' nsrce differs between files, nsrce=',nsrce,nsrcv
        npts=0
        goto 1111
      elseif(nptv.ne.npts) then
        write(*,'(/,a,2i5)')
     &  ' npts differs between files, npts=',npts,nptv
        npts=0
        goto 1111
      elseif(ddt.ne.ddtv) then
        write(*,'(/,a,f8.4,1x,f8.4)')
     &  ' dt differs between files, dt=',ddt,ddtv
        npts=0
        goto 1111
      elseif(jcom.ne.jcomv) then
        write(*,'(/,a,2i5)')
     &  ' jcom differs between files, jcom=',jcom,jcomv
        npts=0
        goto 1111
      elseif(nbran.ne.nbranv) then
        write(*,'(/,a,2i5)')
     &  ' nbran differs between files, nbran=',nbran,nbranv
        npts=0
        goto 1111
      endif

c  check that sample intervals and npts are consistent
      if(dt.ne.sngl(ddt)) then
	write(*,'(a)') ' love file dt is inconsistent!'
	npts=0
        goto 1111
      endif
      if(npts.ne.nlen) then
	write(*,'(a)') ' love file npts is inconsistent!'
	npts=0
        goto 1111
      endif
c  check that component and input files are consistent
      if((jcom.ne.2)) then
        write(*,'(a)') ' ***wrong file type***'
        npts=0
        goto 1111
      endif

c  choose source depth
      if(vary.ne.'m') then
        do  i=1,nsrce
        read(lrec) sddep(i)
        read(lsou) sddep(i)
        sdep(i)=sngl(sddep(i))
        enddo
      else
        do  i=1,nsrce
        read(lvar) sddep(i)
        read(lref) sddep(i)
        sdep(i)=sngl(sddep(i))
        enddo
      endif
c  find nearest source depth in file to event depth
      call srcget(d0,nsdp)
c     call srcget
      if(sourcd.ne.sdep(nsdp)) then
        write(*,'(a)') ' love file source depths inconsistent!'
        npts=0
        goto 1111
      endif

      do 7777 kk=1,nbran
      read(lvar,end=10) nb
      read(lref,end=10) nbr
      if(nb.ne.nbr) then
        write(*,'(/,a,2i5)') 
     &  ' incompatible mode number in files! ',nb,nbr
        goto 1111
      endif
        do 777 j=1,nh1
        iki=2*(nper+2-j)
	ikr=iki-1
c  read in love excitations
        if(vary.eq.'s') then
          read(lsou,err=556) w,ccc,gmmm,b1,gb1
          ccv=sngl(ccc)
          gmv=sngl(gmmm)
          do k=1,nsrce
          read(lsou) py1,py2
          if(k.eq.nsdp)then
            py(1)=sngl(py1)
            py(2)=sngl(py2)
          endif
          enddo
 556      read(lrec,err=555) w,ccc,gmmm,b1,gb1
          ccr=sngl(ccc)
          gmr=sngl(gmmm)
          do k=1,nsrce
          read(lrec) py1,py2
          enddo
        elseif(vary.eq.'r') then
          read(lsou,err=756) w,ccc,gmmm,b1,gb1
          ccr=sngl(ccc)
          gmr=sngl(gmmm)
          do k=1,nsrce
          read(lsou) py1,py2
          if(k.eq.nsdp)then
            py(1)=sngl(py1)
            py(2)=sngl(py2)
          endif
          enddo
 756      read(lrec,err=555) w,ccc,gmmm,b1,gb1
          ccv=sngl(ccc)
          gmv=sngl(gmmm)
          do k=1,nsrce
          read(lrec) py1,py2
          enddo
       else
          read(lvar,err=557) w,ccc,gmmm,b1,gb1
          ccv=sngl(ccc)
          gmv=sngl(gmmm)
          do k=1,nsrce
          read(lvar) py1,py2
          enddo
557       read(lref,err=555) w,ccc,gmmm,b1,gb1
          ccr=sngl(ccc)
          gmr=sngl(gmmm)
          do k=1,nsrce
          read(lref) py1,py2
          if(k.eq.nsdp)then
            py(1)=sngl(py1)
            py(2)=sngl(py2)
          endif
          enddo
       endif

       if(nb.ge.mdmin.and.nb.le.mdmx)then

	xkv=distv/ccv
	xkr=distr/ccr
        wvd=(xkv+xkr)*w
        cwd=cos(wvd)
        swd=sin(wvd)
		
	gl=cmplx( py(2)*(-c1*sol(5)+s1*sol(4)),
     &      py(1)*(c2*sol(6)+s22*(sol(3)-sol(2))) )
	 		
        atn=exp(-(gmr*distr+gmv*distv))
			 		
	v=gl*cmplx(atn*b1,0.)*cmplx(cwd,-swd)*1.e-5
		
	dvdr=-cmplx(0.,xk)*v
	dvdz=0.
	if(b1.ne.0.) dvdz=-.5*gb1*v/b1

        if(ivel.gt.0) then
          strn(ikr,9)=strn(ikr,9)+ dimag(v*w)
          strn(iki,9)=strn(iki,9)- dreal(v*w)
        else
          strn(ikr,9)=strn(ikr,9)+ real(v)*1.e5
          strn(iki,9)=strn(iki,9)+ aimag(v)*1.e5
        endif
						
        strn(ikr,2)=0.
	strn(iki,2)=0.

        strn(ikr,4)=strn(ikr,4)+ real(dvdz)
        strn(iki,4)=strn(iki,4)+ aimag(dvdz)

        erp=.5*dvdr
        strn(ikr,6)=strn(ikr,6)+ real(erp)
        strn(iki,6)=strn(iki,6)+ aimag(erp)

       endif
  777  continue
  555  write(*,333) nb
       if(nb.eq.mdmx) go to 10
 7777  continue
  10   continue
333    format(' mode ',i2,' completed')  
  
c   include travel time,pi/4 phase shifts and instr. response

       do 666 kk=2,nh1
       iki=2*kk
       ikr=iki-1
       w=float(kk-1)*df
       pht=pi*.25+w*to
       cph=cos(pht)
       sph=sin(pht)
       respr=real(resp(kk))
       respi=aimag(resp(kk))
       crl=(respr*cph-respi*sph)*fact
       cri=(respr*sph+respi*cph)*fact
       do jj=1,9
       save=(strn(ikr,jj)*crl-strn(iki,jj)*cri)
       strn(iki,jj)=(strn(ikr,jj)*cri+strn(iki,jj)*crl)
       strn(ikr,jj)=save
       enddo

  666  continue
  
c fft to the time domain
      do j=1,9
      call fftl(strn(1,j),nlen,-2,ier)
      if(ier.ne.0) then
        write(*,'(a)') ' fft error!'
	npts=0
        goto 1111
      endif
      enddo  

1111  close(lref)
      close(lvar)
      lref=lref-2
      lvar=lvar-2
      close(lref)
      close(lvar)
      return
      end
