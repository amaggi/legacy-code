      subroutine lll ( nfreq, delfrq, xre, xim, subtyp, freepd, damp,
     .                 nerr )
c
      include '../../inc/mach'
c
c   .....LLL - for all LLL broadband (analog) seismometers.....
c
      real * 8 delfrq, xre(1), xim(1), twopi
      complex zero(3), pole(8)
      character * 8 subtyp
c
      data twopi / 6.283185307179586 /
c
c
c
c LLL STATION NETWORK AT ELKO, KANAB, LANDERS, MINA.
c THIS SECTION CONSISTS OF THE INSTRUMENT RESPONSE AS WELL AS THE IW
c TERM TO CONVERT DISPLACEMENT INTO VELOCITY.
c THE TRANSFER FUNCTION IS SHIFTED BY PI FROM THE PUBLISHED TRANSFER
c FUNCTION IN ORDER TO CONFORM WITH SIGN CONVENTIONS.
c REFERENCES ARE:
c     JUDD, D. D. (1969), DETERMINATION OF TRANSFER FUNCTION OF
c        SEISMOMETER/TELEMETRY SYSTEM THROUGH USE OF PSEUDO-RANDOM
c        BINARY SEQUENCES, UCRL-50490.
c     DENNY, M. D. (1977), THE INSTALLATION OF HORIZONTAL SEISMOMETERS
c        IN THE LLL SEISMIC NET AND THEIR CALIBRATION, UCRL-52216.
c
      if ( subtyp .eq. 'LV      ' ) then
         t0 = 20.9
         h = 0.800
      else if (subtyp .eq. 'LR      ' ) then
         t0 = 39.2
         h = 0.707
      else if ( subtyp .eq. 'LT      ' ) then
         t0 = 39.6
         h = 0.707
      else if ( subtyp .eq. 'MV      ' ) then
         t0 = 19.3
         h = 0.790
      else if ( subtyp .eq. 'MR      ' ) then
         t0 = 41.5
         h = 0.707
      else if ( subtyp .eq. 'MT      ' ) then
         t0 = 40.6
         h = 0.707
      else if ( subtyp .eq. 'KV      ' ) then
         t0 = 20.0
         h = .800
      else if ( subtyp .eq. 'KR      ' ) then
         t0 = 40.56
         h = 0.707
      else if ( subtyp .eq. 'KT      ' ) then
         t0 = 39.98
         h = 0.707
      else if ( subtyp .eq. 'EV      ' ) then
         t0 = 30.9
         h = 0.770
      else if ( subtyp .eq. 'ER      ' ) then
         t0 = 40.4
         h = 0.707
      else if ( subtyp .eq. 'ET      ' ) then
         t0 = 40.4
         h = 0.707
      else if ( subtyp .eq. 'BB      ' ) then
         t0 = freepd
         h = damp
      else
         nerr=2105
         call setmsg('ERROR',nerr)
         call apcmsg('LLL:')
         call apcmsg(subtyp)
         go to 8888
      end if
c
c
      om0 = twopi / t0
      const = 3.93785011 e12
      nzero = 3
         do 40 i = 1, nzero
         zero(1) = cmplx ( 0.0, 0.0 )
   40    continue
      npole = 8
c    ???????????????????????????????
      RAD=SQRT(1.0-H**2)
      POLE(1)=CMPLX(-OM0*H,OM0*RAD)
      POLE(2)=CMPLX(-OM0*H,-OM0*RAD)
c ????????????????????????????????????
      POLE(3)=CMPLX(-114.28,23.317)
      POLE(4)=CMPLX(-114.28,-23.317)
      POLE(5)=CMPLX(-100.48,70.65)
      POLE(6)=CMPLX(-100.48,-70.65)
      POLE(7)=CMPLX(-67.677,120.85)
      POLE(8)=CMPLX(-67.677,-120.85)
c
c   .....Compute transfer function.....
c
      call getran ( nfreq, delfrq, const, nzero, zero, npole,
     .              pole, xre, xim )
c
 8888 return
      end
