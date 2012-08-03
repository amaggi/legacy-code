      subroutine wiech ( nfreq, delfrq, xre, xim )
c
      include '../../inc/mach'
c
c   .....WIEC - for a Wiechert seismometer.....
c
      real * 8 delfrq, xre(1), xim(1), twopi
      complex zero(2), pole(2)
      real * 4 const
c
      data twopi / 6.283185307179586 /
c
c
c MECHANICAL INSTRUMENT:  AMPLITUDE AND PHASE RESPONSES ARE COMPUTED
c FROM EQNS. (82), (98), AND (13) OF SOHON (1932), SEISMOMETRY, IN
c PART II OF INTRODUCTION TO THEORETICAL SEISMOLOGY, J. B. MACELWANE
c AND F. W. SOHON, JOHN WILEY AND SONS, NEW YORK.
c
      t0 = 9.65
      h = 0.403712752
      const = 188.5
      om0 = twopi / t0
      nzero = 2
      zero(1) = ( 0.0, 0.0 )
      zero(2) = ( 0.0, 0.0 )
c
      npole = 2
c  ???????????????????????????????????????
      RAD=SQRT(1.0-H**2)
      POLE(1)=CMPLX(-OM0*H,OM0*RAD)
      POLE(2)=CMPLX(-OM0*H,-OM0*RAD)
c  ???????????????????????????????????????
c
c   .....Compute transfer function.....
c
      call getran ( nfreq, delfrq, const, nzero, zero, npole,
     .              pole, xre, xim )
c
      return
      end
