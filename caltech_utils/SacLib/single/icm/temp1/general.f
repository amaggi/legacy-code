      subroutine general ( nfreq, delfrq, xre, xim, nzer, t0,
     .           h, const )
c
      include '../../inc/mach'
c
c   .....GENERAL - for a general seismometer.....
c
      real * 8 delfrq, xre(1), xim(1), twopi
      complex zero(3), pole(2), crad
c
      data twopi / 6.283185307179586 /
c
c
c   .....Set poles and zeros.....
c
      nzero = 3
         do 20 i = 1, nzero
	 zero(i) = cmplx ( 0.0, 0.0 )
   20    continue
      nzero = nzer
c
      npole = 2
      omo = twopi / t0
      crad = csqrt ( cmplx ( 1.0 - h**2, 0.0 ) )
      pole(1) = omo * ( -h + cmplx ( 0.0, 1.0 ) * crad )
      pole(2) = omo * ( -h - cmplx ( 0.0, 1.0 ) * crad )
c
c   .....Compute transfer function.....
c
      call getran ( nfreq, delfrq, const, nzero, zero, npole,
     .              pole, xre, xim )
c
      return
      end
