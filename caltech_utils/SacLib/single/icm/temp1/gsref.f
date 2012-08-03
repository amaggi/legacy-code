      subroutine gsref ( nfreq, delfrq, xre, xim )
c
      include '../../inc/mach'
c
      complex zero(5), pole(11)
      real * 8 delfrq, xre(1), xim(1)
c
c   .....USGS Refraction Seismometer.....
c
c
c   .....Set poles and zeros (due to J. Zucca).....
c
      const = 276.46 * 276.46 * 283.177 * 293.349 * 293.349 *
     .        330.873 * 330.873
      nzero = 5
c
	 do 20 i = 1, nzero
	 zero(i) = cmplx(0.0,0.0)
   20    continue
c
      npole = 11
      pole( 1) = - cmplx ( 10.0531, 7.5398 )
      pole( 2) = - cmplx ( 10.0531, -7.5398 )
      pole( 3) = - cmplx ( 0.5969, 0.0 )
      pole( 4) = - cmplx ( 0.5969, 0.0 )
      pole( 5) = - cmplx ( 276.4602, 0.0 )
      pole( 6) = - cmplx ( 276.4602, 0.0 )
      pole( 7) = - cmplx ( 283.1769, 0.0 )
      pole( 8) = - cmplx ( 260.2009, 135.4598 )
      pole( 9) = - cmplx ( 260.2009, -135.4598 )
      pole(10) = - cmplx ( 180.6564, 277.2001 )
      pole(11) = - cmplx ( 180.6564, -277.2001 )
c
c   .....Compute transfer function.....
c
      call getran ( nfreq, delfrq, const, nzero, zero, npole,
     .              pole, xre, xim )
c
      return
      end
