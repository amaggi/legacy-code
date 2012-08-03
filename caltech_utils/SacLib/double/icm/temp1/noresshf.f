      subroutine noresshf ( nfreq, delfrq, xre, xim )
c
      include '../../inc/mach'
c
c   .....NORESS high frequency element (Durham, 6 March 1986).....
c          output in volts/nm
c          A. Smith, 6 May 1986
c
      real * 8 delfrq, xre(1), xim(1), twopi
      complex cs, chhf, chsi, chif, chpa, chia, chhp, chlp
      complex temp1, temp2, temp3
c
      data twopi / 6.283185307179586 /
c
      delomg = twopi * delfrq
c
         do 20 i = 1, nfreq
         omega = float (i-1) * delomg
         cs = cmplx ( 0.0, omega )
	 chsi = 55.73 * cs**2 / ( 1.0 + 0.2387 * cs + 0.02533 * cs**2 )
	 chif = 0.9174 / ( 1.0 + 8.165 e-4 * cs + 7.339 e-7 * cs**2 )
         chpa = 5.544 * cs / ( ( 1.0 + 1.559 e-3 * cs ) *
     .                      ( 1.0 + 0.1598 * cs ) )
         chia = 0.3186 * cs / ( ( 1.0 + 1.593 * cs ) *
     .                      ( 1.0 + 2.496 e-3  * cs ) )

c Use of temporary variables below due to bug in MASSCOMP compiler (870528)

         temp1 = cs * cs
         temp2 = 0.0124 * temp1
         temp3 = ( ( 1.0 + 0.0500 * cs ) * ( 1.0 + 0.124  * cs ) )
         chhp = temp2 / temp3

         chlp = 1.0 / ( 1.0 + ( 2.894 e-3 * cs ) ** 36 )
c
	 chhf = chsi * chif * chpa * chia * chhp * chlp * 1.0 e-9 * cs
c
         xre(i) = real ( chhf )
         xim(i) = aimag ( chhf )
   20    continue
      return
      end
