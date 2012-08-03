      subroutine dss ( nfreq, delfrq, xre, xim )
c
      include '../../inc/mach'
c
c   .....LLL DSS - For LLL digital network.....
c      ( poles and zeros due to P. Rodgers )
c
      real * 8 delfrq, xre(1), xim(1), twopi
      complex zero(3), pole(8), crad
c
      data twopi / 6.283185307179586 /
c
c
c   .....Set poles and zeros.....
c
      omo = twopi / 30.0
      const = 6.152890842 e10
      nzero = 3
	 do 20 i = 1, nzero
	 zero(i) = cmplx ( 0.0, 0.0 )
   20    continue
c
      npole = 8
      crad = csqrt ( cmplx ( 1.0 - 0.707**2, 0.0 ) )
      pole(1) = omo * ( -0.707 + cmplx ( 0.0, 1.0 ) * crad )
      pole(2) = omo * ( -0.707 - cmplx ( 0.0, 1.0 ) * crad )
      pole(3) = cmplx ( -33.80165, 60.35942 )
      pole(4) = cmplx ( -33.80165, -60.35942 )
      pole(5) = cmplx ( -50.43199, 35.28386 )
      pole(6) = cmplx ( -50.43199, -35.28386 )
      pole(7) = cmplx ( -57.07708, 11.65531 )
      pole(8) = cmplx ( -57.07708, -11.65531 )
c
c   .....Compute transfer function.....
c
      call getran ( nfreq, delfrq, const, nzero, zero, npole,
     .              pole, xre, xim )
c
      return
      end
