      subroutine acc ( nfreq, delfrq, xre, xim )
c
      include '../../inc/mach'
c
      complex zero(1), pole(2)
      real * 8 delfrq, xre(1), xim(1)
c
c   .....Acceleration Spectral Operator.....
c
c
c   .....Set poles and zeros.....
c
      const = 1.0
      npole = 0
      nzero = 2
c
	 do 20 i = 1, nzero
	 zero(i) = cmplx(0.0,0.0)
   20    continue
c
c   .....Compute transfer function.....
c
      call getran ( nfreq, delfrq, const, nzero, zero, npole,
     .              pole, xre, xim )
c
      return
      end
