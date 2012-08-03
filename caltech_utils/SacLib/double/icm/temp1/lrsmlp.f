      subroutine lrsmlp ( nfreq, delfrq, xre, xim )
c
      include '../../inc/mach'
c
c   .....LRSM LP - Blacknest specified poles and zeros.....
c
      real * 8 delfrq, xre(1), xim(1)
      complex zero(4), pole(7)
      real * 4 const
c
c
c   .....Set poles and zeros.....
c
      const = 0.65871 e-1
      nzero = 4
c
	 do 20 i = 1, nzero
	 zero(i) = cmplx ( 0.0, 0.0 )
   20    continue
c
      npole = 7
      pole( 1) = cmplx ( -0.19635, 0.24524 )
      pole( 2) = cmplx ( -0.19635, -0.24524 )
      pole( 3) = cmplx ( -0.20942, 0.0 )
      pole( 4) = cmplx ( -0.20942, 0.0 )
      pole( 5) = cmplx ( -0.00628, 0.0 )
      pole( 6) = cmplx ( -0.17593, 0.17948 )
      pole( 7) = cmplx ( -0.17593, -0.17948 )
c
c   .....Compute transfer function.....
c
      call getran ( nfreq, delfrq, const, nzero, zero, npole,
     .              pole, xre, xim )
c
      return
      end
