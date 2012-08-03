      subroutine oldkir ( nfreq, delfrq, xre, xim )
c
      include '../../inc/mach'
c
c   .....OLD KIRN - Blacknest specified poles and zeros.....
c
      real * 8 delfrq, xre(1), xim(1)
      complex zero(3), pole(4)
      real * 4 const
c
c
c   .....Set poles and zeros.....
c
c     const = 8.033295 e1   /*  AMP of 1.0 at FREQ of 1.0 Hz
      const = 8.033295 e1
      nzero = 3
c
	 do 20 i = 1, nzero
	 zero(i) = cmplx ( 0.0, 0.0 )
   20    continue
c
      npole = 4
      pole( 1) = cmplx ( -0.2140, 0.2300 )
      pole( 2) = cmplx ( -0.2140, -0.2300 )
      pole( 3) = cmplx ( -0.3150, 0.0 )
      pole( 4) = cmplx ( -80.0000, 0.0 )
c
c   .....Compute transfer function.....
c
      call getran ( nfreq, delfrq, const, nzero, zero, npole,
     .              pole, xre, xim )
c
      return
      end
