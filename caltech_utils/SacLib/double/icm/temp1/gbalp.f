      subroutine gbalp ( nfreq, delfrq, xre, xim )
c
      include '../../inc/mach'
c
c   .....GBA LP - Blacknest specified poles and zeros.....
c
      real * 8 delfrq, xre(1), xim(1)
      complex zero(4), pole(10)
c
c
c   .....Set poles and zeros.....
c
c                                /*  AMP of 1.0 at FREQ of 0.05 Hz
      const = 0.100844452 e-1
      nzero = 4
         do 20 i = 1, nzero
	 zero(i) = cmplx ( 0.0, 0.0 )
   20    continue
c
      npole = 10
      pole(1) = cmplx ( -0.2140, 0.2300 )
      pole(2) = cmplx ( -0.2140, -0.2300 )
      pole(3) = cmplx ( -0.1340, 0.1605 )
      pole(4) = cmplx ( -0.1340, -0.1605 )
      pole(5) = cmplx ( -0.0312, 0.0 )
c                                          /* KKN
      pole(6) = cmplx ( -2.060, 0.0 )
      pole(7) = cmplx ( -0.1670, 0.2550 )
      pole(8) = cmplx ( -0.1670, -0.2550 )
      pole(9) = cmplx ( -0.1670, 0.2550 )
      pole(10) = cmplx ( -0.1670, -0.2550 )
c
c   .....Compute transfer function.....
c
      call getran ( nfreq, delfrq, const, nzero, zero, npole,
     .              pole, xre, xim )
c
      return
      end
