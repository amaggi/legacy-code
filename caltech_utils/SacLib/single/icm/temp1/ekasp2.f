      subroutine ekasp2 ( nfreq, delfrq, xre, xim )
c
      include '../../inc/mach'
c
c   .....EKA SP2 - Blacknest specified poles and zeros.....
c
      real * 8 delfrq, xre(1), xim(1)
      complex zero(4), pole(4)
c
c
c   .....Set poles and zeros.....
c
c                                /*  AMP of 1.0 at FREQ of 1.0 Hz
      const = 12.9083123695
      nzero = 4
         do 20 i = 1, nzero
	 zero(i) = cmplx ( 0.0, 0.0 )
   20    continue
c
      npole = 4
      pole(1) = cmplx ( -3.7700, 5.0200 )
      pole(2) = cmplx ( -3.7700, -5.0200 )
      pole(3) = cmplx ( -67.2000, 0.0 )
      pole(4) = cmplx ( -0.3300, 0.0 )
c
c   .....Compute transfer function.....
c
      call getran ( nfreq, delfrq, const, nzero, zero, npole,
     .              pole, xre, xim )
c
      return
      end
