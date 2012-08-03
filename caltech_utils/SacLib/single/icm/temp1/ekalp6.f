      subroutine ekalp6 ( nfreq, delfrq, xre, xim )
c
      include '../../inc/mach'
c
c   .....EKA LP6 - Blacknest specified poles and zeros.....
c
      real * 8 delfrq, xre(1), xim(1)
      complex zero(6), pole(9)
c
c
c   .....Set poles and zeros.....
c
c                             /*  AMP of 1.0 at FREQ of 0.05 Hz
      const = 0.1084564
      nzero = 6
         do 20 i = 1, 4
	 zero(i) = cmplx ( 0.0, 0.0 )
   20    continue
      zero(5) = cmplx ( -0.00524, 1.04720 )
      zero(6) = cmplx ( -0.00524, -1.04720 )
c
      npole = 9
      pole(1) = cmplx ( -0.29323, 0.29915 )
      pole(2) = cmplx ( -0.29323, -0.29915 )
      pole(3) = cmplx ( -0.10996, 0.11218 )
      pole(4) = cmplx ( -0.10996, -0.11218 )
      pole(5) = cmplx ( -0.03140, 0.0 )
      pole(6) = cmplx ( -0.22000, 0.22400 )
      pole(7) = cmplx ( -0.22000, -0.22400 )
      pole(8) = cmplx ( -1.04720, 0.0 )
      pole(9) = cmplx ( -1.04720, 0.0 )
c
c   .....Compute transfer function.....
c
      call getran ( nfreq, delfrq, const, nzero, zero, npole,
     .              pole, xre, xim )
c
      return
      end
