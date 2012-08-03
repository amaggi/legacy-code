      subroutine wwsp ( nfreq, delfrq, xre, xim )
c
      include '../../inc/mach'
c
c  WWSSN short period seismometer.
c  Ref:  Luh, P. C. (1977).  A scheme for expressing instrumental
c  responses parametrically, BSSA, 67, 957-969.
c
      real * 8 delfrq, xre(1), xim(1)
      complex zero(3), pole(5)
      real * 4 const
c
c
c   .....Set poles and zeros.....
c
      const = 397.54767
      nzero = 3
      zero(1) = cmplx ( 0.0, 0.0 )
      zero(2) = cmplx ( 0.0, 0.0 )
      zero(3) = cmplx ( 0.0, 0.0 )
c
      npole = 5
      pole(1) = cmplx ( -5.0136607, 6.4615109 )
      pole(2) = cmplx ( -5.0136607, -6.4615109 )
      pole(3) = cmplx ( -8.2981509, 0.0 )
      pole(4) = cmplx ( -8.6940765, -7.1968661 )
      pole(5) = cmplx ( -8.6940765, 7.1968661 )
c
c   .....Compute transfer function.....
c
      call getran ( nfreq, delfrq, const, nzero, zero, npole,
     .              pole, xre, xim )
c
      return
      end
