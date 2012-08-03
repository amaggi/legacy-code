      subroutine ykasp ( nfreq, delfrq, xre, xim )
c
      include '../../inc/mach'
c
c  .....YKA SP - Blacknest specified poles and zeros.....
c
      real * 8 delfrq, xre(1), xim(1)
      complex zero(4), pole(6)
      real * 4 const
c
c
c   .....Set poles and zeros.....
c
c     const = 3.8605143 e5      /*   Amp of 1.0 at Freq of 1.0 Hz
      const = 3.8605143 e5
c
      nzero = 4
         do 20 i = 1, nzero
	 zero(i) = cmplx ( 0.0, 0.0 )
   20    continue
c
      npole = 6
      pole(1) = cmplx ( -3.830, 4.980 )
      pole(2) = cmplx ( -3.830, -4.980 )
      pole(3) = cmplx ( -88.700, 88.700 )
      pole(4) = cmplx ( -88.700, -88.700 )
      pole(5) = cmplx ( -0.628, 0.0 )
      pole(6) = cmplx ( -125.66, 0.0 )
c
c   .....Compute transfer function.....
c
      call getran ( nfreq, delfrq, const, nzero, zero, npole,
     .              pole, xre, xim )
c
      return
      end
