      subroutine ptbllp ( nfreq, delfrq, xre, xim )
c
      include '../../inc/mach'
c
c   .....PTBL LP - Blacknest specified poles and zeros.....
c
      real * 8 delfrq, xre(1), xim(1)
      complex zero(5), pole(8)
      real * 4 const
c
c
c   .....Set poles and zeros.....
c
c     const = 2.554443      /*  AMP of 1.0 at FREQ of O.O5 Hz
      const = 2.554443
      nzero = 5
c
	 do 20 i = 1, nzero
	 zero(i) = cmplx ( 0.0, 0.0 )
   20    continue
c
      npole = 8
      pole(1) = cmplx ( -0.2930, 0.2950 )
      pole(2) = cmplx ( -0.2930, -0.2950 )
      pole(3) = cmplx ( -0.0256, 0.0439 )
      pole(4) = cmplx ( -0.0256, -0.0439 )
      pole(5) = cmplx ( -0.0417, 0.0417 )
      pole(6) = cmplx ( -0.0417, -0.0417 )
      pole(7) = cmplx ( -6.2800, 0.0 )
      pole(8) = cmplx ( -0.5700, 0.0 )
c
c   .....Compute transfer function.....
c
      call getran ( nfreq, delfrq, const, nzero, zero, npole,
     .              pole, xre, xim )
c
      return
      end
