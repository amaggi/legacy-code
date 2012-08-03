      subroutine benbog ( nfreq, delfrq, xre, xim )
c
      include '../../inc/mach'
c
c   .....BEN BOG - Blacknest specified poles and zeros.....
c
      real * 8 delfrq, xre(1), xim(1)
      complex zero(4), pole(8)
c
c
c   .....Set poles and zeros.....
c
c                             /*  AMP of 1.0 at FREQ of 1.0 Hz
      const = 7.5111330 e6
      nzero = 4
	 do 20 i = 1, nzero
	 zero(i) = cmplx ( 0.0, 0.0 )
   20    continue
c
      npole = 8
      pole(1) = cmplx ( -10.539, 0.0 )
      pole(2) = cmplx ( -5.787, 7.407 )
      pole(3) = cmplx ( -5.787, -7.407 )
      pole(4) = cmplx ( -22.21, 22.21 )
      pole(5) = cmplx ( -22.21, -22.21 )
      pole(6) = cmplx ( -0.06283, 0.0 )
      pole(7) = cmplx ( -29.62, 29.62 )
      pole(8) = cmplx ( -29.62, -29.62 )
c
c   .....Compute transfer function.....
c
      call getran ( nfreq, delfrq, const, nzero, zero, npole,
     .              pole, xre, xim )
c
      return
      end
