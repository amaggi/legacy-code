      subroutine bbvel ( nfreq, delfrq, xre, xim )
c
      include '../../inc/mach'
c
c   .....BB VEL - Blacknest specified poles and zeros.....
c
      real * 8 delfrq, xre(1), xim(1)
      complex zero(4), pole(7)
c
c
c   .....Set poles and zeros.....
c
c                             /*  AMP of 1.0 at FREQ of 1.0 Hz
      const = 5.873691 e4
      nzero = 4
	 do 20 i = 1, nzero
	 zero(i) = cmplx ( 0.0, 0.0 )
   20    continue
c
      npole = 7
      pole(1) = cmplx ( -0.2140, 0.2300 )
      pole(2) = cmplx ( -0.2140, -0.2300 )
      pole(3) = cmplx ( -6.4400, 23.4000 )
      pole(4) = cmplx ( -6.4400, -23.4000 )
      pole(5) = cmplx ( -25.0000, 0.0 )
      pole(6) = cmplx ( -25.0000, 0.0 )
      pole(7) = cmplx ( -0.0555, 0.0 )
c
c   .....Compute transfer function.....
c
      call getran ( nfreq, delfrq, const, nzero, zero, npole,
     .              pole, xre, xim )
c
      return
      end
