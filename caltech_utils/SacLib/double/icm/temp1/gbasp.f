      subroutine gbasp ( nfreq, delfrq, xre, xim )
c
      include '../../inc/mach'
c
c   .....GBA SP - Blacknest specified poles and zeros.....
c
      real * 8 delfrq, xre(1), xim(1)
      complex zero(3), pole(4)
c
c
c   .....Set poles and zeros.....
c
c                             /*  AMP of 1.0 at FREQ of 1.6 Hz
      const = 2.5597471 e2
      nzero = 3
         do 20 i = 1, nzero
	 zero(i) = cmplx ( 0.0, 0.0 )
   20    continue
c
      npole = 4
      pole(1) = cmplx ( -4.02, 4.82 )
      pole(2) = cmplx ( -4.02, -4.82 )
      pole(3) = cmplx ( -40.2, 30.2 )
      pole(4) = cmplx ( -40.2, -30.2 )
c
c   .....Compute transfer function.....
c
      call getran ( nfreq, delfrq, const, nzero, zero, npole,
     .              pole, xre, xim )
c
      return
      end
