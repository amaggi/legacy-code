      subroutine lrsmsp ( nfreq, delfrq, xre, xim )
c
      include '../../inc/mach'
c
c   .....LRSM SP - Blacknest specified poles and zeros.....
c
      real * 8 delfrq, xre(1), xim(1)
      complex zero(3), pole(5)
      real * 4 const
c
c
c   .....Set poles and zeros.....
c
c     const = 3.9927709 e3   /*  AMP of 1.0 at FREQ of 1.0 Hz
      const = 3.9927709 e3
      nzero = 3
c
	 do 20 i = 1, nzero
	 zero(i) = cmplx ( 0.0, 0.0 )
   20    continue
c
      npole = 5
      pole( 1) = cmplx ( -5.60893, 7.40733 )
      pole( 2) = cmplx ( -5.60893, -7.40733 )
      pole( 3) = cmplx ( -9.27502, 0.0 )
      pole( 4) = cmplx ( -28.24515, 14.97041 )
      pole( 5) = cmplx ( -28.24515, -14.97041 )
c
c   .....Compute transfer function.....
c
      call getran ( nfreq, delfrq, const, nzero, zero, npole,
     .              pole, xre, xim )
c
      return
      end
