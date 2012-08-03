      subroutine redkir ( nfreq, delfrq, xre, xim )
c
      include '../../inc/mach'
c
c   .....RED KIRN - Blacknest specified poles and zeros.....
c
      real * 8 delfrq, xre(1), xim(1)
      complex zero(3), pole(4)
      real * 4 const
c
c
c   .....Set poles and zeros.....
c
c     const = 83.692405956 e0  /*  AMP of 1.0 at FREQ of 1.0 Hz
      const = 83.692405956
      nzero = 3
c
	 do 20 i = 1, nzero
	 zero(i) = cmplx ( 0.0, 0.0 )
   20    continue
c
      npole = 4
      pole(1) = cmplx ( -0.12759, 0.23031 )
      pole(2) = cmplx ( -0.12759, -0.23031 )
      pole(3) = cmplx ( -0.29915, 0.0 )
      pole(4) = cmplx ( -83.43929, 0.0 )
c
c   .....Compute transfer function.....
c
      call getran ( nfreq, delfrq, const, nzero, zero, npole,
     .              pole, xre, xim )
c
      return
      end
