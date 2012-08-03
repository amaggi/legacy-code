      subroutine hfslpwb ( nfreq, delfrq, xre, xim )
c
      include '../../inc/mach'
c
c   .....HFS LPWB - Blacknest specified poles and zeros.....
c
      real * 8 delfrq, xre(1), xim(1)
      complex zero(4), pole(6)
c
c
c   .....Set poles and zeros.....
c
c                            /*  AMP of 1.0 at FREQ of 0.05 Hz
      const = 0.1761249
      nzero = 4
         do 20 i = 1, nzero
	 zero(i) = cmplx ( 0.0, 0.0 )
   20    continue
c
      npole = 6
      pole(1) = cmplx ( -0.1100, 0.2380 )
      pole(2) = cmplx ( -0.1100, -0.2380 )
      pole(3) = cmplx ( -0.1340, 0.1605 )
      pole(4) = cmplx ( -0.1340, -0.1605 )
      pole(5) = cmplx ( -0.0312, 0.0 )
      pole(6) = cmplx ( -0.6450, 0.0 )
c
c   .....Compute transfer function.....
c
      call getran ( nfreq, delfrq, const, nzero, zero, npole,
     .              pole, xre, xim )
c
      return
      end
