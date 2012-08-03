      subroutine wwlpbn ( nfreq, delfrq, xre, xim )
c
      include '../../inc/mach'
c
c  .....WWSSN LP - Blacknest specified poles and zeros.....
c
      real * 8 delfrq, xre(1), xim(1)
      complex zero(3), pole(4)
      real * 4 const
c
c
c   .....Set poles and zeros.....
c
c     const = 0.5985275       /*    AMP of 1.0 at FREQ of 0.05 Hz
      const = 0.5985275
c
      nzero = 3
      zero(1) = cmplx ( 0.0, 0.0 )
      zero(2) = cmplx ( 0.0, 0.0 )
      zero(3) = cmplx ( 0.0, 0.0 )
c
      npole = 4
      pole(1) = cmplx ( -0.257, 0.3376 )
      pole(2) = cmplx ( -0.257, -0.3376 )
      pole(3) = cmplx ( -0.06283, 0.0 )
      pole(4) = cmplx ( -0.06283, 0.0 )
c
c   .....Compute transfer function.....
c
      call getran ( nfreq, delfrq, const, nzero, zero, npole,
     .              pole, xre, xim )
c
      return
      end
