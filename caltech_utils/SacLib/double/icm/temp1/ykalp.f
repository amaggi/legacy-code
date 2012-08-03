      subroutine ykalp ( nfreq, delfrq, xre, xim )
c
      include '../../inc/mach'
c
c  .....YKA LP - Blacknest specified poles and zeros.....
c
      real * 8 delfrq, xre(1), xim(1)
      complex zero(4), pole(6)
      real * 4 const
c
c
c   .....Set poles and zeros.....
c
c     const = 0.28761294212      /*  AMP of 1.0 at FREQ of 0.05 Hz.
      const = 0.28761294212
c
      nzero = 4
         do 20 i = 1, nzero
         zero(i) = cmplx ( 0.0, 0.0 )
   20    continue
c
      npole = 6
      pole(1) = cmplx ( -0.2010, 0.2415 )
      pole(2) = cmplx ( -0.2010, -0.2415 )
      pole(3) = cmplx ( -0.134, 0.161 )
      pole(4) = cmplx ( -0.134, -0.161 )
      pole(5) = cmplx ( -0.628, 0.0 )
c     pole(6) = cmplx ( -0.0134, 0.0 )              /*  KKN
      pole(6) = cmplx ( -0.0134, 0.0 )
c
c   .....Compute transfer function.....
c
      call getran ( nfreq, delfrq, const, nzero, zero, npole,
     .              pole, xre, xim )
c
      return
      end
