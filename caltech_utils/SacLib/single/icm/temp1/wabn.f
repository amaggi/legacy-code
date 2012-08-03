      subroutine wabn ( nfreq, delfrq, xre, xim )
c
      include '../../inc/mach'
c
c   .....WABN - WOOD-ANDERSON: Blacknest specified poles and zeros.....
c
      real * 8 delfrq, xre(1), xim(1)
      complex zero(2), pole(2)
      real * 4 const
c
c
c   .....Set poles and zeros.....
c
c     const = 2.077594       /*    AMP of 1.0 at FREQ of 1.0 Hz
      const = 2.077594
      nzero = 2
      zero(1) = cmplx ( 0.0, 0.0 )
      zero(2) = cmplx ( 0.0, 0.0 )
c
      npole = 2
      pole(1) = cmplx ( -6.28318, 4.71239 )
      pole(2) = cmplx ( -6.28318, -4.71239 )
c
c   .....Compute transfer function.....
c
      call getran ( nfreq, delfrq, const, nzero, zero, npole,
     .              pole, xre, xim )
c
      return
      end
