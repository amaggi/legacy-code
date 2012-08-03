      subroutine wwspbn ( nfreq, delfrq, xre, xim )
c
      include '../../inc/mach'
c
c  .....WWSSN SP - Blacknest specified poles and zeros.....
c
      real * 8 delfrq, xre(1), xim(1)
      complex zero(3), pole(5)
      real * 4 const
c
c
c   .....Set poles and zeros.....
c
      const = 432.83395
c
      nzero = 3
      zero(1) = cmplx ( 0.0, 0.0 )
      zero(2) = cmplx ( 0.0, 0.0 )
      zero(3) = cmplx ( 0.0, 0.0 )
c
      npole = 5
      pole(1) = cmplx ( -4.04094, 6.47935 )
      pole(2) = cmplx ( -4.04094, -6.47935 )
      pole(3) = cmplx ( -9.25238, 0.0 )
      pole(4) = cmplx ( -7.67430, 0.0 )
      pole(5) = cmplx ( -16.72981, 0.0 )
c
c   .....Compute transfer function.....
c
      call getran ( nfreq, delfrq, const, nzero, zero, npole,
     .              pole, xre, xim )
c
      return
      end
