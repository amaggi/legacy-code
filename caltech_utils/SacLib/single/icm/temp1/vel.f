      subroutine vel ( nfreq, delfrq, xre, xim )
c
      include '../../inc/mach'
c
c   .....VEL - velocity spectral operator.....
c
      real * 8 delfrq, xre(1), xim(1)
      complex zero(1), pole(1)
      real * 4 const
c
c
c   .....Set poles and zeros.....
c
      const = 1.0
      nzero = 1
c
      zero(1) = cmplx ( 0.0, 0.0 )
c
      npole = 0
c
c   .....Compute transfer function.....
c
      call getran ( nfreq, delfrq, const, nzero, zero, npole,
     .              pole, xre, xim )
c
      return
      end
