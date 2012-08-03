      subroutine eyeomg ( nfreq, delfrq, xre, xim, nzer )
c
      include '../../inc/mach'
c
      complex zero(1), pole(30)
      real * 8 delfrq, xre(1), xim(1)
c
c   .....I - Omega.....
c
c
c   .....Set poles and zeros.....
c
      const = 1.0
      npole = 0
c
	 do 30 i = 1, nzer
	 zero(i) = cmplx(0.0,0.0)
   30    continue
c
c   .....Compute transfer function.....
c
      call getran ( nfreq, delfrq, const, nzer, zero, npole,
     .              pole, xre, xim )
c
      return
      end
