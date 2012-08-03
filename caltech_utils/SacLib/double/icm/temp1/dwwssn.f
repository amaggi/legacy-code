      subroutine dwwssn ( nfreq, delfrq, xre, xim )
c
      include '../../inc/mach'
c
c   .....DWWSSN - For a digital WWSSN system response.....
c            ( poles and zeros due to H. Patton )
c
c     digital WWSSN instrumental response.....
c     default response computed assuming poles and zeros for
c     stations LON and JAS supplied on day tapes for day 
c     301, 1983
c
      real * 8 delfrq, xre(1), xim(1)
      complex zero(5), pole(11)
c
c
c   .....Set poles and zeros.....
c
      const = 0.0243 * 500.
      nzero = 5
	 do 20 i = 1, nzero
	 zero(i) = cmplx ( 0.0, 0.0 )
   20    continue
c
      npole = 11
      pole(1) = cmplx ( -0.369, .199 )
      pole(2) = cmplx ( -0.369, -.199 )
      pole(3) = cmplx ( -0.628, 0.0 )
      pole(4) = cmplx ( -0.0209, 0.0 )
      pole(5) = cmplx ( -0.0209, 0.0 )
         do 30 i = 6, npole
         pole(i) = cmplx ( -.273, 0.0 )
   30    continue
c
c   .....Compute transfer function.....
c
      call getran ( nfreq, delfrq, const, nzero, zero, npole,
     .              pole, xre, xim )
c
      return
      end
