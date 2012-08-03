      subroutine llsn ( nfreq, delfrq, xre, xim )
c
      include '../../inc/mach'
c
c   .....LLSN - for an L-4 seismometer.....
c     (poles and zeros due to P. Rodgers)
c
      real * 8 delfrq, xre(1), xim(1)
      complex zero(5), pole(10)
      real * 4 const
c
c
c   .....Set poles and zeros.....
c
c     const = -3.57425 e15   /*  = -A * B
      const = -3.57425 e15
c     nzero = 5              /* checked on 2/2/81
      nzero = 5
c
	 do 20 i = 1, nzero
	 zero(i) = cmplx ( 0.0, 0.0 )
   20    continue
c
      npole = 10
      pole( 1) = cmplx ( -5.026548, 3.769911 )
      pole( 2) = cmplx ( -5.026548, -3.769911 )
      pole( 3) = cmplx ( -0.5969026, 0.0 )
      pole( 4) = cmplx ( -0.5969026, 0.0 )
      pole( 5) = cmplx ( -276.460154, 0.0 )
      pole( 6) = cmplx ( -276.460154, 0.0 )
      pole( 7) = cmplx ( -376.99112, 0.0 )
      pole( 8) = cmplx ( -376.99112, 0.0 )
      pole( 9) = cmplx ( -571.76986, 583.32194 )
      pole(10) = cmplx ( -571.76986, -583.32194 )
c
c   .....Compute transfer function.....
c
      call getran ( nfreq, delfrq, const, nzero, zero, npole,
     .              pole, xre, xim )
c
      return
      end
