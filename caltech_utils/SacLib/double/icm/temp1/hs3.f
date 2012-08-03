      subroutine hs3 ( nfreq, delfrq, xre, xim )
c
      include '../../inc/mach'
c
c   .....S-750 Seismometer Response - output is in units of volts/nm....
c
      real * 8 delfrq, xre(1), xim(1), twopi
      complex cs, cht, chhp, cha1, cha2, cf1, cf2, cha, chb, chs3
c
      data twopi / 6.283185307179586 /
c
c
      delomg = twopi * delfrq
c
c   .....Set poles and zeros.....
c
         do 20 i = 1, nfreq
         omega = float (i-1) * delomg
         cs = cmplx ( 0.0, omega )
	 cht = 2.6 e3 / ( 1.11 e5 + 33.3 * cs + cs**2 )
	 chhp = 16.0 * cs / ( 1.0 + 1.0 e-4 * cs )
         cha1 = ( 21.0 + 1.0 e-4 * cs ) / ( 1.0 + 1.0 e-4 * cs )
         cha2 = -2.0 / ( 1.0 + 2.0 e-5 * cs )
         cf1 = ( 1.0 + 1.7 e-4 * cs ) / ( 1.0 +
     .         1.038 e-3 * cs + 1.47 e-7 * cs**2 + 1.7 e-12 * cs**3 )
         cf2 = 3.0 e-3 * cs
         cha = cht * chhp * cha1 * cha2
         chb = cf1 * cf2
         chs3 = 100.0 * cs *
     .              ( cha / ( 1.0 - cha * chb ) ) * 1.0 e-9 * cs
         xre(i) = real ( chs3 )
         xim(i) = aimag ( chs3 )
   20    continue
c
      return
      end
