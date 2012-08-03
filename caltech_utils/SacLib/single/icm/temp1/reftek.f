      subroutine reftek ( nfreq, delfrq, xre, xim, freepd, damp,
     *                  crfrq, hpfrq )
c
c
c   .....REFTEK - for seismometer-REFTEK box response.....
c             (poles and zeros due to G. Randall, modified version
c              port routine original by H. Patton)
c
c MODIFICATIONS:
c    901012:  Deleted the 1.0e-07 scale factor from the displacement
c             response constant. (wct)
c    900409:  Deleted the gain argument and coding.
c
      real * 8 delfrq, xre(1), xim(1), twopi
      complex zero(4), pole(9)
      real * 4 freepd, damp, crfrq, hpfrq, temp
c
      data twopi / 6.283185307179586 /
c
c     fp = twopi / freepd       /*   (Radians)
      fp = twopi / freepd
      dc = damp
c     cf = twopi * crfrq        /*   Corner Frequency ( Radians )
      cf = twopi * crfrq
c
c
c   .....Six pole Butterworth filter.....
c        Upper left quadrant poles, lower left by symmetry
c
      sfil1r = -cf * sin( twopi/24 )
      sfil1i =  cf * cos( twopi/24 )
      sfil2r = -cf * sin( twopi*3/24 )
      sfil2i =  cf * cos( twopi*3/24 )
      sfil3r = -sfil1i
      sfil3i = -sfil1r
c
c   Seismometer
c    Allow for variable seismometer free period and damping constant.
c    Numerator: S*S ( for velocitygram ).
c    Poles are roots of denominator:
c                 2                 2
c              (S) + 2DC(FP)S + (FP) = 0.0
c         Where:
c                 S = i * omega
c                DC = damping constant.
c                FP = 2 Pi * ( free period in sec. )    /*   (radians)
c
      discrm = dc * dc - 1.0
      if ( discrm .gt. 0.0 ) then
c
c              .....Overdamped.....
c
	 discrm = sqrt ( discrm )
	 s1r = fp * ( -dc + discrm )
	 s1i = 0.0
	 s2r = fp * ( -dc - discrm )
	 s2i = 0.0
c
      else
	 discrm = abs ( discrm )
	 discrm = sqrt ( discrm )
	 s1r = -dc * fp
	 s1i = fp * discrm
	 s2r = s1r
	 s2i = -s1i
      end if
c
c   .....computing displacement response.....
c	 if the hpfrq > 0., then add a single zero at D.C. and
c	 a single pole at the high pass filter freq
c
      const = cf * cf * cf * cf * cf * cf
      nzero = 3
      zero(1) = cmplx ( 0.0, 0.0 )
      zero(2) = cmplx ( 0.0, 0.0 )
      zero(3) = cmplx ( 0.0, 0.0 )
c
      npole = 8
c     Butterworth filter, 6 poles
      pole(1) = cmplx ( sfil1r, sfil1i )
      pole(2) = cmplx ( sfil2r, sfil2i )
      pole(3) = cmplx ( sfil3r, sfil3i )
      pole(4) = cmplx ( sfil1r, -sfil1i )
      pole(5) = cmplx ( sfil2r, -sfil2i )
      pole(6) = cmplx ( sfil3r, -sfil3i )
c     Seismometer poles
      pole(7) = cmplx ( s1r, s1i )
      pole(8) = cmplx ( s2r, s2i )
c     The high pass filter if present.....
      if ( hpfrq .gt. 0.0 ) then
	 nzero = 4
	 npole = 9
	 zero(4) = cmplx( 0.0, 0.0 )
         temp = -twopi * hpfrq
	 pole(9) = cmplx( temp, 0.0 )
       endif
c
      call getran( nfreq, delfrq, const, nzero, zero,
     *             npole, pole, xre, xim )
c
      return
      end
