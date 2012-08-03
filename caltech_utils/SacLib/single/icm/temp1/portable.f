      subroutine portable(nfreq,delfrq,xre,xim,freepd,damp,crfrq)
c
      include '../../inc/mach'
c
c   .....PORTABLE - for portable seismometer-PDR2 response.....
c             (poles and zeros due to H. Patton)
c
      real * 8 delfrq, xre(1), xim(1), twopi
      complex zero(2), pole(4)
      real * 4 freepd, damp, crfrq
c
      data twopi / 6.283185307179586 /
c
c
      sqr202 = 0.707106781
      delomg = twopi * delfrq
c     fp = twopi / freepd       /*   (Radians)
      fp = twopi / freepd
      dc = damp
c     cf = twopi * crfrq        /*   Corner Frequency ( Radians )
      cf = twopi * crfrq
c
c   .....Two pole Butterworth filter.....
c
      comp = sqr202 * cf
c     const = comp * comp       /*   numerator
      const = comp * comp
      sfil1r = -comp
      sfil1i =  comp
      sfil2r = -comp
      sfil2i = -comp
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
c   .....computing velocity response.....
c
      nzero = 2
c     zero(1) = cmplx ( 0.0, 0.0 )              /*  Seismometer.
      zero(1) = cmplx ( 0.0, 0.0 )
      zero(2) = cmplx ( 0.0, 0.0 )
c
      npole = 4
c     pole(1) = cmplx ( sfil1r, sfil1i )        /*  Filter
      pole(1) = cmplx ( sfil1r, sfil1i )
c     pole(2) = cmplx ( sfil2r, sfil2i )        /*  Filter
      pole(2) = cmplx ( sfil2r, sfil2i )
c     pole(3) = cmplx ( s1r, s1i )              /*  Seismometer
      pole(3) = cmplx ( s1r, s1i )
c     pole(4) = cmplx ( s2r, s2i )              /*  Seismometer
      pole(4) = cmplx ( s2r, s2i )
c
      asqrd = 0.0
	 do 50 j = 1, nfreq
	 omega = delomg * float(j-1)
	 trn = 1.0 d0
	 tin = 0.0 d0
c
	    do 30 i = 1, nzero
	    tr = -real ( zero(i) )
	    ti = omega - aimag ( zero(i) )
	    tr0 = trn * tr - tin * ti
	    ti0 = trn * ti + tin * tr
	    trn = tr0
	    tin = ti0
   30       continue
c
	 trd = 1.0 d0
	 tid = 0.0 d0
	    do 40 i = 1, npole
	    tr = -real ( pole(i) )
	    ti = omega - aimag ( pole (i) )
	    tr0 = trd * tr - tid * ti
	    ti0 = trd * ti + tid * tr
	    trd = tr0
	    tid = ti0
   40       continue
	 fac = dble ( const ) / ( trd**2 + tid**2 )
	 xre(j) = fac * ( trn * trd + tin * tid )
	 xim(j) = fac * ( trd * tin - trn * tid )
	 astest = xre(j)**2 + xim(j)**2
	 if ( astest .gt. asqrd ) asqrd = astest
   50    continue
c
c   .....Normalization such that velocity transfer function has
c      maximum equal to unity i. e. displacement transfer function
c      has gain equal to 2*PI*FMAX where FMAX is maximum point on
c      velocity transfer function.....
c
c   .....Also transform velocity response to displacement response
c      - multiply by I * OMEGA.....
c
      anorm = 1.0 / sqrt ( asqrd )
	 do 60 j = 1, nfreq
	 omega = delomg * float ( j - 1 )
	 xxre = -xim(j) * dble( omega * anorm )
	 xxim = xre(j) * dble ( omega * anorm )
	 xre(j) = xxre
	 xim(j) = xxim
   60    continue
c
      return
      end
