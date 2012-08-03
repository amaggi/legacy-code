      subroutine snla3 ( nfreq, delfrq, xre, xim )
c
      include '../../inc/mach'
c
c  Response due to H. Patton
c
c   new system response (1985) of Sandia Network ---
c   sl-250 Geotech Seismometers are set at 20 sec.
c   natural period and critically damped.  Signal is
c   low-pass filtered with an eight-pole Butterworth
c   filter with a 3 db point at 10 Hz.  Sampling rate
c   is 50 samples/sec.
c
      real * 8 delfrq, xre(1), xim(1), twopi, trn, tin, tr, ti,
     .   tr0, ti0, trd, tid, fac, xxre, xxim
      complex zero(2), pole(10)
      real * 4 const
c
      data twopi / 6.283185307179586 /
c
      delomg = twopi * delfrq
      const = 1.0
c
c   .....Eight-pole Butterworth Filter.....
c
      cf3db = twopi * 10.
      prd1 = 0.9808 * cf3db
      prd2 = 0.1951 * cf3db
      prd3 = 0.8315 * cf3db
      prd4 = 0.5556 * cf3db
      pole(1) = cmplx ( -prd1, prd2 )
      pole(2) = cmplx ( -prd1, -prd2 )
      pole(3) = cmplx ( -prd3, prd4 )
      pole(4) = cmplx ( -prd3, -prd4 )
      pole(5) = cmplx ( -prd2, prd1 )
      pole(6) = cmplx ( -prd2, -prd1 )
      pole(7) = cmplx ( -prd4, prd3 )
      pole(8) = cmplx ( -prd4, -prd3 )
c
c   .....Critically damped seismometer with t = 20 sec.....
c
      omega0 = twopi / 20.
      pole(9) = cmplx ( -omega0, 0.0 )
      pole(10) = cmplx ( -omega0, 0.0 )
c
      npole = 10
c
c   .....Computing velocity response.....
c
      nzero = 2
      zero(1) = cmplx ( 0.0, 0.0 )
      zero(2) = cmplx ( 0.0, 0.0 )
c
      asqrd = 0.0
	 do 40 j = 1, nfreq
	 omega = delomg * float ( j - 1 )
	 trn = 1.0 d0
	 tin = 0.0 d0
c
	    do 20 i = 1, nzero
	    tr = -real ( zero(i) )
	    ti = omega - aimag ( zero(i) )
	    tr0 = trn * tr - tin * ti
	    ti0 = trn * ti + tin * tr
	    trn = tr0
	    tin = ti0
   20       continue
c
	 trd = 1.0 d0
	 tid = 0.0 d0
	    do 30 i = 1, npole
	    tr = -real ( pole(i) )
	    ti = omega - aimag ( pole(i) )
	    tr0 = trd * tr - tid * ti
	    ti0 = trd * ti + tid * tr
	    trd = tr0
	    tid = ti0
   30       continue
c
	 fac = dble ( const ) / ( trd**2 + tid**2 )
	 xre(j) = fac * ( trn * trd + tin * tid )
	 xim(j) = fac * ( trd * tin - trn * tid )
	 astest = xre(j)**2 + xim(j)**2
	 if ( astest .gt. asqrd ) asqrd = astest
   40    continue
c
c   .....Normalization such that velocity tranfer function
c      has maximum equal to unity i. e. displacement transfer
c      function has gain equal to 2*pi*fmax where fmax is
c      maximum point on velocity transfer function.....
c
c   .....Also transform velocity response to displacement
c      response -- multiply by i * omega.....
c
      anorm = 1.0 / sqrt ( asqrd )
	 do 50 j = 1, nfreq
	 omega = delomg * float ( j - 1 )
	 xxre = -xim(j) * dble ( omega * anorm )
	 xxim = xre(j) * dble ( omega * anorm )
	 xre(j) = xxre
	 xim(j) = xxim
   50    continue
c
      return
      end
