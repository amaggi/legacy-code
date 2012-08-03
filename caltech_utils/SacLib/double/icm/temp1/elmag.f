      subroutine elmag ( nfreq, delfrq, xre, xim, freepd, mag, nerr )

* for a WWSSN Elecromagnetic Instrument
c
c
      include '../../inc/mach'
c
      real * 8 delfrq, twopi, xre(1), xim(1)
      complex cs, ct, ctd
      integer * 4 nfreq, nmag
      real mag
c
      data twopi / 6.283185307179586 /
c
c ELECTROMAGNETIC INSTRUMENT:  AMPLITUDE AND PHASE RESPONSES ARE
c COMPUTED FROM EQN. (28) AND (29) --- SEE COMMENT PRECEDING STATEMENT
c NUMBER D460 --- OF HAGIWARA (1958) , BULL. EARTHQ. RES. INST., TOKYO
c UNIV., VOL. 36, P.139-164.
c
c READ T1, T2, H1, H2, AND SIGSQ.
c   T1 IS THE PENDULUM PERIOD, T2 IS THE GALVANOMETER PERIOD.
c   H1 IS THE PENDULUM DAMPING FACTOR AND H2 IS THE GALVANOMETER
c    DAMPING FACTOR.
c   SIGSQ IS THE INSTRUMENT COUPLING FACTOR AS DEFINED BY HAGIWARA.
c   ACCORDING TO CHANDRA (1970), BSSA, 60, 539-563,
c       FOR T1=15, T2=100, H1=0.93, H2=1.0
c
c         PEAK MAGNIFICATION      SIGSQ
c
c                 375             0.003
c                 750             0.013
c                1500             0.047
c                3000             0.204
c                6000             0.805
c
c       FOR T1=30, T2=100, H1=1.50, H2=1.0
c
c         PEAK MAGNIFICATION      SIGSQ
c
c                 375             0.003
c                 750             0.012
c                1500             0.044
c                3000             0.195
c                6000             0.767
c
      nerr = 0
      delomg = twopi * delfrq
c
      nmag=ifix(mag+0.01)
      if ( freepd .eq. 15.0 ) then
c
	 t1 = 15.0
	 t2 = 100.0
         h1 = .93
         h2 = 1.0
         if ( nmag .eq. 375 ) then
           const = 712.48076
           sigsq = 0.003
         else if ( nmag .eq. 750 ) then
           const = 1423.9605
           sigsq = 0.013
         else if ( nmag .eq. 1500 ) then
           const = 2840.9091
           sigsq = 0.047
         else if ( nmag .eq. 3000 ) then
           const = 2708.3153
           sigsq = 0.204
         else if ( nmag .eq. 6000 ) then
           const = 10114.803
           sigsq = 0.805
         else
           nerr=2112
           call setmsg('ERROR',nerr)
           go to 8888
         end if
c
      elseif ( freepd .eq. 30.0 ) then
c
         t1 = 30.0
         t2 = 100.0
         h1 = 1.5
         h2 = 1.0
	 if ( nmag .eq. 375 ) then
           const = 1202.5398
           sigsq = 0.003
         else if ( nmag .eq. 750 ) then
           const = 2402.1523
           sigsq = 0.012
         else if ( nmag .eq. 1500 ) then
           const = 4781.9434
           sigsq = 0.044
         else if ( nmag .eq. 3000 ) then
           const = 9272.1372
           sigsq = 0.195
         else if ( nmag .eq. 6000 ) then
           const = 10703.391
           sigsq = 0.767
         else
           nerr=2112
           call setmsg('ERROR',nerr)
           go to 8888
         end if
c
      else
c
        nerr=2112
        call setmsg('ERROR',nerr)
        go to 8888
c
      end if
c
      p = t1 / t2
      om1 = twopi / t1
         do 40 i = 1, nfreq
         omega = float (i-1) * delomg
         cs = cmplx ( 0.0, omega )
         ctd = cs**4 + 2.0 * om1 * ( h1 + p * h2 ) * cs**3
         ctd = ctd + ( ( 1.0 + p**2 ) + 4.0 * h1 * h2 * p *
     .                 ( 1.0 - sigsq ) ) * om1**2 * cs**2
         ctd = ctd + 2.0 * ( p * h1 + h2 ) * p * om1**3 * cs
         ctd = ctd + p**2 * om1**4
         ct = om1 * cs**3 / ctd
c
c THE PHASE RESPONSE IS DEFINED AS NEGATIVE THE PHASE RESPONSE OF
c HAGIWARA IN ORDER TO YIELD UPWARD FIRST MOTION ON THE SEISMOGRAM TRACE
c FOR AN IMPULSE OF GROUND DISPLACEMENT IN THE POSITIVE PHI DIRECTION.
c
	 ct = const * ct
         xre(i) = real ( ct )
         xim(i) = aimag ( ct )
   40    continue
c
 8888 return

      end
