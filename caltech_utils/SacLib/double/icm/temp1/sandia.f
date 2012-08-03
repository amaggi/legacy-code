      subroutine sandia ( nfreq, delfrq, xre, xim, subtyp , nerr)
c
      include '../../inc/mach'
c
c   .....Sandia system 23 instrumental response.....
c
      real * 8 delfrq, xre(1), xim(1)
      real * 4 const
      complex zero(5), pole(7)
      character * 8 subtyp
c
c
c                           /* normalizes max (velocity) of Amplitude
      facnew = 10.37643443 e1
      const = 2.0 e4 * facnew
      nzero = 5
c                           /* Amplifier Transfer Function
      zero(1) = cmplx ( -189.50087, 0.0 )
	 do 20 i = 2, nzero
	 zero(i) = cmplx ( 0.0, 0.0 )
   20    continue
c
      npole = 7
c  poles 1 & 2 are amplifier transfer function
      pole(1) = cmplx ( -0.326726, 0.0 )
      pole(2) = cmplx ( -0.326726, 0.0 )
c  poles 3 to 6 are filter trandfer function
      pole(3) = cmplx ( -125.66371, 0.0 )
      pole(4) = cmplx ( -62.83185, -108.82477 )
      pole(5) = cmplx ( -62.83185, 108.82477 )
c
      if ( subtyp(1:1) .eq. 'O' ) then 
        if ( subtyp(2:2) .eq. 'L' ) then
	    srf = 1.22
	else if ( subtyp(2:2) .eq. 'B' ) then
	    srf = 1.26
	else if ( subtyp(2:2) .eq. 'D' .or. subtyp(2:2) .eq. 'T' ) then
	    srf = 1.27
	else if ( subtyp(2:2) .eq. 'N' .or. subtyp(2:2) .eq. 'E' ) then
	    srf = 1.23
        else
            nerr=2105
            call setmsg('ERROR',nerr)
            call apcmsg('SANDIA:')
            call apcmsg(subtyp)
            go to 8888
	end if
c
      else if ( subtyp(1:1) .eq. 'N' ) then
        if ( subtyp(2:3) .eq. 'NV' ) then
	    srf = 1.29
	else if ( subtyp(2:3) .eq. 'TV' .or. subtyp(2:3) .eq. 'BT'
     .       .or. subtyp(2:3) .eq. 'DV' ) then
	    srf = 1.27
	else if ( subtyp(2:3) .eq. 'TR' .or. subtyp(2:3) .eq. 'LT'
     .       .or. subtyp(2:3) .eq. 'DR' ) then
	    srf = 1.26
	else if ( subtyp(2:3) .eq. 'LR' .or. subtyp(2:3) .eq. 'DT'
     .       .or. subtyp(2:3) .eq. 'NR' ) then
	    srf = 1.25
	else if ( subtyp(2:3) .eq. 'BV' .or. subtyp(2:3) .eq. 'BR'
     .       .or. subtyp(2:3) .eq. 'NT' ) then
	    srf = 1.23
	else if ( subtyp(2:3) .eq. 'TT' .or. subtyp(2:3) .eq. 'LV' )
     .       then
		srf = 1.22
	else
            nerr=2105
            call setmsg('ERROR',nerr)
            call apcmsg('SANDIA:')
            call apcmsg(subtyp)
            go to 8888
	end if
      else
        nerr=2105
        call setmsg('ERROR',nerr)
        call apcmsg('SANDIA:')
        call apcmsg(subtyp)
        go to 8888
      end if
c
c   .....Allow for variable Seismometer Resonant Frequency ( srf ).....
c
c              2                            2
c           (s)  + 2h(2Pi*SRF)*s + (2pi*SRF)  = 0.0
c
c        where:
c              s = i * omega    and  h = seismometer damping factor
c
      s1 = -150.53419 * srf
      s2 = -0.2622555 * srf
c
      pole(6) = cmplx ( s1, 0.0 )
      pole(7) = cmplx ( s2, 0.0 )
c
c   .....Compute transfer function.....
c
      call getran ( nfreq, delfrq, const, nzero, zero, npole,
     .              pole, xre, xim )
c
 8888 return
      end
