      subroutine rstn ( nfreq, delfrq, xre, xim, subtyp, nerr )
c
      include '../../inc/mach'
c
c   .....RSTN Seismometer.....
c
c   .....Poles and zeros due to D. Breding (1983), Data Users' Guide
c        for the Regional Seismic Test Network (RSTN), Sandia Report
c        SAND82-2935.....
c
      real * 8 delfrq, xre(1), xim(1)
      character * 8 subtyp
c
c
      if ( subtyp(1:2) .eq. 'CP' ) then
	 if ( subtyp(3:4) .eq. 'KL' ) then
	    if ( subtyp(5:6) .eq. '.Z' ) then
	       call clz ( nfreq, delfrq, xre, xim )
	    else if ( subtyp(5:6) .eq. '.N' .or. 
     .                subtyp(5:6) .eq. '.E' ) then
	       call clh ( nfreq, delfrq, xre, xim )
	    else
	       nerr = 2105
               call setmsg('ERROR',nerr)
               call apcmsg('RSTN:')
               call apcmsg(subtyp)
               go to 8888
	    end if
	 else if ( subtyp(3:4) .eq. 'KM' ) then
	    if ( subtyp(5:6) .eq. '.Z' ) then
	       call cmz ( nfreq, delfrq, xre, xim )
	    else if ( subtyp(5:6) .eq. '.N' .or. 
     .                subtyp(5:6) .eq. '.E' ) then
	       call cmh ( nfreq, delfrq, xre, xim )
	    else
	       nerr = 2105
               call setmsg('ERROR',nerr)
               call apcmsg('RSTN:')
               call apcmsg(subtyp)
               go to 8888
	    end if
	 else if ( subtyp(3:4) .eq. 'KS' ) then
	    if ( subtyp(5:6) .eq. '.Z' ) then
	       call csz ( nfreq, delfrq, xre, xim )
	    else if ( subtyp(5:6) .eq. '.N' .or. 
     .                subtyp(5:6) .eq. '.E' ) then
	       call csh ( nfreq, delfrq, xre, xim )
	    else
	       nerr = 2105
               call setmsg('ERROR',nerr)
               call apcmsg('RSTN:')
               call apcmsg(subtyp)
               go to 8888
	    end if
	 else
            nerr = 2105
            call setmsg('ERROR',nerr)
            call apcmsg('RSTN:')
            call apcmsg(subtyp)
            go to 8888
	 end if
      else if ( subtyp(1:2) .eq. 'NY' .or. subtyp(1:2) .eq. 'NT'
     .   .or. subtyp(1:2) .eq. 'ON' .or. subtyp(1:2) .eq. 'SD' ) then
	 if ( subtyp(3:4) .eq. 'KL' ) then
	    call rsl ( nfreq, delfrq, xre, xim )
	 else if ( subtyp(3:4) .eq. 'KM' ) then
	    call rsm ( nfreq, delfrq, xre, xim )
	 else if ( subtyp(3:4) .eq. '7S' ) then
	    call rs7 ( nfreq, delfrq, xre, xim, subtyp )
	 else if ( subtyp(3:4) .eq. 'KS' ) then
	    call rsk ( nfreq, delfrq, xre, xim, subtyp )
	 else
            nerr = 2105
            call setmsg('ERROR',nerr)
            call apcmsg('RSTN:')
            call apcmsg(subtyp)
            go to 8888
	 end if
      else
         nerr = 2105
         call setmsg('ERROR',nerr)
         call apcmsg('RSTN:')
         call apcmsg(subtyp)
         go to 8888
      end if
c
 8888 return
      end
