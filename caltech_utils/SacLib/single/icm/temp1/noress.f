      subroutine noress ( nfreq, delfrq, xre, xim, subtyp )
c
      implicit complex ( c )
      include '../../inc/mach'
c
c   .....NORESS - for NORESS instrument....
c          (poles and zeros due to H. Durham)
c
      real * 8 delfrq, xre(1), xim(1), twopi
      character * 8 subtyp
c
      data twopi / 6.283185307179586 /
c
      nerr=0

      delomg = twopi * delfrq
c
c   .....computation based on LP, IP, or SP sub-type.....
c
      if ( subtyp .eq. 'LP      ' ) then
         do 20 i = 1, nfreq
         omega = float (i-1) * delomg
         cs = cmplx ( 0.0, omega )
	 cht = 12.66 / ( 1.0 + 0.02256 * cs + 0.6333 * cs**2 )
	 cha = -3360. / ( 1.0 + 0.003789 * cs )
         chb = 1.4215 e-4 * ( 1.0 + 0.01393 * cs ) /
     .                      ( 1.0 + 0.0009936 * cs )
         ckb = 1370.
         chka = cht * cha / ( 1. - cht * cha * chb * ckb )
         chkv = 100.0 * chka * cs
         chhp1 = -390.6 * cs / ( 1.0 + 7.813 * cs + 0.4505 * cs**2 )
         chbp1 = 324.0 * cs / ( 1.0 + 16.90 * cs + 113.3 * cs**2 )
         chdf = -1.0 / ( 1.0 + 0.001592 * cs )
         cha = 1.0 / ( 1.0 + 0.00159 * cs )
	 chb = 1.0 / ( 1.0 + 4.576 * cs + 11.66 * cs**2 )
	 chc = cs / ( ( 1.0 + 1504 * cs ) * ( 1.0 + 32.02 * cs ) )
	 che = 1.0 / ( 1.0 + 0.9477 * cs + 0.2428 * cs**2 )
	 chf = 1.0 / ( 1.0 + 0.6927 * cs + 0.2399 * cs**2 )
	 chg = 1.0 / ( 1.0 + 0.1913 * cs + 0.2403 * cs**2 )
	 chh = cs / ( ( 1.0 + 0.1913 * cs ) * ( 1.0 + 32.02 * cs ) )
	 gklp = 6589.6
	 chlpf = -gklp * cha * chb * chc * chb * che * chf * chg * chh
	 ch = chkv * chhp1 * chbp1 * chdf * chlpf * cs * 1.0 e-9
         xre(i) = real ( ch )
         xim(i) = aimag ( ch )
   20    continue
      return
c
      else if ( subtyp .eq. 'IP      ' ) then
         do 30 i = 1, nfreq
         omega = float (i-1) * delomg
         cs = cmplx ( 0.0, omega )
	 cht = 12.66 / ( 1.0 + 0.02256 * cs + 0.6333 * cs**2 )
	 cha = -3360. / ( 1.0 + 0.003789 * cs )
         chb = 1.4215 e-4 * ( 1.0 + 0.01393 * cs ) /
     .                      ( 1.0 + 0.0009936 * cs )
         ckb = 1370.
         chka = cht * cha / ( 1. - cht * cha * chb * ckb )
         chkv = 100.0 * chka * cs
         chbp2 = -10.17 * cs * ( 127.8 + cs ) /
     .                  ( 1.0 + 16.90 * cs + 113.1 * cs**2 )
         chdf = -1.0 / ( 1.0 + 0.001592 * cs )
         cha = 1.0 / ( 1.0 + 0.00159 * cs )
	 chb = 1.0 / ( 1.0 + 0.195 * cs + 0.00990 * cs**2 )
	 chc = cs / ( ( 1.0 + 0.00153 * cs ) * ( 1.0 + 5.11 * cs ) )
	 chd = 1.0 / ( 1.0 + 0.154 * cs + 0.00888 * cs**2 )
	 che = 1.0 / ( 1.0 + 0.0817 * cs + 0.0069 * cs**2 )
	 chf = cs / ( ( 1.0 + 0.00419 * cs ) * ( 1.0 + 5.11 * cs ) ) 
	 gkip = 129.5
	 chipf = -gkip * cha * chb * chc * chd * che * chf
	 ch = chkv * chbp2 * chdf * chipf * cs * 1.0 e-9
         xre(i) = real ( ch )
         xim(i) = aimag ( ch )
   30    continue
      return
c
      else if ( subtyp .eq. 'SP      ' ) then
         do 40 i = 1, nfreq
         omega = float (i-1) * delomg
         cs = cmplx ( 0.0, omega )
	 chs3 = 55.73 * cs**2 / ( 1.0 + 0.2387 * cs + 0.02533 * cs**2 )
	 chif = 0.9174 / ( 1.0 + 8.165 e-4 * cs + 7.339 e-7 * cs**2 )
         chpa = 5.544 * cs / ( ( 1.0 + 1.559 e-3 * cs ) *
     .                      ( 1.0 + 0.1598 * cs ) )
         go = 50.0
         chi = 1.0 / ( 1.0 + 7.950 e-4 * cs )
         ch1 = 1.0 / ( 1.0 + 0.0221 * cs + 0.000125 * cs**2 )
         ch2 = 1.0 / ( 1.0 + 0.0195 * cs + 0.000119 * cs**2 )
         ch3 = 1.0 / ( 1.0 + 0.0143 * cs + 0.000104 * cs**2 )
         ch4 = 1.0 / ( 1.0 + 0.00739 * cs + 0.0000820 * cs**2 )
         ch5 = 1.0 / ( 1.0 + 0.0150 * cs + 0.000113 * cs**2 )
         cho = 0.03186 * cs / ( ( 1.0 + 0.1593 * cs ) *
     .                          ( 1.0 + 1.416 e-3 * cs ) )
	 chsp = -go * chi * ch1 * ch2 * ch3 * ch4 * ch5 * cho
	 chdf = -1.0 / ( 1.0 + 0.001592 * cs )
	 ch = chs3 * chif * chpa * chsp * chdf * cs * 1.0 e-9
         xre(i) = real ( ch )
         xim(i) = aimag ( ch )
   40    continue
      return
      else
         nerr = 2105
         call setmsg('ERROR',nerr)
         call apcmsg('NORESS:')
         call apcmsg(subtyp)
	 return
      end if
c
      end
