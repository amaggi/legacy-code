      subroutine sro ( nfreq, delfrq, xre, xim, subtyp )
c
      include '../../inc/mach'
c
c   .....SRO - for an SRO seismometer.....
c
      real * 8 delfrq, xre(1), xim(1)
      complex zero(8), pole(17)
      real * 4 const
      character * 8 subtyp
c
c
c   .....Set poles and zeros.....
c
c
c SRO SEISMOMETER.  THERE ARE THREE SYSTEMS THAT CAN BE COMPUTED
c HERE:  BROADBAND (BB), SHORT PERIOD (SP), OR LONG PERIOD (LP).
c REFERENCES ARE:
c MCCOWAN, D. W. AND R. T. LACOSS (1978).  TRANSFER FUNCTION
c FOR THE SEISMIC RESEARCH OBSERVATORY SEISMOGRAPH SYSTEM, BULL. SEIS.
c SOC. AM., 68, 501-512.
c BERGER, J., D. W. MCCOWAN, W. E. FARRELL, AND R. T. LACOSS (1978).
c COMMENTS ON 'TRANSFER FUNCTIONS FOR THE SEISMIC RESEARCH OBSERVATORY
c SEISMOGRAPH SYSTEM' BY DOUGLAS W. MCCOWAN AND RICHARD T. LACOSS, BULL.
c SEIS. SOC. AM., 68, 1537-1538.
c
c
c  .....SRO BROADBAND     'BB'.....
c
      if ( subtyp .eq. 'BB      ' ) then
         const = -394.0
c
         nzero = 4
         zero(1) = cmplx ( -0.125, 0.0 )
         zero(2) = cmplx ( -50.0, 0.0 )
         zero(3) = cmplx ( 0.0, 0.0 )
         zero(4) = cmplx ( 0.0, 0.0 )
c
	 npole = 4
         pole(1) = cmplx ( -0.13, 0.0 )
         pole(2) = cmplx ( -6.02, 0.0 )
         pole(3) = cmplx ( -8.66, 0.0 )
         pole(4) = cmplx ( -35.2, 0.0 )
c
c   .....SRO SHORT PERIOD.      'SP'.....
c INCLUDES SEISMOMETER AND SHAPING FILTER TRANSFER FUNCTIONS.
c
      else if ( subtyp .eq. 'SP    ' ) then
         const = 5.08233208 e13
c
	 nzero = 5
         zero(1) = cmplx ( -50.0, 0.0 )
         zero(2) = cmplx ( 0.0, 0.0 )
         zero(3) = cmplx ( 0.0, 0.0 )
         zero(4) = cmplx ( 0.0, 0.0 )
         zero(5) = cmplx ( 0.0, 0.0 )
c
         npole = 9
         pole(1) = cmplx ( -0.13, 0.0 )
         pole(2) = cmplx ( -6.02, 0.0 )
         pole(3) = cmplx ( -8.66, 0.0 )
         pole(4) = cmplx ( -35.2, 0.0 )
         pole(5) = cmplx ( -100.0, 0.0 )
         pole(6) = cmplx ( -17.97, 0.0 )
         pole(7) = cmplx ( -17.97, 0.0 )
         pole(8) = cmplx ( -63.29, 0.0 )
         pole(9) = cmplx ( -63.29, 0.0 )
c
c   .....SRO LONG PERIOD         'LP'.....
c INCLUDES SEISMOMETER, SHAPING FILTERS, AND ANTI-ALIASING FILTERS.
c
      else if ( subtyp .eq. 'LPDE    ' ) then
	 const = 16892.61226
c
	 nzero = 8
         zero(1) = cmplx ( -50.0, 0.0 )
         zero(2) = cmplx ( 0.0, 0.0 )
         zero(3) = cmplx ( 0.0, 0.0 )
         zero(4) = cmplx ( 0.0, 0.0 )
         zero(5) = cmplx ( 0.0, 1.05 )
         zero(6) = cmplx ( 0.0, -1.05 )
         zero(7) = cmplx ( 0.0, 0.0 )
         zero(8) = cmplx ( 0.0, 0.0 )
c
         npole = 17
         pole( 1) = cmplx ( -0.13, 0.0 )
         pole( 2) = cmplx ( -6.02, 0.0 )
         pole( 3) = cmplx ( -8.66, 0.0 )
         pole( 4) = cmplx ( -35.2, 0.0 )
         pole( 5) = cmplx ( -100.0, 0.0 )
         pole( 6) = cmplx ( -3.93, 0.0 )
         pole( 7) = cmplx ( -0.282, 0.0 )
         pole( 8) = cmplx ( -0.201, 0.241 )
         pole( 9) = cmplx ( -0.201, -0.241 )
         pole(10) = cmplx ( -0.134, 0.1 )
         pole(11) = cmplx ( -0.134, -0.1 )
         pole(12) = cmplx ( -0.0251, 0.0 )
         pole(13) = cmplx ( -0.00942, 0.0 )
         pole(14) = cmplx ( -0.24, 0.58 )
         pole(15) = cmplx ( -0.58, 0.24 )
         pole(16) = cmplx ( -0.58, -0.24 )
         pole(17) = cmplx ( -0.24, -0.58 )
c
      else
         nerr=2105
         call setmsg('ERROR',nerr)
         call apcmsg('SRO:')
         call apcmsg(subtyp)
         go to 8888
      end if
c
c   .....Compute transfer function.....
c
      call getran ( nfreq, delfrq, const, nzero, zero, npole,
     .              pole, xre, xim )
c
 8888 return
      end
