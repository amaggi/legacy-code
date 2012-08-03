c
c*******************************************************************************
c
c    Subroutine setbg
c
c    Author - Ben Kohl
c
c*******************************************************************************
c
      subroutine setbg (color)
c
      character*(*) color
c
c    Subroutine setbg will set the background color for subsequent plotting.
c
c    Inputs  -	color	= The name of the color of the background.  This
c			  name follows the naming convention of X11, e.g.
c			  "black", "white", "lightBlue", "#2fa".
c
      common /npcolr/ fhue, flight, fsat, bhue, blight, bsat
      real*4 hue
      real*4 light
      real*4 sat
      integer err
c
      err = hdbg (color, hue, light, sat)
      if ( err .eq. 1 ) return
      call hdfore( hue, light, sat )
      bhue = hue
      blight = light
      bsat = sat
      return
      end
