c
c*******************************************************************************
c
c    Subroutine setfg
c
c    Author - Ben Kohl
c
c*******************************************************************************
c
      subroutine setfg (color)
c
      character*(*) color
c
c    Subroutine setfg will set the foreground color for subsequent plotting.
c
c    Inputs  -	color	= The name of the color of the foreground.  This
c			  name follows the naming convention of X11, e.g.
c			  "black", "white", "lightBlue", "#2fa".
c
      common /npcolr/ fhue, flight, fsat, bhue, blight, bsat
      real*4 hue
      real*4 light
      real*4 sat
      integer err
c
      err = hdfg (color, hue, light, sat)
      if ( err .eq. 1 ) return
      call hdfore( hue, light, sat )
      fhue = hue
      flight = light
      fsat = sat
      return
      end
