c
c*******************************************************************************
c
c    Subroutine setfor
c
c    Author - Danny Harvey
c
c*******************************************************************************
c
      subroutine setfor (hue, light, sat)
c
      real*4 hue
      real*4 light
      real*4 sat
c
c    Subroutine setfor will set the foreground color for subsequent plotting.
c
c    Inputs  -	hue	= The hue of the foreground color. (0.0 to 360.0
c			  where 0.0 = red, 120.0 = green, 240.0 = blue)
c		light	= The lightness of the foreground color. (0.0 to 1.0
c			  where 0.0 is always black, 1.0 is always white and
c			  0.5 is a pure color when sat = 1.0. light is a gray
c			  scale when sat = 0.0)
c		sat	= The saturation of the foreground color. (0.0 to 1.0
c			  where 0.0 means no color (gray scale) and 1.0 means
c			  pure color when light = 0.5)
c
      common /npcolr/ fhue, flight, fsat, bhue, blight, bsat
c
      call hdfore (hue, light, sat)
      fhue = hue
      flight = light
      fsat = sat
      return
      end
c
      subroutine getfor (hue, light, sat)
c
      real*4 hue
      real*4 light
      real*4 sat
c
      common /npcolr/ fhue, flight, fsat, bhue, blight, bsat
c
      hue = fhue
      light = flight
      sat = fsat
c
      return
      end
