c
c*******************************************************************************
c
c    Subroutine setbac
c
c    Author - Danny Harvey
c
c*******************************************************************************
c
      subroutine setbac (hue, light, sat)
c
      real*4 hue
      real*4 light
      real*4 sat
c
c    Subroutine setbac will set the background color for subsequent plotting.
c
c    Inputs  -	hue	= The hue of the background color. (0.0 to 360.0
c			  where 0.0 = red, 120.0 = green, 240.0 = blue)
c		light	= The lightness of the background color. (0.0 to 1.0
c			  where 0.0 is always black, 1.0 is always white and
c			  0.5 is a pure color when sat = 1.0. light is a gray
c			  scale when sat = 0.0)
c		sat	= The saturation of the background color. (0.0 to 1.0
c			  where 0.0 means no color (gray scale) and 1.0 means
c			  pure color when light = 0.5)
c
      common /npcolr/ fhue, flight, fsat, bhue, blight, bsat
c
      call hdback (hue, light, sat)
      bhue = hue
      blight = light
      bsat = sat
      return
      end
      subroutine getbac (hue, light, sat)
      real*4 hue, light, sat
c
      common /npcolr/ fhue, flight, fsat, bhue, blight, bsat
c
      hue = bhue
      light = blight
      sat = bsat
      return
      end
