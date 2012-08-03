c
c*******************************************************************************
c
c    Subroutine setctb
c
c    Author - Danny Harvey
c
c*******************************************************************************
c
      subroutine setctb (nc, h, l, s, ictg, ipixel, ierr)
c
      integer*4          nc
      real*4                 h(nc)
      real*4                    l(nc)
      real*4                       s(nc)
      integer*4                       ictg
      integer*4                             ipixel(nc)
      integer*4                                     ierr
c
c    setctb will attempt to define a color table. See setfor for definitions
c    of hue lightness and saturation.
c
c    Inputs  -	nc	= The number of entries in the color table.
c		h(nc)	= The hues of each color in the table.
c		l(nc)	= The lightnesses of each color in the table.
c		s(nc)	= The saturations of each color in the table.
c		ictg	= Contiguity flag.
c			  = 0 - The returned pixels need not be contiguous.
c			  = 1 - The returned pixels must be contiguous.
c
c    Outputs -	ipixel(nc)
c			= The pixel values corresponding to each requested
c			  color.
c		ierr	= An error flag.
c			  = 0 - No errors.
c			  = 1 - Too many colors requested.
c
      return
      end
