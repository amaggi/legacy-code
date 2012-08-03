c
c*******************************************************************************
c
c    Subroutine chrdir
c
c    Author - Danny Harvey
c
c*******************************************************************************
c
      subroutine chrdir(angle)
c
c    Subroutine CHRDIR will set the direction or orientation of text
c    labelling.
c
c    input  - angle  = the angle in degrees measured positive
c                      counter-clockwise from a horizontal line which
c                      defines the orientation of a line of text
c                      labelling
c
      common /partxt/ sinang,cosang,tansln,hite,rat
c
      data  rad     /  0.0174532925  /
c
      a = angle*rad
      sinang = sin(a)
      cosang = cos(a)
c
      return
c
c    end of Subroutine CHRDIR
c
      end
