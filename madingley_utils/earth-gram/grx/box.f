c
c*******************************************************************************
c
c    Subroutine box
c
c    Author - Danny Harvey
c
c*******************************************************************************
c
      subroutine box(xleft,xright,ybot,ytop,thick,ithick,iclip)
c
c    routine box draws a box
c
c    inputs  - xleft  = horizontal location of left edge of box in user units
c	       xright = horizontal location of right edge of box in user units
c	       ybot   = vertical location of bottom edge of box in user units
c	       ytop   = vertical location of top edge of box in user units
c	       thick  = thickness of lines making up box in inches
c              ithick = thickness flag
c	       iclip  = clip flag see routine plot1
c
      dimension x(9),y(9)
c
      x(1) = (xright+xleft)/2.
      y(1) = ytop
      x(2) = xright
      y(2) = ytop
      x(3) = xright
      y(3) = (ytop+ybot)/2.
      x(4) = xright
      y(4) = ybot
      x(5) = x(1)
      y(5) = ybot
      x(6) = xleft
      y(6) = ybot
      x(7) = xleft
      y(7) = y(3)
      x(8) = xleft
      y(8) = ytop
      x(9) = x(1)
      y(9) = ytop
      call nplot(9,x,y,0,iclip,thick,ITHICK,' ')
c
      return
      end
