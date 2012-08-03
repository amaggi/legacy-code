c
c*******************************************************************************
c
c    Subroutine cursor
c
c    Author - Danny Harvey
c
c*******************************************************************************
c
      SUBROUTINE CURSOR(X, Y, CHR)
C
C    DJH  3/29/83
C
C    Subroutine CURSOR will enable the cross-hair graphics input mode and
C    read in the position of the cursor. Note that the appearance of the
C    cursor is controlled by a prior call to setcursor() and cursor tracking
C    can be done with a prior call to trackcursor().
C
C    inputs -  X      = X-coordinate in user units for prepositioning cursor
C              Y      = Y-coordinate in user units for prepositioning cursor
C
C    outputs - X      = X-coordinate in user units of input cursor position
C              Y      = Y-coordinate in user units of input cursor position
C              CHR    = character which was typed to initiate read
C
      CHARACTER*(*) CHR
c
      common /pdim/ xdim,ydim,xlow,ylow,rxdim,rydim,rxlow,rylow,
     1              xbl,ybl,xbh,ybh,xbm,ybm,itran,tangle,
     2              ca,sa,cellht,cellwd,ixtype,iytype
c
      common /pscl/ xmin,xmax,ymin,ymax,xscale,yscale,xrange,yrange
c
      common /spc/ xll,yll,xur,yur,ASPECT,xc,xcl,fl
c
      common /ccom/ itype, width, height
c
c
c    Pre-position the pointer
c
      xpre = xmap(x) + xll
      ypre = ymap(y) + yll
      ixmin = rxlow + 0.5
      ixmax = rxlow + rxdim + 0.5
      iymax = yur - (rylow + 0.5) + 1
      iymin = yur - (rylow + rydim + 0.5) + 1
      ypre = yur - ypre
      if (itype .eq. 1) then
        iwidth = width*abs(xscale) + 0.5
        iheight = 0
      else if (itype .eq. 2) then
        iwidth = 0
        iheight = height*abs(yscale) + 0.5
      else if (itype .eq. 3) then
        iwidth = width*abs(xscale) + 0.5
        iheight = height*abs(yscale) + 0.5
      else
        iwidth = 0
        iheight = 0
      end if
      if (xpre-iwidth/2 .lt. xll) xpre = xll+iwidth/2
      if (ypre-iheight/2 .lt. yll) ypre = yll+iheight/2
      if (xpre+iwidth/2 .gt. xur) xpre = xur-iwidth/2
      if (ypre+iheight/2 .gt. yur) ypre = yur-iheight/2
      ixpre = xpre + 0.5
      iypre = ypre + 0.5
      call hdmptr (ixmin,ixmax,iymin,iymax,ixpre,iypre)
c
c    Get the new pointer position.
c
      xsc = 1.0/xscale
      ysc = -1.0/yscale
      xmn = xmin
      ymn = ymax
      call hdgptr (itype,iwidth,iheight,ixmin,ixmax,xsc,xmn,
     +             iymin,iymax,ysc,ymn,ix,iy,chr)
      xpos = ix
      ypos = iy
      ypos = yur - ypos
      xras = xpos - xll
      yras = ypos - yll
      if (ixtype .eq. 0) then
	x = (xras-rxlow)/xscale + xmin
      else
	x = 10.0**((xras-rxlow)/xscale + xmin)
      end if
      if (iytype .eq. 0) then
	y = (yras-rylow)/yscale + ymin
      else
	y = 10.0**((yras-rylow)/yscale + ymin)
      end if
C
C    Normal exit
C
      RETURN
C
C    End of subroutine CURSOR
C
      END
c
c*******************************************************************************
c
c    Subroutine setcursor
c
c    Author - Danny Harvey
c
c*******************************************************************************
c
      subroutine setcursor(type, param1, param2)
C
C    DJH  3/29/83
C
C    Subroutine setcursor will set the appearance of the graphics cursor.
C    This should be called prior to cursor().
C
C    inputs -  type   = Cursor type.
C			- crosshairs	= Normal crosshairs (default).
C			- x_window	= Two vertical lines centered about
C					  the cursor defining a
C					  horizontal window.
C			- y_window	= Two horizontal lines centered about
C					  the cursor defining a
C					  vertical window.
C			- box		= A rectangle centered on the cursor.
C	       param1, param2
C		      = Associated parameters according to type.
C			type		definitions
C			crosshairs	Not used.
C			x_window	param1 = horizontal window width.
C			y_window	param1 = vertical window width.
C			box		param1 = horizontal box width.
C					param2 = vertical box width.
C
      CHARACTER*(*) type
c
      common /ccom/ itype, width, height
c
      data  itype / 0 /
      data  width / 0.0 /
      data  height / 0.0 /
c
      if (type .eq. 'x_window') then
	itype = 1
	width = param1
      else if (type .eq. 'y_window') then
	itype = 2
	height = param1
      else if (type .eq. 'box') then
	itype = 3
	width = param1
	height = param2
      else
	itype = 0
      end if
C
C    Normal exit
C
      RETURN
C
C    End of subroutine setcursor
C
      END
c
c*******************************************************************************
c
c    Subroutine trackcursor
c
c    Author - Danny Harvey
c
c*******************************************************************************
c
      subroutine trackcursor(callback)
C
C    DJH  3/29/83
C
C    Subroutine trackcursor will cause a user supplied callback subroutine
C    to be called whenever the cursor is moved or drawn.  This should be 
C    called prior to cursor(). A call to trackcursoroff() should be made
C    to disable cursor tracking.
C
C    inputs -  callback   = Callback subroutine which must be declared
C			    as an EXTERNAL in the calling program.
C
      call hdtptr (callback)
C
C    Normal exit
C
      RETURN
C
C    End of subroutine trackcursor
C
      END
c
c*******************************************************************************
c
c    Subroutine trackcursoroff
c
c    Author - Danny Harvey
c
c*******************************************************************************
c
      subroutine trackcursoroff
C
C    DJH  3/29/83
C
C    Subroutine trackcursoroff will disable cursor tracking.
C
      call hdtptroff
C
C    Normal exit
C
      RETURN
C
C    End of subroutine trackcursoroff
C
      END
