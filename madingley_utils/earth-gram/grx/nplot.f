c
c*******************************************************************************
c
c    Subroutine nplot
c
c    Author - Danny Harvey
c
c*******************************************************************************
c
      subroutine nplot(n,x,y,igraf,iclip,thick,ithick,asymb)
c
c    routine nplot will plot a set of points specified by the arrays x
c    and y
c
c    inputs  - n      = number of points to plot
c              x      = plot data horizontal coordinate array
c              y      = plot data vertical coordinate array
c              igraf  = plot flag
c			< 0 - plot only symbols with no lines
c			= 0 - plot only lines with no symbols
c			> 0 - plot both lines and symbols
c	       iclip  = clipping flag see routine plot1
c              thick  = thickness of line in inches on the versatek
c              ithick = flag which controls whether line is thickened
c			symeterically about the center or asymetrically
c			= 0 - symetrically
c			.gt. 0 - thicken to the right in the direction of
c				 plotting
c			.lt. 0 - thicken to the left
c              asymb  = plotting symbol.  note:  if the number of points,
c			n, is negative then it as assumed that asymb
c			is an array of dimension n
c
      common /spc/ xll,yll,xur,yur,ASPECT,xc,xcl,fl
c
      common /pdim/ xdim,ydim,xlow,ylow,rxdim,rydim,rxlow,rylow,
     1              xbl,ybl,xbh,ybh,xbm,ybm,itran,tangle,
     2              ca,sa,cellht,cellwd,ixtype,iytype
c
      common /pscl/ xmin,xmax,ymin,ymax,xscale,yscale,xrange,yrange
c
      common /filcom/ ifill, jfill
c
      character*1 asymb(*)
      dimension x(*),y(*)
c
      common  /ocflag/  imflag,iocflg,iltp
      integer*4  ltP,ICHR
c
c    rmin = minimum spacing in standard rasters for line thickening
c
      data  rmin  / 1.0 /
c
      nn = iabs(n)
c
c    line plot
c
      imfl = imflag
      imflag = 0
      if (imfl .eq. 1)  CALL NOPEN
      ltp = iltp
      if (igraf .lt. 0)  go to 200
      rthick = thick*xc
      ilines = rthick/rmin
      if (ilines .eq. 0) ilines = 1
      rthick = (ilines-1)*rmin/2.
      if (ithick .gt. 0) rthick = 0.
      if (ithick .lt. 0) rthick = (ilines-1)*rmin
      ilabs = ilines
      if (ilines .lt. 0) ilabs = -thick
c
      call nlstyl(ltp)
      rlthick = thick*xcl/ilabs
      nth = rlthick + 0.5
      if (ifill .eq. 1) then
        call hdthik (0)
        roff = 0.0
        jfill = 1
        call plot1(nn, x, y, roff, iclip)
        jfill = 0
        call hdfill (1)
        call hdfilll (1)
      end if
      call hdthik (nth)
      do 100  i = 1,ilabs
      roff = rmin*(i-1) - rthick
      if (ilines .lt. 0) roff = 0.
c
c    loop for data gaps
c
      j1 = 1
      do 110  j = 1,nn
      if (y(j) .lt. 1.e30)  go to 110
      if (j1 .eq. j)  go to 115
      np = j - j1
      call plot1(np, x(j1), y(j1), roff, iclip)
      j1 = j
  115 j1 = j1 + 1
  110 continue
      np = nn + 1 - j1
  100 call plot1(np,x(j1),y(j1),roff,iclip)
c
c    symbol plot
c
  200 if (igraf .eq. 0)  go to 900
C      call jlstyl(0)
C      CALL JJUST(2,2)
C      CALL JSIZE(100.,70.)
      do 300  i = 1,nn
      j = 1
      if (n .lt. 0) j = i
      xx = xmap(x(I))
      yy = ymap(y(I))
c
c    clip if desired
c
      if (iclip .ne. 0)  go to 310
      if (xx .lt. xbl) GO TO 300
      if (xx .gt. xbh) GO TO 300
      if (yy .lt. ybl) GO TO 300
      if (yy .gt. ybh) GO TO 300
  310 if (xx .lt. 0.) GO TO 300
      if (xx .gt. xbm) GO TO 300
      if (yy .lt. 0.) GO TO 300
      if (yy .gt. ybm) GO TO 300
      if (asymb(j) .eq. ' ')  go to 300
      CALL MOVEV(XX, YY, ITRAN)
      ICHR = ICHAR(ASYMB(J))
C      CALL J2TEXT(1,ICHR)
      call text(X(I),Y(I),0.,4,asymb(j),iclip)
  300 continue
      call nlstyl(ltp)
c
  900 if (imfl .eq. 1)  call NCLOSE
c     IF (IMFL .EQ. 1)  CALL HOME
      imflag = imfl
      call hdstrk
      call hdstrkl
      call hdthik (0)
      return
      end
