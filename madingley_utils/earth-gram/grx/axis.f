c
c*******************************************************************************
c
c    Subroutine axis
c
c    Author - Danny Harvey
c
c*******************************************************************************
c
      subroutine axis(xdim,ydim,xmarg,ymarg,xlow,ylow,xmax,xmin,
     1  ymax,ymin,dxsmal,dxnumb,dysmal,dynumb,
     2  fmtx,fmty,labelx,labely,title,iclear)
c
c    routine axis will draw a box with tic marks, label the x and y axes
c    and label the plot with a title.
c
c    inputs  - xdim   = x dimension of the box in inches
c              ydim   = y dimension of the box in inches
c              xmarg  = margin below the x-axis in inches - determines
c			the vertical spacing of the x-axis label
c              ymarg  = margin to the left of the y-axis in inches - determines
c			the horizontal spacing of the y-axis label
c              xlow   = x-location of the lower left hand corner of the box
c			in inches from the lower left hand corner of the
c			plot
c              ylow   = y-location etc.
c              xmax   = x value of the right edge of the box in user units
c              xmin   = x value of the left edge of the box in user units
c              ymax   = y value of the top of the box in user units
c              ymin   = y value of the bottom of the box in user units
c              dxsmal = increment between small tic marks along the
c			x-axis without number labelling
c              dxnumb = increment between large tic marks along the
c                       x-axis with number labelling
c              dysmal,dynumb = etc. for y-axis
c              fmtx   = character string format specification including
c			paranthesis which determines the x-axis numerical
c			labeling format
c              fmty   = etc. for the y-axis
c              labelx = character string label for the x-axis
c              labely = etc. for y-axis
c              title  = character string title which is placed on top
c			of the plot
c              iclear = clear flag.
c			= 0 - Dont clear plot area before drawing axes.
c			= 1 - Clear plot area before drawing axes.
c
      character*(*) labelx,labely,title,fmtx,fmty
c
      common  /ocflag/  imflag,iocflg,iltp
      integer*4  ltp
c
      common /npcolr/ fhue, flight, fsat, bhue, blight, bsat
c
      common /axstuf/ ifont
c
      common /ngrcom/ xsth, xsh, xsl, xss, ixss,
     +                xnth, xnh, xnl, xns, ixns,
     +                ysth, ysh, ysl, yss, iyss,
     +                ynth, ynh, ynl, yns, iyns
c
      data  ifont / 112 /
c
      data  xsth / -1.0 /
      data  xnth / -1.0 /
      data  ysth / -1.0 /
      data  ynth / -1.0 /
c
      imfl = imflag
      if (imflag .eq. 1) call nopen
      ltp = iltp
      CALL NLSTYL(0)
      call ntype('LIN','LIN')
c
c    clear region
c
      x = xdim   + .03
      y = ydim  + .03
      xl = xlow  - .015
      yl = ylow  - .015
      call setdim(x,y,xl,yl)
      call setscl(0.,1.,0.,1.)
      if (iclear .eq. 1) call clrrgn(0.,1.,0.,1.)
c
      call cfont(ifont)
      call chrsiz(.10,1.,0.)
c
c    x-axis
c
      call setdim(xdim,ydim,xlow,ylow)
      call setscl(xmin,xmax,0.,ydim)
C
C    SMALL TICS
C
      DDXBIG = DXNUMB
      if (dxsmal .le. 0.)  go to 101
C
C    FIRST DETERMINE "GOOD" INCREMENT
C
      DDXSML = DXSMAL
  109 NTICS = ABS(XMAX-XMIN)/DDXSML
      IF (NTICS .LT. 60)  GO TO 110
      DDXSML = DDXSML*2.
      CALL DNICE(DDXSML,DDXSML,DDXBIG)
      GO TO 109
  110 if (xmax .gt. xmin) then
        xmx = xmax
        xmn = xmin
        iref = 1
      else
        xmx = xmin
        xmn = xmax
        iref = -1
      end if
      if (xmn .ge. 0.) then
        xx = (xmn/DDXSML)*1.000001 + 1.
        nmin = xx
      else
        xx = (-xmn/DDXSML)*0.999999
        nmin = xx
        nmin = -nmin
      end if
      if (xmx .ge. 0.) then
        xx = (xmx/DDXSML)*0.999999
        nmax = xx
      else
        xx = (-xmx/DDXSML)*1.000001 + 1.
        nmax = xx
        nmax = -nmax
      end if
      ntics = nmax - nmin + 1
      if (ntics .le. 0)  go to 101
      stic = nmin*DDXSML - xmn
      if (xsth .ge. 0.0) then
        call getfor (xh, xl, xs)
        call setfor (xsh, xsl, xss)
        do 105  i = nmin, nmax
          x1 = i*DDXSML
          call line (x1, 0.0, x1, ydim, xsth, 0, 0)
  105   continue
        call setfor (xh, xl, xs)
      end if
      call tics(xmn,0.,xmx,0.,stic,ntics,DDXSML,.10,0.,iref)
      iref = -iref
      call tics(xmn,ydim,xmx,ydim,stic,ntics,DDXSML,
     1  .10,0.,iref)
  101 if (dxnumb .le. 0.)  go to 100
      if (xmax .gt. xmin) then
        xmx = xmax
        xmn = xmin
        iref = 1
      else
        xmx = xmin
        xmn = xmax
        iref = -1
      end if
      if (xmn .ge. 0.) then
        xx = (xmn/DDXBIG)*1.000001 + 1.
        nmin = xx
      else
        xx = (-xmn/DDXBIG)*0.999999
        nmin = xx
        nmin = -nmin
      end if
      if (xmx .ge. 0.) then
        xx = (xmx/DDXBIG)*0.999999
        nmax = xx
      else
        xx = (-xmx/DDXBIG)*1.000001 + 1.
        nmax = xx
        nmax = -nmax
      end if
      ntics = nmax - nmin + 1
      if (ntics .le. 0)  go to 90
      stic = nmin*DDXBIG - xmn
      if (xnth .ge. 0.0) then
        call getfor (xh, xl, xs)
        call setfor (xnh, xnl, xns)
        do 106  i = nmin, nmax
          x1 = i*DDXBIG
          call line (x1, 0.0, x1, ydim, xnth, 0, 0)
  106   continue
        call setfor (xh, xl, xs)
      end if
      call tics(xmn,0.,xmx,0.,stic,ntics,DDXBIG,.17,0.,iref)
      iref = -iref
      call tics(xmn,ydim,xmx,ydim,stic,ntics,DDXBIG,
     1  .17,0.,iref)
      if (abs(stic-DDXBIG)/DDXBIG .lt. .00001) nmin = nmin - 1
      if (abs(xmx-nmax*DDXBIG-DDXBIG)/DDXBIG .lt. .00001)
     1  nmax = nmax + 1
      do 50  i = nmin,nmax
      x = DDXBIG*i
      call number(x,-.1,0.,5,x,fmtx)
   50 continue
      go to 100
   90 call number(xmin,-.1,0.,5,xmin,fmtx)
      call number(xmax,-.1,0.,5,xmax,fmtx)
      ix = xmax
  100 xx = -xmarg
      yy = (xmin+xmax)/2.
      call text(yy,xx,0.,3,labelx,1)
c
c    y-axis
c
      call setscl(0.,xdim,ymin,ymax)
C
C    SMALL TICS
C
      DDYBIG = DYNUMB
      if (dysmal .le. 0.)  go to 201
C
C    FIRST DETERMINE "GOOD" INCREMENT
C
      DDYSML = DYSMAL
  209 NTICS = ABS(YMAX-YMIN)/DDYSML
      IF (NTICS .LT. 60)  GO TO 210
      DDYSML = DDYSML*2.
      CALL DNICE(DDYSML,DDYSML,DDYBIG)
      GO TO 209
  210 if (ymax .gt. ymin) then
        xmx = ymax
        xmn = ymin
        iref = -1
      else
        xmx = ymin
        xmn = ymax
        iref = 1
      end if
      if (xmn .ge. 0.) then
        xx = (xmn/DDYSML)*1.000001 + 1.
        nmin = xx
      else
        xx = (-xmn/DDYSML)*0.999999
        nmin = xx
        nmin = -nmin
      end if
      if (xmx .ge. 0.) then
        xx = (xmx/DDYSML)*0.999999
        nmax = xx
      else
        xx = (-xmx/DDYSML)*1.000001 + 1.
        nmax = xx
        nmax = -nmax
      end if
      ntics = nmax - nmin + 1
      if (ntics .le. 0)  go to 201
      stic = nmin*DDYSML - xmn
      if (ysth .ge. 0.0) then
        call getfor (xh, xl, xs)
        call setfor (ysh, ysl, yss)
        do 205  i = nmin, nmax
          x1 = i*DDYSML
          call line (0.0, x1, xdim, x1, ysth, 0, 0)
  205   continue
        call setfor (xh, xl, xs)
      end if
      call tics(0.,xmn,0.,xmx,stic,ntics,DDYSML,.10,0.,iref)
      iref = -iref
      call tics(xdim,xmn,xdim,xmx,stic,ntics,DDYSML,
     1  .10,0.,iref)
  201 if (dynumb .le. 0.)  go to 200
      if (ymax .gt. ymin) then
        xmx = ymax
        xmn = ymin
        iref = -1
      else
        xmx = ymin
        xmn = ymax
        iref = 1
      end if
      if (xmn .ge. 0.) then
        xx = (xmn/DDYBIG)*1.000001 + 1.
        nmin = xx
      else
        xx = (-xmn/DDYBIG)*0.999999
        nmin = xx
        nmin = -nmin
      end if
      if (xmx .ge. 0.) then
        xx = (xmx/DDYBIG)*0.999999
        nmax = xx
      else
        xx = (-xmx/DDYBIG)*1.000001 + 1.
        nmax = xx
        nmax = -nmax
      end if
      ntics = nmax - nmin + 1
      if (ntics .le. 0)  go to 290
      stic = nmin*DDYBIG - xmn
      if (ynth .ge. 0.0) then
        call getfor (xh, xl, xs)
        call setfor (ynh, ynl, yns)
        do 206  i = nmin, nmax
          x1 = i*DDYBIG
          call line (0.0, x1, xdim, x1, ynth, 0, 0)
  206   continue
        call setfor (xh, xl, xs)
      end if
      call tics(0.,xmn,0.,xmx,stic,ntics,DDYBIG,.17,0.,iref)
      iref = -iref
      call tics(xdim,xmn,xdim,xmx,stic,ntics,DDYBIG,
     1  .17,0.,iref)
      if (abs(stic-DDYBIG)/DDYBIG .lt. .00001) nmin = nmin - 1
      if (abs(xmx-nmax*DDYBIG-DDYBIG)/DDYBIG .lt. .00001)
     1  nmax = nmax + 1
      do 250  i = nmin,nmax
      x = DDYBIG*i
      call number(-.1,x,0.,7,x,fmty)
  250 continue
      go to 200
  290 call number(ymin,-.1,0.,5,ymin,fmty)
      call number(ymax,-.1,0.,5,ymax,fmty)
  200 xx = -ymarg
      yy = (ymin+ymax)/2.
      call text(xx,yy,90.,5,labely,1)
c
c    box and title
c
      x = xdim   + .03
      y = ydim  + .03
      xl = xlow  - .015
      yl = ylow  - .015
      call setdim(x,y,xl,yl)
      call setscl(0.,1.,0.,1.)
      call box(0.,1.,0.,1.,.03,0,1)
c
      call cfont(ifont)
      call setdim(xdim,ydim,xlow,ylow)
      call setscl(0.,1.,0.,y)
      call chrsiz(.15,1.,0.)
      y2 = y + .1
      call text(.5,y2,0.,3,title,1)
      call chrsiz(.10,1.,0.)
c
      call setscl(xmin,xmax,ymin,ymax)
c
      if (imfl .eq. 1) call nclose
      call hdstrk
      call hdstrkl
c
      return
      end
c
c*******************************************************************************
c
c    Subroutine number
c
c    Author - Danny Harvey
c
c*******************************************************************************
c
      subroutine number(x,y,ang,iref,xnum,fmt)
c
c    routine number causes a number to be printed on a plot
c
      character*(*) fmt
      character*40 numb, numb2
c
      if (fmt .eq. '(none)') return
      if (xnum .lt. 0.0) then
        xn = xnum*1.000001
        ix = xn - 0.5
      else
        xn = xnum*1.000001
        ix = xn + 0.5
      end if
      if (fmt .eq. '(*)') then
	write (numb,'(f20.10)') xn
	j = 0
	numb2 = ' '
	jcount = 6
	jdot = 0
	do 10  i = 1, 40
	  if (numb(i:i) .ne. ' ') then
	    j = j + 1
	    numb2(j:j) = numb(i:i)
	    if (numb2(j:j) .eq. '-') jcount = jcount + 1
	    if (numb2(j:j) .eq. '.') then
	      jdot = j
	      jcount = jcount + 1
	    end if
	    if (j .gt. jcount .and. j .ne. jdot) numb2(j:j) = '0'
	  end if
   10   continue
	jlast = j
	do 20  i = jlast, jdot, -1
	  if (numb2(i:i) .ne. '0') go to 21
	  numb2(i:i) = ' '
   20   continue
   21   continue
	if (numb2(i:i) .eq. '.') numb2(i:i) = ' '
	if (numb2 .eq. ' ') numb2 = '0'
        call text(x,y,ang,iref,numb2,1)
	return
      end if
      if (fmt(2:2) .eq. 'i' .or. fmt(2:2) .eq. 'I')
     1   write (numb,fmt) ix
      if (fmt(2:2) .ne. 'i' .and. fmt(2:2) .ne. 'I')
     1   write (numb,fmt) xn
      j = 0
      do 100  k = 1,20
      if (numb(k:k) .ne. ' ')  go to 110
      j = j + 1
  100 continue
  110 if (j .eq. 0)  go to 120
      do 130  k = 1,j
      do 130  m = 2,20
  130 numb(m-1:m-1) = numb(m:m)
  120 continue
      call text(x,y,ang,iref,numb,1)
c
      return
      end
      SUBROUTINE DNICE(DX, DS, DB)
C
C    THIS SUBROUTINE WILL DETERMINE "NICE" LOOKING INCREMENTS
C
      DDX = ABS(DX)
      DMULT = 1.
      DDDX = DDX
   10 IF (DDDX .LT. 10.)  GO TO 20
      DMULT = DMULT*.1
      DDDX = DDX*DMULT
      GO TO 10
   20 IF (DDDX .GE. 1.)  GO TO 30
      DMULT = DMULT*10.
      DDDX = DDX*DMULT
      GO TO 20
   30 IF (DDDX .LT. 7.5)  GO TO 40
      DS = 10./DMULT
      DB = 10.*DS
      GO TO 100
   40 IF (DDDX .LT. 3.5)  GO TO 50
      DS = 5./DMULT
      DB = 5.*DS
      GO TO 100
   50 IF (DDDX .LT. 1.5)  GO TO 60
      DS = 2./DMULT
      DB = 5.*DS
      GO TO 100
   60 DS = 1./DMULT
      DB = 10.*DS
  100 RETURN
      END
      subroutine setaxf (jfont)
      common /axstuf/ ifont
      ifont = jfont
      return
      end
