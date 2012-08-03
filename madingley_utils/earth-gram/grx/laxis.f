c
c*******************************************************************************
c
c    Subroutine laxis
c
c    Author - Danny Harvey
c
c*******************************************************************************
c
      subroutine laxis(xdim,ydim,xmarg,ymarg,xlow,ylow,xmax,xmin,
     1                  ymax,ymin,dysmal,dynumb,fmty,
     2                  labelx,labely,title,iclear)
c
c    routine laxis will draw a box with tic marks, label the x and y axes
c    and label the plot with a title. This routine is similar to axis
c    except that llaxis will make a axes for a log-lin plot.
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
c              dysmal = increment between small tic marks along the
c			y-axis without number labelling. If <= 0, then
c			no small tics are drawn.
c              dynumb = increment between large tic marks along the
c                       y-axis with number labelling. If <= 0, then
c			no large tics are drawn.
c              fmty   = character string format specification including
c			paranthesis which determines the y-axis numerical
c			labeling format. If '(*)', then the format is
c			automatically determined.
c              labelx = character string label for the x-axis
c              labely = etc. for y-axis
c              title  = character string title which is placed on top
c			of the plot
c              iclear = clear flag for tektronix (if .ne. 0 then clear)
c
      character*(*) labelx,labely,title,fmty
c
      common  /ocflag/  imflag,iocflg,iltp
      integer*4  ltp
      common /axstuf/ ifont
c
      common /ngrcom/ xsth, xsh, xsl, xss, ixss,
     +                xnth, xnh, xnl, xns, ixns,
     +                ysth, ysh, ysl, yss, iyss,
     +                ynth, ynh, ynl, yns, iyns
c
c    clear if appropriate
c
      imfl = imflag
      if (imflag .eq. 1) call nopen
      ltp = iltp
      CALL NLSTYL(0)
      call ntype('LIN','LIN')
c
c    box first
c
      x = xdim   + .03
      y = ydim  + .03
      xl = xlow  - .015
      yl = ylow  - .015
      call setdim(x,y,xl,yl)
      call setscl(0.,1.,0.,1.)
      if (iclear .eq. 1) call clrrgn(0.,1.,0.,1.)
      call box(0.,1.,0.,1.,.03,0,1)
c
c    title next
c
      call cfont(ifont)
      call setdim(xdim,ydim,xlow,ylow)
      call setscl(0.,1.,0.,y)
      call chrsiz(.15,1.,0.)
      y2 = y + .2
      call text(.5,y2,0.,3,title,1)
      call chrsiz(.10,1.,0.)
c
c    x-axis
c
      if (xmax .eq. xmin) go to 555
      call ntype('LIN','LIN')
      call setdim(xdim,ydim,xlow,ylow)
      call setscl(alog10(xmin),alog10(xmax),0.,ydim)
      if (xmax .gt. xmin) then
        xmxx = xmax
        xmnn = xmin
        xmx = alog10(xmax)
        xmn = alog10(xmin)
        iref = 1
      else
        xmxx = xmin
        xmnn = xmax
        xmx = alog10(xmin)
        xmn = alog10(xmax)
        iref = -1
      end if
      if (xmn .ge. 0.) then
        xx = xmn*1.000001 + 1.
        nmin = xx
      else
        xx = (-xmn)*0.999999
        nmin = xx
        nmin = -nmin
      end if
      if (xmx .ge. 0.) then
        xx = xmx*0.999999
        nmax = xx
      else
        xx = (-xmx)*1.000001 + 1.
        nmax = xx
        nmax = -nmax
      end if
      ntics = nmax - nmin + 1
      if (ntics .le. 0)  go to 90
      stic = nmin - xmn
      if (xnth .ge. 0.0) then
        call getfor (xh, xl, xs)
        call setfor (xnh, xnl, xns)
        do 106  i = nmin, nmax
          x1 = i
          call line (x1, 0.0, x1, ydim, xnth, 0, 0)
  106   continue
        call setfor (xh, xl, xs)
      end if
      call tics(xmn,0.,xmx,0.,stic,ntics,1.,.17,0.0,iref)
      iref = -iref
      call tics(xmn,ydim,xmx,ydim,stic,ntics,1.,
     1  .17,0.0,iref)
      iref = -iref
      if (abs(stic-1.) .lt. .00001) nmin = nmin - 1
      if (abs(xmx-nmax-1.) .lt. .00001)
     1  nmax = nmax + 1
      do 50  i = nmin,nmax
      x = i
      call text(x,-.15,0.,8,'10',1)
      call chrsiz(.06,1.,0.)
      call number(x,-.16,0.,0,x,'(i4)')
      call chrsiz(.10,1.,0.)
   50 continue
      call ntype('LOG','LIN')
      call setscl(xmin,xmax,0.,ydim)
      dx10 = 10.**(nmin-1)
      do 52  i = nmin,nmax+1
      dx = dx10
      dx10 = 10.*dx
      ntics = xmxx/dx
      if (ntics .gt. 9) ntics = 9
      dxs = dx
      if (dx .ge. xmnn) then
        stic = 0.
      else
        mtics = xmnn/dx
        stic = mtics
        ntics = ntics - mtics
        stic = stic*dx
      end if
      dxe = dxs*10.
      if (xsth .ge. 0.0) then
        call getfor (xh, xl, xs)
        call setfor (xsh, xsl, xss)
        do 105  j = 2, 9
          x1 = j*dxs
          call line (x1, 0.0, x1, ydim, xsth, 0, 0)
  105   continue
        call setfor (xh, xl, xs)
      end if
      call tics(dxs,0.,dxe,0.,stic,ntics,dx,.10,0.,iref)
      iref = -iref
      call tics(dxs,ydim,dxe,ydim,stic,ntics,dx,
     1  .10,0.,iref)
      iref = -iref
   52 continue
      go to 100
   90 call number(xmin,-.1,0.,5,xmin,'(e12.5)')
      call number(xmax,-.1,0.,5,xmax,'(e12.5)')
      ix = xmax
  100 xx = -xmarg
      call ntype('LIN','LIN')
      call setscl(xmin,xmax,0.,ydim)
      yy = (xmin+xmax)/2.
      call text(yy,xx,0.,3,labelx,1)
c
c    y-axis
c
  555 if (ymax .eq. ymin) go to 556
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
      if (xmin .ne. xmax .and. ymin .ne. ymax) then
        call ntype('LOG','LIN')
        call setscl(xmin,xmax,ymin,ymax)
      end if
c
  556 if (imfl .eq. 1) call nclose
c
      return
      end
