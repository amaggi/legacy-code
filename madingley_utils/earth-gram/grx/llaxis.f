c
c*******************************************************************************
c
c    Subroutine llaxis
c
c    Author - Danny Harvey
c
c*******************************************************************************
c
      subroutine llaxis(xdim,ydim,xmarg,ymarg,xlow,ylow,xmax,xmin,
     1                  ymax,ymin,labelx,labely,title,iclear)
c
c    routine llaxis will draw a box with tic marks, label the x and y axes
c    and label the plot with a title. This routine is similar to axis
c    except that llaxis will make a axes for a log-log plot.
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
c              labelx = character string label for the x-axis
c              labely = etc. for y-axis
c              title  = character string title which is placed on top
c			of the plot
c              iclear = clear flag for tektronix (if .ne. 0 then clear)
c
      character*(*) labelx,labely,title
c
      common  /ocflag/  imflag,iocflg,iltp
      common /axstuf/ ifont
      integer*4  ltp
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
      call chrsiz(.2,1.,0.)
      y2 = y + .2
      call text(.5,y2,0.,3,title,1)
      call chrsiz(.14,1.,0.)
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
      call tics(xmn,0.,xmx,0.,stic,ntics,1.,.17,0.02,iref)
      iref = -iref
      call tics(xmn,ydim,xmx,ydim,stic,ntics,1.,
     1  .17,0.02,iref)
      iref = -iref
      if (abs(stic-1.) .lt. .00001) nmin = nmin - 1
      if (abs(xmx-nmax-1.) .lt. .00001)
     1  nmax = nmax + 1
      do 50  i = nmin,nmax
      x = i
      call text(x,-.15,0.,8,'10',1)
      call chrsiz(.09,1.,0.)
      call number(x,-.16,0.,0,x,'(i4)')
      call chrsiz(.14,1.,0.)
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
      call setscl(0.,xdim,alog10(ymin),alog10(ymax))
      if (ymax .gt. ymin) then
        xmxx = ymax
        xmnn = ymin
        xmx = alog10(ymax)
        xmn = alog10(ymin)
        iref = -1
      else
        xmxx = ymin
        xmnn = ymax
        xmx = alog10(ymin)
        xmn = alog10(ymax)
        iref = 1
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
      if (ntics .le. 0)  go to 92
      stic = nmin - xmn
      call tics(0.,xmn,0.,xmx,stic,ntics,1.,.17,0.02,iref)
      iref = -iref
      call tics(xdim,xmn,xdim,xmx,stic,ntics,1.,
     1  .17,0.02,iref)
      iref = -iref
      if (abs(stic-1.) .lt. .00001) nmin = nmin - 1
      if (abs(xmx-nmax-1.) .lt. .00001)
     1  nmax = nmax + 1
      do 54  i = nmin,nmax
      x = i
      call text(-.33,x,0.,8,'10',1)
      call chrsiz(.09,1.,0.)
      call number(-.33,x,0.,0,x,'(i4)')
      call chrsiz(.14,1.,0.)
   54 continue
      call ntype('LIN','LOG')
      call setscl(0.,xdim,ymin,ymax)
      dx10 = 10.**(nmin-1)
      do 56  i = nmin,nmax+1
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
      dxe = 10.*dxs
      call tics(0.,dxs,0.,dxe,stic,ntics,dx,.10,0.,iref)
      iref = -iref
      call tics(xdim,dxs,xdim,dxe,stic,ntics,dx,
     1  .10,0.,iref)
      iref = -iref
   56 continue
      go to 200
   92 call number(-.1,ymin,0.,7,ymin,'(e12.5)')
      call number(-.1,ymax,0.,7,ymax,'(e12.5)')
  200 xx = -ymarg
      call ntype('LIN','LIN')
      call setscl(0.,xdim,ymin,ymax)
      yy = (ymin+ymax)/2.
      call text(xx,yy,90.,5,labely,1)
c
      if (xmin .ne. xmax .and. ymin .ne. ymax) then
        call ntype('LOG','LOG')
        call setscl(xmin,xmax,ymin,ymax)
      end if
c
  556 if (imfl .eq. 1) call nclose
c
      return
      end
