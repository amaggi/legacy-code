c
c*******************************************************************************
c
c    Subroutine initt
c
c    Author - Danny Harvey
c
c*******************************************************************************
c
      subroutine initt(itran, plotfile, display, program, 
     1                        ssize, xwin, ywin)
c
c    routine initt is an initialization routine which MUST BE CALLED prior
c    to any of the other plot calls.  It need be called only once to
c    properly initialize the connection to the X server and the
c    postscript plot file.
c
c    inputs  - itran  	= flag which indicates which mode for plotting
c			  = 0 - Portrait mode
c			  .ne. 0 - Landscape mode
c              plotfile	= The name of a postscript plot file. If this string
c			  is empty, then a default file named "plotfile"
c			  will be created. If this string is "none", then
c			  no postcript output will be produced.
c              display	= The name of the X server to connect to. If this
c			  string is empty, then the the environment
c			  variable DISPLAY is used to define the X server.
c              program	= A program name that is displyed in the X window
c			  title bar and the icon.
c              ssize	= The initial size of the window in screen units
c			  (0.0 - 1.0). In portrait mode this is the height of
c			  the screen and in landscape mode this is the width
c			  of the window.
c              xwin,ywin
c			= The initial x-y location of the window in screen
c			  units. 
c
      character*(*) plotfile, display, program
c
      common /pdim/ xdim,ydim,xlow,ylow,rxdim,rydim,rxlow,rylow,
     1              xbl,ybl,xbh,ybh,xbm,ybm,iitran,tangle,
     2              ca,sa,cellht,cellwd,ixtype,iytype
c
      common /spc/ xll,yll,xur,yur,ASPECT,xc,xcl,fl
c
      common  /ocflag/  imflag,iocflg,iltp
c
      common /xplpid/ ipid
c
      common /filcom/ ifill, jfill
c
      character*80 pscolor
c
      data  imflag  /  1  /
      data  iocflg  /  0  /
      data  iltp  /  0  /
c     data  scrw, scrh / 1280.0, 1024.0 /
c     data  scrw, scrh / 1365.0, 1024.0 /
      data  scrw, scrh / 1363.0, 1024.0 /
      data  rlaser / 300.0 /
      data  ifirst / 1 /
      data  ifill  / 0 /
      data  jfill  / 0 /
c
      hsize = ssize
      if (itran .eq. 0) then
        wsize = ssize*(scrh/scrw)
      else
        wsize = ssize*(scrw/scrh)
      end if
      call hdopen(plotfile, display, program, xwin, ywin,
     2                        wsize, hsize, wd, ht, ipid)
      xll = 0.0
      yll = 0.0
      xur = wd - 1.0
      yur = ht - 1.0
      iitran = itran
      tangle = 0.
      if (itran .ne. 0) tangle = 90.
      if (itran .eq. 0) then
        xc = (yur - yll) / 10.0
      else
        xc = (xur - xll) / 10.0
      end if
      xcl = rlaser
      if (plotfile .eq. 'none') xcl = 0.0
      fl = xcl / xc
      CALL CFONT(1)
      call ntype('LIN','LIN')
      call chrdir(0.)
      call setdim(7.,7.,.5,.5)
      call setscl(0.,7.,0.,7.)
      call chrsiz(.14,1.,0.)
      call hdinit(itran, plotfile, program)
      call hdfore (0.0, 0.0, 0.0)
      call hdback (0.0, 1.0, 0.0)
      call getenv ('GRX_PSCOLOR', pscolor)
      if (pscolor .eq. 'full') then
	call hdpscolorset (1)
      else if (pscolor(1:4) .eq. 'fore') then
	call hdpscolorset (3)
      else if (pscolor .eq. ' ') then
	call hdpscoloroff
      else
	call hdpscolorset (2)
      end if
      if (ifirst .eq. 0) call clear
      ifirst = 0
c
      return
      end
