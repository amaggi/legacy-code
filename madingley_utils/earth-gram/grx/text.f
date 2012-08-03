c
c*******************************************************************************
c
c    Subroutine text
c
c    Author - Danny Harvey
c
c*******************************************************************************
c
      subroutine text(x, y, angle, iref, string, iclip)
c
c    Subroutine TEXT will cause a text string to be plotted.
c
c    inputs  - x      = reference X-coordinate in user units
c              y      = reference Y-coordinate in user units
c              angle  = polar orientation of character string about point
c                       (x,y) measured in degrees counter-clockwise about
c                       the X-axis
c              iref   = flag which defines the reference point (x,y) relative
c                       to the character string
c                       = 0 - (x,y) is at lower left hand corner of string
c                       = 1 - (x,y) is at center of left hand edge of string
c                       = 2 - (x,y) is at top left hand corner of string
c                       = 3 - (x,y) is at bottom center of string
c                       = 4 - (x,y) is at center of string
c                       = 5 - (x,y) is at top center of string
c                       = 6 - (x,y) is at lower right hand corner of string
c                       = 7 - (x,y) is at center of right hand edge of string
c                       = 8 - (x,y) is at top right hand corner of string
c              string = character string to be plotted
c              iclip  = clip flag
c                       = 0 - clip
c                       .ne. 0 - do not clip
c
      character*(*) string
c
      dimension ix(150),iy(150),ioff(27),npoint(27),iwidth(100)
      dimension xoff(100)
c
      character*1 chold, nchar
      common /chcom/ chold
c
      common /pdim/   xdim,ydim,xlow,ylow,rxdim,rydim,rxlow,rylow,
     1                xbl,ybl,xbh,ybh,xbm,ybm,itran,tangle,
     2                ca,sa,cellht,cellwd,ixtype,iytype
c
      common /partxt/ sinang,cosang,tansln,hite,rat
c
      common /ocflag/  imflag,iocflg,iltp
      integer*4 ltp
c
      data  sinang  / 0. /
      data  cosang  / 1. /
      data  tansln  / 0. /
      data  hite    / 400. /
      data  rat     / 1. /
c
      data  chold  /  '|'  /
c
c    compute reference point in rasters
c
      xx = xmap(x)
      yy = ymap(y)
c
c    clip if desired
c
      if (iclip .ne. 0)  go to 10
      if (xx .lt. xbl) return
      if (xx .gt. xbh) return
      if (yy .lt. ybl) return
      if (yy .gt. ybh) return
   10 if (xx .lt. 0.) return
      if (xx .gt. xbm) return
      if (yy .lt. 0.) return
      if (yy .gt. ybm) return
c
c    delete trailing blanks and determine number of characters in string
c
      j = len(string)
      lstr = 0
      do 20  i = 1,j
      if (string(i:i) .ne. ' ') lstr = i
   20 continue
      if (lstr .eq. 0) return
c
c    set angle of labelling
c
      call chrdir(angle)
c
c    determine overall height and width of string
c
      iiref = iref
      if (iref .lt. 0) iiref = 0
      if (iref .gt. 8) iiref = 0
      iswid = 0
      do 30  i = 1,lstr
      nchar = string(i:i)
      call getsiz(string(i:i), ishite, iwidth(i))
   30 iswid = iswid + iwidth(i)
c
c    compute scaling factors based on string height and ratio
c
      vscale = hite/ishite
      hscale = vscale*rat
      i1 = iiref/3
      i2 = iiref - i1*3
      yoff = float(ishite)*i2/2.
      xxoff = float(iswid)*i1/2.
      xoff(1) = xxoff
      if (lstr .eq. 1)  go to 41
      do 40 i = 2,lstr
   40 xoff(i) = xoff(i-1) - iwidth(i-1)
   41 continue
c
c    single character plot loop starts here
c
      call hdthik (2)
      call hdtext (x, y, angle, iref, string(1:lstr), iclip)
      if (imflag .eq. 1)  call NOPEN
      if (imflag .eq. 1)  call Nlstyl(0)
      do 100  i = 1,lstr
c
c    first get the character strokes
c
      if (string(i:i) .ne. chold)
     1call getchr(string(i:i), nstrks, npoint, ioff, ix, iy)
      chold = string(i:i)
c
c    start of stroke loop
c
      if (nstrks .eq. 0)  go to 111
      do 110  j = 1,nstrks
      np = npoint(j)
      joff = ioff(j)
c
c    start of single stroke plot loop
c
      do 120  k = 1,np
      k1 = joff + k - 1
      xp = ix(k1)
      yp = iy(k1)
c
c    slant point
c
      xp = xp + yp*tansln
c
c    offset point relative to reference point
c
      xp = xp - xoff(i)
      yp = yp - yoff
c
c    scale point to rasters
c
      xp = xp*hscale
      yp = yp*vscale
c
c    rotate point about reference point
c
      xt = xp*cosang - yp*sinang
      yp = xp*sinang + yp*cosang
      xp = xt
c
c    add point to reference point
c
      xp = xp + xx
      yp = yp + yy
c
c    plot point
c
      if (k .eq. 1) call movev(xp, yp, itran)
      call drawv(xp, yp, itran)
c
c    end of single stroke plot loop
c
  120 continue
c
c    end of stroke loop
c
  110 continue
  111 continue
c
c    end of character loop
c
  100 continue
      if (imflag .eq. 1) call Nclose
      if (imflag .eq. 1)  iocflg = 0
      call hdstrk
      call hdstrkl
      call hdthik (0)
      call hdtxtf
c
c    normal exit
c
      return
c
c    end of Subroutine TEXT
c
      end
