      subroutine plot1(n,x,y,roff,iclip)
c
c    routine plot1 will cause one single curve to be plotted AS A SOLID
c    LINE ONLY.
c
c    inputs  - n      = number of points to plot must be greater than one
c			(if point plots or symbol plots are desired use
c			the appropriate routines)
c              x      = plot data horizontal coordinate array
c              y      = plot data vertical coordinate array
c              roff   = the offset in standard raster units from the
c			true curve for this particular curve
c			(used to thicken lines)
c	       iclip  = flag to turn on or off the soft clipping
c			defined by the physical plot dimensions
c			= 0 - clip
c			.ne. 0 - do not clip
c
      common /pdim/ xdim,ydim,xlow,ylow,rxdim,rydim,rxlow,rylow,
     1              xxbl,yybl,xxbh,yybh,xbm,ybm,itran,tangle,
     2              ca,sa,cellht,cellwd,ixtype,iytype
c
      common /pscl/ xmin,xmax,ymin,ymax,xscale,yscale,xrange,yrange
c
      common /filcom/ ifill, jfill
c
      dimension x(*),y(*)
c
      real*8 x8,rx,ry,sn,cs,sno,cso,s,rxn,ryn
c
      k = 1
      do 10  i = 2,n
      if (x(i) .eq. x(k) .and. y(i) .eq. y(k))  go to 10
      k = k + 1
      x(k) = x(i)
      y(k) = y(i)
   10 continue
      n = k
      if (n .lt. 1) return
      iltold = iltp
      if (iclip .ne. 0)  go to 5
      xbh = xxbh
      xbl = xxbl
      ybh = yybh
      ybl = yybl
      go to 6
    5 xbl = 0.
      ybl = 0.
      xbh = xbm
      ybh = ybm
    6 continue
      J = 0
      IFLAG = 0
      ifirst = 1
      xxo = xmap(x(1))
      yyo = ymap(y(1))
      if (n .eq. 1)  go to 50
      xt = xmap(x(1))
      rx = xt
      yt = ymap(y(1))
      ry = yt
      xtn = xmap(x(2))
      rxn = xtn
      ytn = ymap(y(2))
      ryn = ytn
      cs = rxn - rx
      sn = ryn - ry
      s = dsqrt(cs**2 + sn**2)
      cs = cs/s
      sn = sn/s
      cso = cs
      sno = sn
      call offset(roff,rx,ry,sn,cs,sno,cso,xxo,yyo)
   50 IBXO = 0
      IBYO = 0
      IF (XXO .GT. XBH) IBXO = 2
      IF (XXO .LT. XBL) IBXO = 1
      IF (YYO .GT. YBH) IBYO = 2
      IF (YYO .LT. YBL) IBYO = 1
      IF (IBXO+IBYO .NE. 0)  GO TO 90
      ifirst = 1
      J = J + 1
      IFLAG = 1
      call tdraw(xxo,yyo,ifirst,itran)
      if (n .eq. 1) go to 900
   90 DO 100  I = 2,n
      rx = rxn
      ry = ryn
      if (i .eq. n)  go to 91
      xtn = xmap(x(i+1))
      rxn = xtn
      ytn = ymap(y(i+1))
      ryn = ytn
      cs = rxn - rx
      sn = ryn - ry
      s = dsqrt(cs**2 + sn**2)
      cs = cs/s
      sn = sn/s
   91 call offset(roff,rx,ry,sn,cs,sno,cso,xx,yy)
      IBX = 0
      IBY = 0
      IF (XX .GT. XBH)  IBX = 2
      IF (XX .LT. XBL)  IBX = 1
      IF (YY .GT. YBH)  IBY = 2
      IF (YY .LT. YBL)  IBY = 1
      IF (IBX+IBY .EQ. 0 .AND. IFLAG .EQ. 1)  GO TO 195
      IF (IBX+IBY .EQ. 0 .AND. IGRAF .LT. 0)  GO TO 195
      IF (IBX+IBY .NE. 0 .AND. IFLAG .EQ. 0)  GO TO 198
      IF (IBX+IBY .NE. 0 .AND. IGRAF .LT. 0)  GO TO 198
      IF (IFLAG .EQ. 1)  GO TO 150
      IBX = IBXO
      IBY = IBYO
  150 IF (IBX .NE. 2)  GO TO 155
      XINT = XBH
      GO TO 160
  155 IF (IBX .NE. 1)  GO TO 170
      XINT = XBL
  160 YINT = YYO + (XINT-XXO)*(YY-YYO)/(XX-XXO)
      IF (YINT .GT. YBH .OR. YINT .LT. YBL)  GO TO 170
      GO TO 190
  170 IF (IBY .NE. 2)  GO TO 175
      YINT2 = YBH
      GO TO 180
  175 IF (IBY .NE. 1)  GO TO 190
      YINT2 = YBL
  180 XINT2 = XXO + (YINT2-YYO)*(XX-XXO)/(YY-YYO)
      IF (XINT2 .GT. XBH .OR. XINT2 .LT. XBL)  GO TO 190
      XINT = XINT2
      YINT = YINT2
  190 J = J + 1
      if (j .eq. 1 .and. jfill .eq. 0) ifirst = 1
      call tdraw(xint,yint,ifirst,itran)
      IF (IFLAG .NE. 1)  GO TO 195
      J = 0
      GO TO 198
  195 J = J + 1
      IFLAG = 1
      if (j .eq. 1 .and. jfill .eq. 0) ifirst = 1
      call tdraw(xx,yy,ifirst,itran)
      GO TO 199
  198 IFLAG = 0
      IBXO = IBX
      IBYO = IBY
  199 XXO = XX
  100 YYO = YY
c
  900 continue
      return
      end
c
      subroutine offset(roff,rx,ry,sn,cs,sno,cso,xx,yy)
c
c    routine offset computes the offseted point in standard rasters
c    for line thickening
c
      implicit real*8 (a-h,o-z)
      real*4 roff,xx,yy
c
      if (roff .eq. 0.)  go to 300
      rxn = rx - roff*sn
      ryn = ry + roff*cs
      del = -sn*cso + cs*sno
      if (del .eq. 0.)  go to 500
      rxo = rx - roff*sno
      ryo = ry + roff*cso
      x = rxn*sn-ryn*cs
      y = rxo*sno-ryo*cso
      x2 = (-cso*x+cs*y)/del
      y2 = (-sno*x+sn*y)/del
      sno = sn
      cso = cs
      xx = x2
      yy = y2
c
      return
  500 xx = rxn
      yy = ryn
      sno = sn
      cso = cs
      return
c
  300 xx = rx
      yy = ry
      return
c
      end
