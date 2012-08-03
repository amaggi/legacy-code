c
c*******************************************************************************
c
c    Subroutine ncontour
c
c    Author - Danny Harvey
c
c*******************************************************************************
c
      subroutine ncontour (nx, nxmax, ny, x, y, z, ctype, zmin, zmax,
     1                     dz)
c
      integer              nx
      integer                  nxmax
      integer                      ny
      real*4                           x(nx)
      real*4                              y(ny)
      real*4                                 z(nxmax,ny)
      character*(*)                             ctype
      real*4                                           zmin
      real*4                                                 zmax
      real*4               dz
c
c    Subroutine ncontour will render line and/or color contours for displaying
c    z = z(x,y) in "map" view. Limits and scales for the x-y axes must be
c    previously established with calls to setdim, setscl, or axis. z(x,y)
c    must be specified as a regular rectangular grid of points.
c
c    Inputs  -	nx	= The number of x-coordinates in the grid.
c    		ny	= The number of y-coordinates in the grid.
c    		x(nx)	= The x-coordinates. Repeated values are not allowed
c			  and the points must either increase or decrease
c			  monotonically.
c    		y(ny)	= The y-coordinates. Repeated values are not allowed
c			  and the points must either increase or decrease
c			  monotonically.
c		z(nxmax,ny)	
c			= The z(x,y) values for each x-y grid point. These
c			  must be given in x order first.
c		ctype	= The contour type flag.
c			  = 'lines'  - Make line contours only.
c			  = 'colors' - Make color contours only.
c			  = 'both'   - Make line contours on top of color
c				       contours.
c		zmin	= The minimum z-value for contour clipping.
c		zmax	= The maximum z-value for contour clipping. If 
c			  zmax = zmin, then no clipping will take place and
c			  internal values of zmin and zmax will be computed
c			  automatically from the range of the data.
c		dz    	= The z contour line increment. If dz <= 0.0, then
c			  a "nice" increment will be chosen.
c
c    Output -	dz	= Chosen "nice" increment in input dz <= 0.0
c
      parameter  (NXYM  = 20000)
c
      real*4 xout(NXYM), yout(NXYM)
c
      common /legend/ nlevel, zlevel, dzlevl
c
      if (nx .lt. 2) return
      if (ny .lt. 2) return
      if (x(nx) .gt. x(1)) then
	do 5  i = 2, nx
	  if (x(i) .le. x(i-1)) return
    5   continue
      else
	do 6  i = 2, nx
	  if (x(i) .ge. x(i-1)) return
    6   continue
      end if
      if (y(ny) .gt. y(1)) then
	do 7  i = 2, ny
	  if (y(i) .le. y(i-1)) return
    7   continue
      else
	do 8  i = 2, ny
	  if (y(i) .ge. y(i-1)) return
    8   continue
      end if
      if (zmin .eq. zmax) then
	zzmin = z(1,1)
	zzmax = z(1,1)
	do 10  j = 1, ny
	do 10  i = 1, nx
	  if (z(i,j) .lt. zzmin) zzmin = z(i,j)
	  if (z(i,j) .gt. zzmax) zzmax = z(i,j)
   10   continue
      else
	if (zmin .lt. zmax) then
	  zzmin = zmin
	  zzmax = zmax
	else
	  zzmin = zmax
	  zzmax = zmin
	end if
      end if
      if (zzmax .eq. zzmin) return
      if (dz .le. 0.0) then
	ddz = (zzmax-zzmin)/20.0
	call dnice (ddz, ds, db)
	ddz = ds
      else
	ddz = dz
      end if
      if (ddz .eq. 0.0) return
      if (zzmin .lt. 0.0) then
        nz = zzmin / ddz
      else
        nz = zzmin / ddz + 1.0
      end if
      zstart = nz*ddz
      nz = (zzmax - zstart) / ddz + 1.0
      if (nz .lt. 1) return
      nptsm = NXYM
c     j1 = 1
c     i1 = 1
c     do 80  j = 1, ny
c     do 80  i = 1, nx
c       z(i1, j1) = z(i,j)
c       i1 = i1 + 1
c       if (i1 .gt. nxmax) then
c         i1 = 1
c         j1 = j1 + 1
c       end if
c  80 continue
      nlevel = nnz
      zlevel = zstart
      dzlevl = ddz
      if (ctype .eq. 'colors' .or. ctype .eq. 'both') then
	nnz = nz + 1
        if (nnz .gt. 100) nnz = 100
	if (zmax .ne. 0.0) zzmax = zmax
	call ncrimg (nnz, image)
	call ninimg (image, nx, nxmax, ny, x, y, z, zstart, ddz)
	call npnimg (image)
	call hdfrimg (image)
      end if
      if (ctype .eq. 'lines' .or. ctype .eq. 'both') then
        call contrd (nx,x,'LIN', ny, y, 'LIN', nxmax, z, 'LIN', zstart,
     1            ddz, nz, nptsm, xout, yout)
      end if
c
c    Normal exit.
c
      dz = ddz
      return
      end
c
c*******************************************************************************
c
c    Subroutine contr
c
c    Author - Danny Harvey
c
c*******************************************************************************
c
      subroutine contr (nx, x, xtype, ny, y, ytype,
     +                  nxmax, z, ztype, zstart, dz, nz, nptsm,
     +                  xout, yout, 
     +                  izxy, zs, zxmin, zxmax, zymin, zymax)
c
c    Subroutine contr will compute the X-Y coordinates of constant
c    Z value contours for some Z = Z(X,Y).  Z(X,Y) must be sampled
c    along some rectangular grid of X-Y values although the grid
c    need not have equal spacing everywhere.
c
c    inputs  - nx     = number of different X values in the grid.
c              x(nx)  = X sampling values.  These must be ordered
c                       but can be either increasing or decreasing.
c              xtype  = CHARACTER*3 flag which indicates whether to
c                       return linear ('LIN') or logarithmic ('LOG')
c                       contour values of X in xout
c              ny     = number of different Y values in the grid.
c              y(ny)  = Y sampling values.  These must be ordered
c                       but can be either increasing or decreasing.
c              ytype  = CHARACTER*3 flag which indicates whether to
c                       return linear ('LIN') or logarithmic ('LOG')
c                       contour values of Y in yout
c              z(nx*ny) = sampled Z values for each X-Y sampling point.
c                         This array is indexed by X values first.
c              ztype  = CHARACTER*3 flag which indicates whether to
c                       make the contour spacing linear ('LIN') or
c                       logarithmic ('LOG')
c              zstart = first contour Z value.
c              dz     = contour Z value increment.
c              nz     = number of contours to be computed.
c              nptsm  = maximum number of values in the xout and yout
c                       arrays.
c              xout(nptsm)
c		      = Array used for storing X-contour plot points.
c              yout(nptsm)
c		      = Array used for storing Y-contour plot points.
c
      parameter (MAX  = 200000)
      parameter (MAX2 = 500000)
c
      integer*4  nx, ny, nz
      integer*4  nptsm
      real*4     x(*), y(*), z(*), zstart, dz
      real*4     xout(*), yout(*)
      character*(*)  xtype, ytype, ztype
c
      common /clims/ xmin, xmax, ymin, ymax
      common /ctype/ xt, yt
      character*3    xt, yt
c
c     integer*2  izxy(MAX2)
c     real*4     zs(MAX), zxmin(MAX), zxmax(MAX),
c    +           zymin(MAX), zymax(MAX)
      integer*2  izxy(*)
      real*4     zs(*), zxmin(*), zxmax(*),
     +           zymin(*), zymax(*)
      logical itest
c
      xt = xtype
      yt = ytype
c
c    look at sizes
c
c     if (nx .lt. 2 .or. nx .gt. MAX) then
c       write (6, '(a)') 'CONTR: wrong # of X values'
c       return
c     end if
c     if (ny .lt. 2 .or. ny .gt. MAX) then
c       write (6, '(a)') 'CONTR: wrong # of Y values'
c       return
c     end if
      ixlim = nx
      iylim = ny
      if (nz .lt. 1) return
c     if (nz .gt. MAX) then
c       write (6, '(a)') 'CONTR: too many contours'
c       return
c     end if
      nxy = nx
      nxy = nxy*ny
c     if (nxy .gt. MAX2/4) then
c       write (6, '(a)') 'CONTR: too many Z values'
c       return
c     end if
c
c    determine min and max X,Y,Z values
c
      if (x(nx) .eq. x(1)) then
        write (6, '(a)') 'CONTR: repeated X values'
        return
      end if
      if (x(nx) .gt. x(1)) then
        do 5000 i = 2, nx
          if (x(i) .eq. x(i-1)) then
            write (6, '(a)') 'CONTR: repeated X values'
            return
          end if
          if (x(i) .lt. x(i-1)) then
            write (6, '(a)') 'CONTR: X array out of order'
            return
          end if
 5000   continue
        xmin = x(1)
        xmax = x(nx)
      else
        do 5010 i = 2, nx
          if (x(i) .eq. x(i-1)) then
            write (6, '(a)') 'CONTR: repeated X values'
            return
          end if
          if (x(i) .gt. x(i-1)) then
            write (6, '(a)') 'CONTR: X array out of order'
            return
          end if
 5010   continue
        xmin = x(nx)
        xmax = x(1)
      end if
c
      if (y(ny) .eq. y(1)) then
        write (6, '(a)') 'CONTR: repeated Y values'
        return
      end if
      if (y(ny) .gt. y(1)) then
        do 5020 i = 2, ny
          if (y(i) .eq. y(i-1)) then
            write (6, '(a)') 'CONTR: repeated Y values'
            return
          end if
          if (y(i) .lt. y(i-1)) then
            write (6, '(a)') 'CONTR: Y array out of order'
            return
          end if
 5020   continue
        ymin = y(1)
        ymax = y(ny)
      else
        do 5030 i = 2, ny
          if (y(i) .eq. y(i-1)) then
            write (6, '(a)') 'CONTR: repeated Y values'
            return
          end if
          if (y(i) .gt. y(i-1)) then
            write (6, '(a)') 'CONTR: Y array out of order'
            return
          end if
 5030   continue
        ymin = y(ny)
        ymax = y(1)
      end if
c
      zmin = 2.e30
      zmax = -2.e30
      do 5040 i = 1, nx
        zxmin(i) = 2.e30
        zxmax(i) = -2.e30
        do 5050 j = 1, ny
	  k = i + (j-1)*nxmax
          zz = z(k)
          if (zz .gt. zmax) zmax = zz
          if (zz .lt. zmin) zmin = zz
          if (zz .gt. zxmax(i)) zxmax(i) = zz
          if (zz .lt. zxmin(i)) zxmin(i) = zz
 5050   continue
 5040 continue
      do 5060 j = 1, ny
        zymin(j) = 2.e30
        zymax(j) = -2.e30
        do 5070 i = 1, nx
	  k = i + (j-1)*nxmax
          zz = z(k)
          if (zz .gt. zymax(j)) zymax(j) = zz
          if (zz .lt. zymin(j)) zymin(j) = zz
 5070   continue
 5060 continue
c
c    determine contour Z values
c
      if (ztype(2:2) .eq. 'o' .or. ztype(2:2) .eq. 'O') then
        do 5080 i = 1, nz
          zs(i) = zstart*(dz**(i-1))
 5080   continue
      else
        do 5090 i = 1, nz
          zs(i) = zstart + (i-1)*dz
 5090   continue
      end if
c     print *,' xmin,xmax = ',xmin,xmax
c     print *,' ymin,ymax = ',ymin,ymax
c     print *,' zmin,zmax = ',zmin,zmax
c
c    main processing loop indexed on contours
c
      do 5100 i = 1, nz
        zzs = zs(i)
c
c    check contour value against overall Z data range
c
        if (zzs .ge. zmin .and. zzs .le. zmax) then
c
c    clear out the izxy array
c
          do 5110 j = 1, 4*nxy
            izxy(j) = 0
 5110     continue
c
c    first look along bottom and top boundaries
c
          is = 1
          do 5120 j = 0, iylim-1, iylim-1
            j1 = j + 1
c
c    check contour value against Z data range for this Y
c
            if (zzs .ge. zymin(j1) .and. zzs .le. zymax(j1)) then
c
c    range over X values
c
              do 5130 k = 0, ixlim-2
                k1 = k + 1
                i1 = k1 + j*nxmax
                i2 = i1 + 1
                i3 = 2*k1 + 2*j*(2*ixlim-1)
                if (mod(izxy(i3),2) .eq. 0 .and.
     +              itest(z(i1),z(i2),zzs)) then
c                 print *,' bottom/top j,k = ',j,k
c                 print *,' x,y,z1 = ',x(k1),y(j1),z(i1)
c                 print *,' x,y,z2 = ',x(k1+1),y(j1),z(i2)
                  izxy(i3) = izxy(i3) + 1
                  xpnt = (x(k1)*(z(i2)-zzs)-x(k1+1)*(z(i1)-zzs))/
     +                   (z(i2)-z(i1))
                  ypnt = y(j1)
                  call cscale (xpnt, ypnt)
                  call ctrace (j, k, is, 0, xpnt, ypnt, x, y, z,
     +                         zzs, izxy, ixlim-1, iylim-1,
     +                         nptsm, zout, npts,
     +                         xout, yout)
                end if
 5130         continue
            end if
            is = -1
 5120     continue
c
c    look along left and right boundaries
c
          is = 1
          do 5140 k = 0, ixlim-1, ixlim-1
            k1 = k + 1
c
c    check contour value against Z data range for this X
c
            if (zzs .ge. zxmin(k1) .and. zzs .le. zxmax(k1)) then
c
c    range over Y values
c
              do 5150 j = 0, iylim-2
                j1 = j + 1
                i1 = k1 + j*nxmax
                i2 = i1 + nxmax
                i3 = 2*k + 1 + (2*j+1)*(2*ixlim-1)
                if (mod(izxy(i3),2) .eq. 0 .and.
     +              itest(z(i1),z(i2),zzs)) then
c                 print *,' left/right j,k = ',j,k
c                 print *,' x,y,z1 = ',x(k1),y(j1),z(i1)
c                 print *,' x,y,z2 = ',x(k1),y(j1+1),z(i2)
                  izxy(i3) = izxy(i3) + 1
                  ypnt = (y(j1)*(z(i2)-zzs)-y(j1+1)*(z(i1)-zzs))/
     +                   (z(i2)-z(i1))
                  xpnt = x(k1)
                  call cscale (xpnt, ypnt)
                  call ctrace (j, k, 0, is, xpnt, ypnt, x, y, z,
     +                         zzs, izxy, ixlim-1, iylim-1,
     +                         nptsm, zout, npts,
     +                         xout, yout)
                end if
 5150         continue
            end if
            is = -1
 5140     continue
c
c    now look along inside Y values
c
          do 5160 j = 1, iylim-2
            j1 = j + 1
c
c    check contour value against Z data range for this Y
c
            if (zzs .ge. zymin(j1) .and. zzs .le. zymax(j1)) then
c
c    range over X values
c
              do 5170 k = 0, ixlim-2
                k1 = k + 1
                i1 = k1 + j*nxmax
                i2 = i1 + 1
                i3 = 2*k1 + 2*j*(2*ixlim-1)
                if (mod(izxy(i3),2) .eq. 0 .and.
     +              itest(z(i1),z(i2),zzs)) then
c                 print *,' Y j,k = ',j,k
c                 print *,' x,y,z1 = ',x(k1),y(j1),z(i1)
c                 print *,' x,y,z2 = ',x(k1+1),y(j1),z(i2)
                  izxy(i3) = izxy(i3) + 1
                  xpnt = (x(k1)*(z(i2)-zzs)-x(k1+1)*(z(i1)-zzs))/
     +                   (z(i2)-z(i1))
                  ypnt = y(j1)
                  call cscale (xpnt, ypnt)
                  call ctrace (j, k, 1, 0, xpnt, ypnt, x, y, z,
     +                         zzs, izxy, ixlim-1, iylim-1,
     +                         nptsm, zout, npts,
     +                         xout, yout)
                end if
 5170         continue
            end if
 5160     continue
c
c    finally look along inside X values
c
          do 5180 k = 1, ixlim-2
            k1 = k + 1
c
c    check contour value against Z data range for this X
c
            if (zzs .ge. zxmin(k1) .and. zzs .le. zxmax(k1)) then
c
c    range over Y values
c
              do 5190 j = 0, iylim-2
                j1 = j + 1
                i1 = k1 + j*nxmax
                i2 = i1 + nxmax
                i3 = 2*k + 1 + (2*j+1)*(2*ixlim-1)
                if (mod(izxy(i3),2) .eq. 0 .and.
     +              itest(z(i1),z(i2),zzs)) then
c                 print *,' X j,k = ',j,k
c                 print *,' x,y,z1 = ',x(k1),y(j1),z(i1)
c                 print *,' x,y,z2 = ',x(k1),y(j1+1),z(i2)
                  izxy(i3) = izxy(i3) + 1
                  ypnt = (y(j1)*(z(i2)-zzs)-y(j1+1)*(z(i1)-zzs))/
     +                   (z(i2)-z(i1))
                  xpnt = x(k1)
                  call cscale (xpnt, ypnt)
                  call ctrace (j, k, 0, 1, xpnt, ypnt, x, y, z,
     +                         zzs, izxy, ixlim-1, iylim-1,
     +                         npts, zout, npts,
     +                         xout, yout)
                end if
 5190         continue
            end if
 5180     continue
        end if
c       write (6, 200) i, zzs, nsegs, nptr(nsegs)
  200   format (' contour # ',I3,', value = ',E12.6,
     +          ', nsegs = ',I4,', pointer at ',I5)
c
c    end of contour loop
c
 5100 continue
c
c    normal exit
c
      return
c
      end
      subroutine ctrace (ji, ki, jstepi, kstepi, xpnt, ypnt,
     +                   x, y, z, zs, izxy, ixlim, iylim,
     +                   nptsm, zout, npts,
     +                   xout, yout)
c
      integer*4  i, ji, ki, jstepi, kstepi, ixlim, iylim
      integer*2  izxy(*)
      real*4     xpnt, ypnt, x(*), y(*), z(*), zs
c
      logical itest
c
      integer*4  npts
      integer*4  nptsm
      real*4     xout(*), yout(*), zout
c
      npts = 1
      xout(npts) = xpnt
      yout(npts) = ypnt
      zout = zs
      j = ji
      j1 = j + 1
      k = ki    
      k1 = k + 1
      jstep = jstepi
      kstep = kstepi
      if (kstep .eq. 0) then
        kstep = 1
  100   j2 = j + jstep
        if (j2 .gt. iylim .or. j2 .lt. 0) go to 900
        i1 = k + j2*(ixlim + 1) + 1
        i2 = i1 + kstep
        i3 = 2*k + kstep + 2*j2*(2*ixlim + 1) + 1
        if (itest(z(i1),z(i2),zs)) then
          xpnt = (x(k1)*(z(i2)-zs)-x(k1+kstep)*(z(i1)-zs))/
     +           (z(i2)-z(i1))
          j = j + jstep
          j1 = j + 1
          ypnt = y(j1)
          call cscale (xpnt, ypnt)
          npts = npts + 1
          if (npts .gt. nptsm) then
            write (6, '(a)') 'CTRACE: attempt to exceed max # of'
     +                   // ' output points'
	    npts = npts - 1
            go to 900
          end if
          xout(npts) = xpnt
          yout(npts) = ypnt
          if (mod (izxy(i3), 2) .ne. 0) go to 900
          izxy(i3) = izxy(i3) + 1
          go to 100
        end if
        i1 = i2 - jstep*(ixlim+1)
        if (itest(z(i1),z(i2),zs)) then
          k = k + kstep
          k1 = k + 1
        else
          i2 = i2 - kstep
          i1 = i2 - jstep*(ixlim+1)
          if (itest(z(i1),z(i2),zs)) then
            kstep = -kstep
          else
            go to 900
          end if
        end if
        xpnt = x(k1)
        ypnt = (y(j1)*(z(i2)-zs)-y(j1+jstep)*(z(i1)-zs))/
     +         (z(i2)-z(i1))
        call cscale (xpnt, ypnt)
        npts = npts + 1
        if (npts .gt. nptsm) then
          write (6, '(a)') 'CTRACE: attempt to exceed max # of'
     +                   // ' output points'
	    npts = npts - 1
          go to 900
        end if
        xout(npts) = xpnt
        yout(npts) = ypnt
        i3 = i3 + kstep - jstep*(2*ixlim+1)
        if (mod (izxy(i3), 2) .ne. 0) go to 900
        izxy(i3) = izxy(i3) + 1
        go to 200
      else
        jstep = 1
  200   k2 = k + kstep
        if (k2 .gt. ixlim .or. k2 .lt. 0) go to 900
        i1 = k2 + j*(ixlim + 1) + 1
        i2 = i1 + jstep*(ixlim+1)
        i3 = 2*k2 + (2*j+jstep)*(2*ixlim + 1) + 1
        if (itest(z(i1),z(i2),zs)) then
          ypnt = (y(j1)*(z(i2)-zs)-y(j1+jstep)*(z(i1)-zs))/
     +           (z(i2)-z(i1))
          k = k + kstep
          k1 = k + 1
          xpnt = x(k1)
          call cscale (xpnt, ypnt)
          npts = npts + 1
          if (npts .gt. nptsm) then
            write (6, '(a)') 'CTRACE: attempt to exceed max # of'
     +                   // ' output points'
	    npts = npts - 1
            go to 900
          end if
          xout(npts) = xpnt
          yout(npts) = ypnt
          if (mod (izxy(i3), 2) .ne. 0) go to 900
          izxy(i3) = izxy(i3) + 1
          go to 200
        end if
        i1 = i2 - kstep
        if (itest(z(i1),z(i2),zs)) then
          j = j + jstep
          j1 = j + 1
        else
          i2 = i2 - jstep*(ixlim+1)
          i1 = i2 - kstep
          if (itest(z(i1),z(i2),zs)) then
            jstep = -jstep
          else
            go to 900
          end if
        end if
        ypnt = y(j1)
        xpnt = (x(k1)*(z(i2)-zs)-x(k1+kstep)*(z(i1)-zs))/
     +         (z(i2)-z(i1))
        call cscale (xpnt, ypnt)
        npts = npts + 1
        if (npts .gt. nptsm) then
          write (6, '(a)') 'CTRACE: attempt to exceed max # of'
     +                   // ' output points'
	  npts = npts - 1
          go to 900
        end if
        xout(npts) = xpnt
        yout(npts) = ypnt
        i3 = i3 - kstep + jstep*(2*ixlim+1)
        if (mod (izxy(i3), 2) .ne. 0) go to 900
        izxy(i3) = izxy(i3) + 1
        go to 100
      end if
c
  900 continue
      thick = 0.0
      call nplot (npts, xout, yout, 0, 0, thick, 0, ' ')
      return
c
      end
      subroutine cscale (xpnt, ypnt)
c
      real*4  xpnt, ypnt, xmin, xmax, ymin, ymax
      character*3 xtype, ytype
c
      common /clims/ xmin, xmax, ymin, ymax
      common /ctype/ xtype, ytype
c
c     print *,' xpnt,ypnt = ',xpnt,ypnt
      if ((xpnt .gt. xmax .and. xpnt .gt. xmin) .or.
     +    (xpnt .lt. xmax .and. xpnt .lt. xmin)) then
        write (6, '(a)') 'CSCALE: xpnt out of range'
        stop
      end if
      if (xtype(2:2) .eq. 'o' .or. xtype(2:2) .eq. 'O') then
        xpnt = alog10(xpnt)
      end if
c
      if ((ypnt .gt. ymax .and. ypnt .gt. ymin) .or.
     +    (ypnt .lt. ymax .and. ypnt .lt. ymin)) then
        write (6, '(a)') 'CSCALE: ypnt out of range'
        stop
      end if
      if (ytype(2:2) .eq. 'o' .or. ytype(2:2) .eq. 'O') then
        ypnt = alog10(ypnt)
      end if
c
      return
c
      end
      logical function itest (a, b, c)
c
      real*4  a, b, c
c
      if (a .ge. 1.e29 .or. b .ge. 1.e29) then
        itest = .false.
      else
        x1 = a - c
        if (x1 .ge. 0.) then
          x1 = 1.
        else
          x1 = -1.
        end if
        x2 = b - c
        if (x2 .ge. 0.) then
          x2 = 1.
        else
          x2 = -1.
        end if
        if (x1*x2 .lt. 0.) then
          itest = .true.
        else
          itest = .false.
        end if
      end if
c
      return
c
      end
