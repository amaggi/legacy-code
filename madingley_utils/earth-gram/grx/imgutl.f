	subroutine resamp( rdat, image, w, h, ww, hh,
     +		xoff, yoff, depth, hd_depth )

	implicit none

	integer w, h
	integer ww, hh
	integer xoff, yoff
	integer rdat(w,h)
	character image(ww,hh)
	integer depth, hd_depth

	integer i, j, k, l
	integer ii, jj
	integer sum
	real*4 factor

	if ( hd_depth .eq. 1 ) then
	  factor = 17.0 / ( 16.0 * depth )
	else
	  factor = 49.0 / ( 16.0 * depth )
	endif

	do 100 i = 1, w
	  ii = xoff + 1 + 4 * ( i - 1 )
	  do 200 j = 1, h
	    jj = xoff + 1 + 4 * ( j - 1 )
	    sum = 0
	    do 300 k = ii, ii + 3
	      do 400 l = jj, jj + 3
		sum = sum + ichar( image(k,l) )
400	      continue
300	    continue
	    rdat(i,j) = 1 + int( 0.499 + factor * sum )
200	  continue
100	continue

	end



	subroutine setmono( x, out, bpl, w, h, col, rv )

	implicit	none

	integer		bpl
	integer		w, h
	integer		x(w,h)
c	byte		out(bpl,4*h)
	character	out(bpl,4*h)
	integer		col(4,17)
	integer		rv

	integer		cx

	integer		i, j, k, l

	do 10 i = 1, w, 2
	  k = ( i + 1 ) / 2
	  do 20 j = 1, 4*h-3, 4
	    l = ( j + 3 ) / 4

	    if ( rv .eq. 0 ) then
	      cx = x(i,l)
	    else
	      cx = 18 - x(i,l)
	    endif
	    cx = max( cx, 1 )
	    cx = min( cx, 17 )

c
c	using color conversion table to set up image data "out"
c
c	    out(k,j)   = col( 1, cx ) * 16
c	    out(k,j+1) = col( 2, cx ) * 16
c	    out(k,j+2) = col( 3, cx ) * 16
c	    out(k,j+3) = col( 4, cx ) * 16
	    call cset( out(k,j)   , col( 1, cx ) * 16 )
	    call cset( out(k,j+1) , col( 2, cx ) * 16 )
	    call cset( out(k,j+2) , col( 3, cx ) * 16 )
	    call cset( out(k,j+3) , col( 4, cx ) * 16 )


20	  continue
10	continue

	do 30 i = 2, w, 2
	  k = i / 2
	  do 40 j = 1, 4*h-3, 4
	    l = ( j + 3 ) / 4

	    if ( rv .eq. 0 ) then
	      cx = x(i,l)
	    else
	      cx = 18 - x(i,l)
	    endif
	    cx = max( cx, 1 )
	    cx = min( cx, 17 )

c	    out(k,j)   = out(k,j)   + col( 1, cx )
c	    out(k,j+1) = out(k,j+1) + col( 2, cx )
c	    out(k,j+2) = out(k,j+2) + col( 3, cx )
c	    out(k,j+3) = out(k,j+3) + col( 4, cx )
	    call cadd( out(k,j)   , col( 1, cx ) )
	    call cadd( out(k,j+1) , col( 2, cx ) )
	    call cadd( out(k,j+2) , col( 3, cx ) )
	    call cadd( out(k,j+3) , col( 4, cx ) )

40	  continue
30	continue

	return

	end



	subroutine setgrey( x, out, bpl, w, h, col, rv )

	implicit	none

	integer		bpl
	integer		w, h
	integer		x(w,h)
c	byte		out(bpl,4*h)
	character	out(bpl,4*h)
	integer		col(4,17)
	integer		rv

	integer		cx, cy

	integer		i, j, k, l

	do 10 i = 1, w
	  k = i
	  do 20 j = 1, 4*h-3, 4
	    l = ( j + 3 ) / 4

	    if ( rv .eq. 0 ) then
	      cx = x(i,l)
	    else
	      cx = 50 - x(i,l)
	    endif
	    cx = max( cx, 1 )
	    cx = min( cx, 49 )
	    cy = ( cx - 1 ) / 16
	    cx = cx - 16 * cy
	    cy = cy * 85

c
c	using color conversion table to set up image data "out"
c
c	    out(k,j)   = col( 1, cx ) + cy
c	    out(k,j+1) = col( 2, cx ) + cy
c	    out(k,j+2) = col( 3, cx ) + cy
c	    out(k,j+3) = col( 4, cx ) + cy
	    call cset( out(k,j)   , col( 1, cx ) + cy )
	    call cset( out(k,j+1) , col( 2, cx ) + cy )
	    call cset( out(k,j+2) , col( 3, cx ) + cy )
	    call cset( out(k,j+3) , col( 4, cx ) + cy )

20	  continue
10	continue

	return
	end



	subroutine setmncol( gray, fg )

c
c	setmncol
c
c	Set the z to color conversion table for the case of a monochrome
c	display.
c
c	For a monochrome screen to get more than just two color shades, must
c	use a dithering technique.  Setmncol uses an n = 4 dither and packs
c	the results such that the pixels values can be specified using an
c	BYTE declared array.
c
c	Here col(1,i) is the top row of the dither matrix for an intensity of i.
c	E.g. :  for and intensity of 7 the dither matrix would look like:
c
c	1 0 1 0		col(1,7) = binary number 1010 = decimal number 10
c	0 1 0 1
c	1 0 1 0
c	0 0 0 1		col(4,7) = binary number 0001 = decimal number 1
c
c	All the values for col(i,j) have been pre-computed and are therefore
c	exactly specified below.
c

c	implicit	none

	integer		gray(4,17)
	integer		fg

	integer		col(4,17)
	integer		scale(17)

	real*4		r, g, b, gr

	integer		i

	integer		ncolors


	col(1,1)  = 0
	col(1,2)  = 8
	col(1,3)  = 8
	col(1,4)  = 10
	col(1,5)  = 10
	col(1,6)  = 10
	col(1,7)  = 10
	col(1,8)  = 10
	col(1,9)  = 10
	col(1,10) = 14
	col(1,11) = 14
	col(1,12) = 15
	col(1,13) = 15
	col(1,14) = 15
	col(1,15) = 15
	col(1,16) = 15
	col(1,17) = 15

	col(2,1)  = 0
	col(2,2)  = 0
	col(2,3)  = 0
	col(2,4)  = 0
	col(2,5)  = 0
	col(2,6)  = 4
	col(2,7)  = 4
	col(2,8)  = 5
	col(2,9)  = 5
	col(2,10) = 5
	col(2,11) = 5
	col(2,12) = 5
	col(2,13) = 5
	col(2,14) = 13
	col(2,15) = 13
	col(2,16) = 15
	col(2,17) = 15

	col(3,1)  = 0
	col(3,2)  = 0
	col(3,3)  = 2
	col(3,4)  = 2
	col(3,5)  = 10
	col(3,6)  = 10
	col(3,7)  = 10
	col(3,8)  = 10
	col(3,9)  = 10
	col(3,10) = 10
	col(3,11) = 11
	col(3,12) = 11
	col(3,13) = 15
	col(3,14) = 15
	col(3,15) = 15
	col(3,16) = 15
	col(3,17) = 15

	col(4,1)  = 0
	col(4,2)  = 0
	col(4,3)  = 0
	col(4,4)  = 0
	col(4,5)  = 0
	col(4,6)  = 0
	col(4,7)  = 1
	col(4,8)  = 1
	col(4,9)  = 5
	col(4,10) = 5
	col(4,11) = 5
	col(4,12) = 5
	col(4,13) = 5
	col(4,14) = 5
	col(4,15) = 7
	col(4,16) = 7
	col(4,17) = 15

	open( 1, file = 'xg.mono', status = 'old', err = 100 )

	read(1,*,err=90,end=90) ncolors

	if ( ncolors .ne. 17 ) then
	  call cewrite(
     +'contplotd:  setmncol:  ncolors not equal 17, using default')
	  close(1)
	  goto 100
	endif

	do 10 i = 1, 17
	  read(1,*,err=90,end=90) r, g, b
	  gr = ( r + g + b ) / 3.0
	  scale(i) = int( gr * 16.9999 + 1.0 )
	  if ( scale(i) .lt. 1 ) scale(i) = 1
	  if ( scale(i) .gt. 17 ) scale(i) = 17
	  if ( fg .ne. 1 ) scale(i) = 18 - scale(i)
10	continue

	close(1)

110	do 20 i = 1, 17
	  gray(1,i) = col(1,scale(i))
	  gray(2,i) = col(2,scale(i))
	  gray(3,i) = col(3,scale(i))
	  gray(4,i) = col(4,scale(i))
20	continue

	return

90	continue
	call cewrite(
     +'contplotd:  setmncol:  error reading xg.mono using deefault')
	close(1)

100	continue
	do 30 i = 1, 17
	  scale(i) = i
	  if ( fg .ne. 1 ) scale(i) = 18 - scale(i)
30	continue
	goto 110

	end



	subroutine setgrcol( col )

c
c	setgrcol
c
c	Set the z to color conversion table for the case of a monochrome
c	display.
c
c	For a monochrome screen to get more than just two color shades, must
c	use a dithering technique.  Setgrcol uses an n = 4 dither and packs
c	the results such that the pixels values can be specified using an
c	BYTE declared array.
c

c	implicit	none

	integer		col(4,17)

	col(1,1)  = 0
	col(1,2)  = 64
	col(1,3)  = 64
	col(1,4)  = 68
	col(1,5)  = 68
	col(1,6)  = 68
	col(1,7)  = 68
	col(1,8)  = 68
	col(1,9)  = 68
	col(1,10) = 84
	col(1,11) = 84
	col(1,12) = 85
	col(1,13) = 85
	col(1,14) = 85
	col(1,15) = 85
	col(1,16) = 85
	col(1,17) = 85

	col(2,1)  = 0
	col(2,2)  = 0
	col(2,3)  = 0
	col(2,4)  = 0
	col(2,5)  = 0
	col(2,6)  = 16
	col(2,7)  = 16
	col(2,8)  = 17
	col(2,9)  = 17
	col(2,10) = 17
	col(2,11) = 17
	col(2,12) = 17
	col(2,13) = 17
	col(2,14) = 81
	col(2,15) = 81
	col(2,16) = 85
	col(2,17) = 85

	col(3,1)  = 0
	col(3,2)  = 0
	col(3,3)  = 4
	col(3,4)  = 4
	col(3,5)  = 68
	col(3,6)  = 68
	col(3,7)  = 68
	col(3,8)  = 68
	col(3,9)  = 68
	col(3,10) = 68
	col(3,11) = 69
	col(3,12) = 69
	col(3,13) = 85
	col(3,14) = 85
	col(3,15) = 85
	col(3,16) = 85
	col(3,17) = 85

	col(4,1)  = 0
	col(4,2)  = 0
	col(4,3)  = 0
	col(4,4)  = 0
	col(4,5)  = 0
	col(4,6)  = 0
	col(4,7)  = 1
	col(4,8)  = 1
	col(4,9)  = 17
	col(4,10) = 17
	col(4,11) = 17
	col(4,12) = 17
	col(4,13) = 17
	col(4,14) = 17
	col(4,15) = 21
	col(4,16) = 21
	col(4,17) = 85

	return

	end



	subroutine lgnd( n, z, dz )
	common /legend/ nlevel, zlevel, dzlevl
	n = nlevel
	z = zlevel
	dz = dzlevl
	return
	end
