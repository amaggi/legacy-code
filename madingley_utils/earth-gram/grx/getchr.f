c
c*******************************************************************************
c
c    Subroutine getchr
c
c    Author - Danny Harvey
c
c*******************************************************************************
c
      subroutine getchr(c, nstrks, npoint, ioff, ix, iy)
c
c    Subroutine GETCHR will return the plot points
c    for a particular character of a particular font.  Subroutine CFONT
c    must be called prior to a call to this subroutine in order to specify
c    the font.  An error exit will occur if this is not the case.
c
c    input   - c      = character (CHARACTER*1)
c
c    outputs - nstrks = number of strokes in the character
c              npoint(nstrks) = number of plot points per stroke
c              ioff(nstrks) = beginning point in IX and IY arrays per stroke
c              ix,iy = plot point arrays
c
c    parameter NPMAX defined in subroutine CFONT
c
      parameter  (NPMAX = 10000)
c
      character*1 c
      dimension npoint(27),ioff(27),ix(150),iy(150)
c
      common /fonts/ ifoff,ipoint(NPMAX)
      integer*2 ipoint
      integer*2 ifoff
c
c    check to see if character in range
c
      ic = ichar(c)
      ic = ic - 31
      if (ic .lt. 1 .or. ic .gt. 95)  go to 900
c
c    check to see if font buffer empty
c
      if (ifoff .eq. -1)  go to 910
c
c    fill in output variables and return
c
      joff = ipoint(ifoff+95+ic)
      koff = ifoff + 190 + joff
      nstrks = ipoint(koff)
      if (nstrks .eq. 0)  return
      j = 0
      do 100  i = 1,nstrks
      koff = koff + 1
      np = ipoint(koff)
      np = np/2
      npoint(i) = np
      ioff(i) = j + 1
      do 110  k = 1,np
      j = j + 1
      koff = koff + 2
      ix(j) = ipoint(koff-1)
  110 iy(j) = ipoint(koff)
  100 continue
c
      return
c
c    error exits
c
c    character outside range
c
  900 write (6,901) c,ic
  901 format (' GETCHR: character ',a1,' (integer ',i3,') ',
     1  'outside printable range - run stopped')
      stop
c
c    font buffer empty
c
  910 write (6,911)
  911 format (' GETCHR: font buffer empty (no prior call to CFONT)',
     1  ' - run stopped')
      stop
c
c    end of Subroutine GETCHR
c
      end
