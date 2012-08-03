c
c*******************************************************************************
c
c    Subroutine axis
c
c    Author - Danny Harvey
c
c*******************************************************************************
c
      subroutine getsiz(c, ihite, iwidth)
c
c    Subroutine GETSIZ will return the height and width
c    for a particular character of a particular font.  Subroutine CFONT
c    must be called prior to a call to this subroutine in order to specify
c    the font.  An error exit will occur if this is not the case.
c
c    input   - c      = character (CHARACTER*1)
c
c    outputs - ihite  = character height
c              iwidth = character width
c
c    parameter NPMAX defined in subroutine CFONT
c
      parameter  (NPMAX = 10000)
c
      character*1 c
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
      ihite = ipoint(ifoff)
      iwidth = ipoint(ifoff+ic)
c
      return
c
c    error exits
c
c    character outside range
c
  900 write (6,901) c,ic
  901 format (' GETSIZ: character ',a1,' (integer ',i3,') ',
     1  'outside printable range - run stopped')
      stop
c
c    font buffer empty
c
  910 write (6,911)
  911 format (' GETSIZ: font buffer empty (no prior call to CFONT)',
     1  ' - run stopped')
      stop
c
c    end of Subroutine GETSIZ
c
      end
