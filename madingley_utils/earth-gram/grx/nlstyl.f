c
c*******************************************************************************
c
c    Subroutine nlstyl
c
c    Author - Danny Harvey
c
c*******************************************************************************
c
      SUBROUTINE NLSTYL(I)
C
C    This routine will set the graphics line style
C
      dimension ipat(10)
c
      data iold / -1 /
c
      if (i .eq. iold) return
      iold = i
      if (i .eq. 0) then
	call hddash (0, ipat)
      else if (i .eq. 1) then
	ipat(1) = 5
	ipat(2) = 8
	call hddash (2, ipat)
      else if (i .eq. 2) then
	ipat(1) = 40
	ipat(2) = 40
	ipat(3) = 280
	ipat(4) = 40
	call hddash (4, ipat)
      else if (i .eq. 3) then
	ipat(1) = 150
	ipat(2) = 50
	call hddash (2, ipat)
      else if (i .eq. 4) then
	ipat(1) = 350
	ipat(2) = 50
	call hddash (2, ipat)
      else if (i .eq. 5) then
	call hddash (0, ipat)
      else if (i .eq. 6) then
	ipat(1) = 20
	ipat(2) = 40
	call hddash (2, ipat)
      else if (i .eq. 7) then
	ipat(1) = 40
	ipat(2) = 40
	ipat(3) = 280
	ipat(4) = 40
	call hddash (4, ipat)
      else if (i .eq. 8) then
	ipat(1) = 150
	ipat(2) = 50
	call hddash (2, ipat)
      else if (i .eq. 9) then
	ipat(1) = 350
	ipat(2) = 50
	call hddash (2, ipat)
      else if (i .eq. 10) then
	call hddash (0, ipat)
      end if
      return
C
      END
