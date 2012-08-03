        Subroutine calendar (year, month, iday, jday)
c
c	Subroutine to convert year, month, and day into julian day, taking
c	account of leap year.
c
        integer year, iday

        leap=1
        itest=mod(year,4)
        if (itest .ne. 0) leap=0
c
	if (jday .le. 31) then
	   iday = jday
	   month = 1
	elseif (jday .le. (59 + leap)) then
	   iday = jday - 31
	   month = 2
	elseif (jday .le. (90 + leap)) then
	   iday = jday - (59 + leap)
	   month = 3
	elseif (jday .le. (120 + leap)) then
	   iday = jday - (90 + leap)
	   month = 4
	elseif (jday .le. (151 + leap)) then
	   iday = jday - (120 + leap)
	   month = 5
	elseif (jday .le. (181 + leap)) then
	   iday = jday - (151 + leap)
	   month = 6
	elseif (jday .le. (212 + leap)) then
	   iday = jday - (181 + leap)
	   month = 7
	elseif (jday .le. (243 + leap)) then
	   iday = jday - (212 + leap)
	   month = 8
	elseif (jday .le. (273 + leap)) then
	   iday = jday - (243 + leap)
	   month = 9
	elseif (jday .le. (304 + leap)) then
	   iday = jday - (273 + leap)
	   month = 10
	elseif (jday .le. (334 + leap)) then
	   iday = jday - (304 + leap)
	   month = 11
	else
	   iday = jday - (334 + leap)
	   month = 12
	endif
	return
	end
