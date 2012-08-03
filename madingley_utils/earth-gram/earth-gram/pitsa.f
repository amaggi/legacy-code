c
c $Id: pitsa.f,v 1.1.1.1 2002/07/12 11:15:20 maggi Exp $
c $Log: pitsa.f,v $
c Revision 1.1.1.1  2002/07/12 11:15:20  maggi
c
c
c Revision 1.1  2002/05/23 10:28:34  maggi
c Initial revision
c
c
          subroutine prnt_hdr
      include '../include/sizes.inc'
      include '../include/commons.inc'

          common/pitsa/selev,xphy,xdist

          write(*,'(/,a)') ' Trace Information'
          write(*,'(3(a,a4))')
     &          '  Station=',sta,' Channel=',chn,' Type=',typ
          write(*,'(3(a,f8.3))')
     &        '  Latitude=',t1,' Longitude=',p1,' Elevation=',selev
          write(*,'(a,f6.3,a,i8)')
     &        '  Sample interval=',dt,' Npts=',nscan
          write(*,'(a,i4,a1,i3,2(a1,i2),a1,f6.3)')
     &           '  Start time=',iy,':',id,':',ih,':',im,':',ss
          write(*,'(/,a)') ' Event Information'
          write(*,'(a,i4,a1,i3,2(a1,i2),a1,f6.3)')
     &           '  Origin time=',jy,':',jd,':',jh,':',jm,':',sec
          write(*,'(a,f9.3,a,f8.3)') 
     &           '  Source-reciver distance=',xdist,' Azimuth=',xphy
          write(*,'(a,f8.3,/)') '  Source depth=',d0
          return
          end

      subroutine juldat(year,month,day,julday)
c******************************************************************
c   given the usual year, month, and day (as integer numbers) subroutine
c   juldat returns integer julian day julday.  it properly accounts for
c   the leap years through 2099.  if the month or day are illegal dates,
c   then julday is set to zero.  

c subroutines called
c       no subroutines called

c functions used             
c       no user defined functions used

c entries into subroutine
c       gredat     

c----- written originally by ray buland
c----- last date modified: prior to jan. 15, 1992 
c******************************************************************

c variables
      integer*4 year, ! years ad 
     &       month, ! month of year
     &         day, ! day of month 
     &      julday, ! julian day of year
     &     dpm(12), ! array of days per month
     &         feb  ! number of days in feb.

c data statements
      data dpm/31,28,31,30,31,30,31,31,30,31,30,31/,feb/28/
                                                      

c---- calculate the proper number of days for february
      if(mod(year,4).eq.0) then
         dpm(2) = feb + 1
      else
         dpm(2) = feb
      endif
                
c---- if a valid month or day then continus else set julian day to 0
      if(month.le.0.or.month.gt.12.or.day.le.0.or.day.gt.dpm(month))then
           julday = 0
      else 
           ! add up the number of days in the months up to input month
           julday = day
           if(month .gt .1) then
              im = month - 1
              do i = 1,im
                julday = julday + dpm(i)
              enddo
           endif
      endif
      return 
	 
      entry gredat(year,month,day,julday)
c******************************************************************
c   given the integer year and julian day, entry gredat returns the
c   correct integer month and day of the month.  gredat correctly treats
c   leap years through 2099.  if the julian day is illegal, then a zero
c   month and day are returned.

c----- written originally by ray buland
c----- last date modified: prior to jan. 15, 1992 
c******************************************************************
                                                        
c---- test for valid julian day. if invalid set month/day to 0
c     otherwise if valid calculate month and day
      if(julday.le.0.or.julday.gt.366-min0(mod(year,4),1)) then
             day = 0
           month = 0             

      else    
c.......  adjust number of days in february if leap year
          if(mod(year,4).eq.0) then
            dpm(2) = feb + 1
          else 
            dpm(2) = feb
          endif

c.......  initialize output variables
          day   = julday
          month = 1

c.......  loop though the months subtracting days and 
c         increasing the months
          do i=1,12
             if(day. gt .dpm(i)) then
                  day   = day - dpm(i)
                  month = month + 1
             endif
          enddo
      endif 

      return
      end
