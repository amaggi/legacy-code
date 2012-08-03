       program call_jday
       integer year, month, day, julday
       write(*,*) 'Input year, month, day'
       read(*,*) year, month, day
       call julian(year,month,day,julday)
       write(*,'("JULIAN DAY = ",i5)') julday
       end program call_jday
