       program call_distaz
       real evla, evlo, stla, stlo, azm, bzm, ddg, dkm
       write(*,*) 'Input coordinates (lat lon) of event:'
       read(*,*) evla, evlo
       write(*,*) 'Input coordinates (lat lon) of station:'
       read(*,*) stla, stlo
       call distaz(evla,evlo,stla,stlo,azm,bzm,ddg,dkm)
       write(*,'("DDG = ",f6.2)') ddg
       write(*,'("DKM = ",f8.2)') dkm
       write(*,'("AZM = ",f6.2)') azm
       write(*,'("BZM = ",f6.2)') bzm
       end program call_distaz