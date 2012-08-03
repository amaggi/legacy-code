c     SAC2MT5 input.txt File Creater  (Onur TAN)

      integer yy,mo,day,ho,mi,idep
      real    lat,lon,dep,sec,mw
      character*20 file,null

      print *, "SAC2MT5 input.txt File Creater (by Onur TAN)"
      
  10  write(*,11)
  11  format("MT5 travel time file :  ",$)
      read(*,'(A20)') file
      open(1,file=file,status="old",err=99)
      goto 15

  99  print *, "Travel time file error. Sellect from list:"
      call system ('ls')
      print *, " "
      goto 10


  15  write(*,16)
  16  format("Mw :  ",$)
      read(*,*) mw

      read(1,*)null
      read(1,20)yy,mo,day,lat,lon,dep,ho,min,sec
   20 format(i2,i2,i2,6X,f7.3,3X,f7.3,3X,f7.3,3X,i2,1X,i2,1X,f4.1)

      idep = dep    ! integer depth 

      open(2,file="input.txt")

      write(2,*) " "

      write(2,30)yy,mo,day,ho,min,sec,lat,lon,idep,mw,
     & yy,mo,day,ho,min,sec,mw
   30 format(i4.4," (000) ",4(i2.2,1X),f4.1,1X,2f7.3,1X,i3,1X,
     & f3.1,1X,5i2.2,f4.1," HEADER Mw=",f3.1 )

      write(2,*) " "



      stop
      end
