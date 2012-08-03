c
c $Id: datetime.f,v 1.1.1.1 2002/07/12 11:15:19 maggi Exp $
c $Log: datetime.f,v $
c Revision 1.1.1.1  2002/07/12 11:15:19  maggi
c
c
c Revision 1.1  2002/05/23 10:28:28  maggi
c Initial revision
c
c
      subroutine datetime(iy,id,ih,im,sec,tshift)
C     given a date and a time shift in second, caluculate the new date.

      sec=sec+tshift
      ism=int(sec/60.)
      sec=mod(sec,60.)
       if(sec.lt.0.) then
       sec=sec+60
       ism=ism-1
       endif
      im=im+ism
      imh=int(im/60)
      im=mod(im,60)
       if(im.lt.0) then
       im=im+60
       imh=imh-1
       endif
      ih=ih+imh
      ihd=int(ih/24)
      ih=mod(ih,24)
       if(ih.lt.0) then
       ih=ih+24
       ihd=ihd-1
       endif
      id=id+ihd-1
      idy=0
 10   if(id.ge.nday(iy+idy)) then
       id=id-nday(iy+idy)
       idy=idy+1
       goto 10
      endif
 20   if(id.lt.0) then
       idy=idy-1
       id=id+nday(iy+idy)
       goto 20
      endif
      id=id+1
      iy=iy+idy

      return
      end


      function nday(iyear)
C     How many days in a year?
      nday=365
      if(mod(iyear,4).eq.0) nday=366
      return
      end
