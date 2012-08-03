c 
c $Id: irisresp.f,v 1.1.1.1 2002/07/12 11:15:20 maggi Exp $
c $Log: irisresp.f,v $
c Revision 1.1.1.1  2002/07/12 11:15:20  maggi
c
c
c Revision 1.1  2002/05/23 10:28:32  maggi
c Initial revision
c
c
       complex function irisresp(w)
       double precision pi
       character stat*5,cha*5,units*4,file*31,opt*4,datime*31
       character char_int*4,date*4
       include '../include/sizes.inc'
       include '../include/commons.inc'

       stat=sta
       cha=chn

c  output displacement responses with unit counts/meter
       units='dis'
       file=''
       opt='-v'
       date=char_int(iy)
       datime(1:4)=date
       datime(5:5)=','
       date=char_int(id)
       datime(6:8)=date(2:4)
       datime(9:9)=','
       date=char_int(ih)
       datime(10:11)=date(3:4)
       datime(12:12)=':'
       date=char_int(im)
       datime(13:14)=date(3:4)
       datime(15:15)=':'
       date=char_int(int(ss))
       datime(16:18)=date(3:4)

       pi=dacos(-1.d0)
       fr=w/(2*pi)
       prl=1.
       pim=0.

       call respc(fr,stat,cha,units,file,opt,datime,prl,pim)
       irisresp=cmplx(prl,pim)


       return
       end
