c Code by Sylvana Pilidou, March 5, 2002  to pick right response for an event from a response file formatted like

c PABB LHZ 1992,294-> 1998,014
c -2
c ZEROS  3
c 0.  0.
c 0.  0.
c 0.  0.
c POLES  4
c -0.01234  0.01234
c -0.01234 -0.01234
c -39.18  49.12
c -39.18 -49.12
c CONSTANT 0.388564E+7 muM
c PABB LHZ 1998,014-> 1999,100
c -2
c ZEROS  3
c 0.  0.
c 0.  0.
c 0.  0.
c POLES  4
c -0.01234  0.01234
c -0.01234 -0.01234
c -39.18  49.12
c -39.18 -49.12
c CONSTANT 0.388564E+15 muM

c and create another response file containing only one and correct response to be used by synt.lr.f code which calculates the synthetic seismograms for the CL inversions.

c INPUT FILES:
c   1. data file name e.g.   /home/pilidou/1999.182.02.01.35.0000.IU.PABB....LHZ.dat
c   2. original response (like above) e.g. "/home/pilidou/IU.PABB....lhz.orig"

c OUTPUT files:
c  response correct for event time e.g. ./IU.PABB....lhz :

c WORKS on its own, but will be embeded in synt.lr.f for the autowfm analysis.

c -----------------------------------------------------------------------------

       program choose_resp
c       implicit none          must comment it for code to compile on saturn without complaining about the lnblnk function , everything else is OK
       integer eventyear,eventday
       character dataname*150, yearchar*4, daychar*3
       character respname*150,outresp*150, sta*4, comp*3
       integer startyear, startday, endyear, endday
       character key*10, units*4
       real foo, poler, poleim, zeror, zeroim, constant
       integer i, inst, number, length, charstart
       integer charstart2, length2

c      read data file name
       write(*,*) 'data file name?'
       read(*,'(a)') dataname
c       dataname='2000.100'
       write(*,*) ' ************************************************'
       write(*,*) 'given data file name is ',dataname
       write(*,*) ' ************************************************'

       
c      get year and day from of file name
c      dataname will be of the form
c   /home/pilidou/jdhflksjhfljdh/lshflhfl/1999.182.....
c    data file name has fixed length of 41 characters.
c      get number of characters
       length=lnblnk(dataname)
c      get number of first character to start from (1999....)       
       charstart=length-40
c       write(*,*) length, charstart

c      get year, julian day
       yearchar=dataname(charstart:charstart+3)
c       yearchar=dataname(1:4)
       daychar=dataname(charstart+5:charstart+7)
c       daychar=dataname(6:8)
       write(*,*) yearchar, daychar

       read(yearchar,'(i4)') eventyear
       read(daychar,'(i3)') eventday

       write(*,*) 'event year is ', eventyear
       write(*,*) 'event day is ', eventday

c      read response file name (e.g. /home/pilidou/.../IU.PABB....lhz.orig
       write(*,*) 'resp PZ file name?'
       read(*,'(a)') respname
       write(*,*) ' ************************************************'
       write(*,*) 'given response file name is ',respname
c      respname='IU.PABB....lhz.orig'

c       do the same as for event file
       length2=lnblnk(respname)
c      output resp file  (= IU.PABB....lhz)     
       charstart2=length2-13-5
       outresp=respname(charstart2:charstart2+13)

c      read(respname,'(A14)') outresp

       write(*,*) 'Output respone file will be ', outresp
       write(*,*) ' ************************************************'

c      open input resp file 
       open(10,file=respname,status='unknown')
       
c      read header line with dates as characters      
111    read(10,*,end=17) sta,comp, startyear, startday, endyear, endday
       write(*,*) sta,comp,startyear, startday, endyear, endday


c      compare event times and resp times
c      first case: year lies within response range
       if(eventyear.gt.startyear.and.eventyear.lt.endyear)then
       write(*,*) '0. found response for event time'
       go to 30
c      second case: year same with start year of resp
       elseif(eventyear.eq.startyear.and.eventyear.lt.endyear)then
           if(eventday.ge.startday)then
	   write(*,*) 'found response for event time'
	   go to 30
	   else 
	   write(*,*) "haven't found response for event time"
	   go to 40
	   endif
c      third case: year same with end year of resp
       elseif(eventyear.gt.startyear.and.eventyear.eq.endyear)then	   
            if(eventday.le.endday)then
	    write(*,*) 'found response for event time'
	    go to 30
	    else
	    write(*,*) "haven't found response for event time"
	    go to 40
	    endif
c      fourth case: start and end time of resp same
       elseif(eventyear.eq.startyear.and.eventyear.eq.endyear)then	    
            if(eventday.ge.startday.and.eventday.le.endday)then
	    write(*,*) 'found response for event time'
	    go to 30
	    else
	    write(*,*) "haven't found response for event time"
	    go to 40
	    endif
       else
       write(*,*) "haven't found response for event time"
       go to 40
       endif

40     write(*,*) 'will read next response period...' 
c      read zeros which won't use....
       read(10,*) inst
       read(10,*) key,number
       do 12 i=1,number
       read(10,*) foo
   12  continue
c      read poles which won't use....
       read(10,*) key,number
       do 13 i=1,number
       read(10,*) foo
   13  continue
c      read constant which won't use....
       read(10,*) key
c      go back to read next response section       
       go to 111
       

17     write(*,*) "Reached end of resp file and"
       write(*,*) "haven't found response for event's time"
       write(*,*) 'WILL EXIT'
       go to 121
       

30     continue
       write(*,*)
       write(*,*) 'Will now read poles and zeros...'
       write(*,*)
       open(20,file=outresp,status='unknown')
       rewind 20
       write(20,*)  sta, startyear, startday, endyear, endday
       read(10,*) inst
       write(20,*) inst
       read(10,*) key,number
       write(*,*) 'key= ',key,'number= ',number
       write(20,'(''ZEROS '',i2)') number 
       do 23 i=1,number
       read(10,*) zeror,zeroim
       write(20,*) zeror,zeroim
   23  continue
       read(10,*) key,number
       write(*,*) 'key= ',key,'number= ',number
       write(20,'(''POLES '',i2)') number 
       do 33 i=1,number
       read(10,*) poler, poleim
       write(20,*) poler, poleim
   33  continue
       read(10,*) key, constant, units
       write(20,'(''CONSTANT '',e12.6,'' muM'')') constant
       
       close(10)
       close(20)
      
121    end

