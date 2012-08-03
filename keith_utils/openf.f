      subroutine openf(lun,fname,stat)
c-----------------------------------------------------------------------
c
c   Opens SEQUENTIAL ACCESS, FORMATTED files with decent error handling.
c
c   Input via call:
c
c   lun    -=- unit number to open
c   fname  -=- file specification (assumed 50 char length)
c   stat  -=- desired open status ('new' or 'old)
c
c   Non-system routines: remove
c
c-----------------------------------------------------------------------
      logical there
      character blank*80, fname*80, state*3, ans*1, stat*3
      blank=' '
      state=stat
      if (fname .eq. 'prn') then
       fname='LPT1:'
       state='old'
      elseif (fname .eq. 'lpt1:') then
       fname='LPT1:'
       state='old'
      elseif (fname .eq. 'lpt2:') then
       fname='LPT2:'
       state='old'
      elseif (fname .eq. 'con:') then
       state='old'
      endif
c
c   Determine if fname exists or not
c
10    INQUIRE(file=fname,exist=there)
c
c   OLD files
c
      if (state .eq. 'old' .and. .not. there) then
	write(*,15) state, fname
15	format(/,' Open error for "',a3,'" file: ',a50)
	write(*,'(a)')   ' File does not exist.  Try again: '
	read(*,'(a)') fname
	  if (fname .ne. blank) then
	   goto 10
	  else
	   stop 'OPENF stop'
	  endif
      elseif (state .eq. 'old' .and. there) then
	open(unit=lun,file=fname,status='old',iostat=icode)
	if (icode .eq. 0) return
	write(*,20) state, fname, icode
20	format(' Open error for "',a3,'" file: ',a50,/,'Error code =',i5)
c
c   NEW files
c
      elseif(state .eq. 'new' .and. there) then
      if (fname .eq. 'con'   .or. fname .eq. 'LPT1:' .or.
     &	  fname .eq. 'LPT2:' .or. fname .eq. 'LPT3:') then
	open(unit=lun,file=fname)
	return
      else
 	write(*,25) fname
25	format(' Already have a file: ',a)
26	write(*,'(a)') ' Delete (d), append (a), or new file (n)? '
	read(*,'(a)') ans
	if (ans .eq. blank) stop 'OPENF stop'
	if (ans .ne. 'd' .and. ans .ne. 'a' .and. ans .ne. 'n') goto 26
	if (ans .eq. 'd') then
c         itmp = access(fname,'w')
	  call remove(fname)
	  goto 10
	elseif (ans .eq. 'a') then
	  open(unit=lun,file=fname,status='old',iostat=icode)
	  if (icode .eq. 0) then
27	   read(lun,'(a)',end=28) ans
	   goto 27
28	   backspace lun
	   return
	  else
	   write(*,20) state, fname, icode
	  endif
	else
	  write(*,'(a)') ' New file name: '
	  read(*,'(a)') fname
	  if (fname .ne. blank) then
	   goto 10
	  else
	   stop 'OPENF stop'
	  endif
	endif
       endif
      elseif (state .eq. 'new' .and. .not. there) then
	open(unit=lun,file=fname,status='new',iostat=icode)
	if (icode .eq. 0) return
	write(*,20) state, fname, icode
c
c   Unrecognized status
c
      else
	write(*,30) state
30	format(' OPENF error. Unrecognized status: ',a3)
      endif
      stop 'OPENF stop'
      end
