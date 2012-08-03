      subroutine remove(fname)
c----------------------------------------------------------------------
c
c   Deletes file "fname" (if it exists).
c   Assumes "fname" is dimensioned 50 characters long.
c
c   Non-system routines: none
c
c
c   The first call to 'access()' returns junk.
c   I don't know why.
c
c
c
c----------------------------------------------------------------------
      logical nothere
      character fname*80
      integer acc
c
c   Do nothing if there is no such file
c
      acc = 0
      inquire(file=fname,exist=nothere)
c     acc = access(fname,'w')
      if (acc .ne. 0) then
	write(*,*) 'You do not have permission to delete ',fname
        return
      else
	call unlink(fname)
        return
      endif
c
c
      end
