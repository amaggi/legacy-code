	subroutine asktxt(quest,answer)
c
c	Interactive i-o for character strings -- string returned may 
c	be a maximun of 32 characters
c
	character answer*32,quest*(*)
	write(*,'(80(a1,$))') (quest(j:j),j=1,len(quest))
	read(*,'(a32)') answer
	return
	end
