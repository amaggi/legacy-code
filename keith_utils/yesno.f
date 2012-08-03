	logical function yesno(quest)
c
c	Interactive i-o for logical variables -- yesno must be declared 
c	logical in calling program
c
	character quest*(*), answer*1
	logical lanswr
	write(*,'(80(a1,$))') (quest(j:j),j=1,len(quest))
	read(*,'(a1)') answer
	lanswr = .false.
	if(answer .eq. 'y') lanswr = .true.
	yesno = lanswr
	return
	end
