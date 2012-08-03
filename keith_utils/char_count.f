	subroutine char_count(root,ncharacter)
	character root*50
	ncharacter = 0
10	ncharacter = ncharacter + 1
	if(root(ncharacter:ncharacter) .eq. ' ') then
	   ncharacter = ncharacter - 1
	   return
	else
	   go to 10
	endif
	return
	end

