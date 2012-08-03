c
c $Id: char_int.f,v 1.1.1.1 2002/07/12 11:15:18 maggi Exp $
c $Log: char_int.f,v $
c Revision 1.1.1.1  2002/07/12 11:15:18  maggi
c
c
c Revision 1.1  2002/05/23 10:28:27  maggi
c Initial revision
c
c       
	character function char_int(nm)
c      input an interger and return the last 4 digits as a character
c      character char_int*4

       n=nm
       do 10 i=4,1,-1
       m=int(n/10)
       id=n-10*m
       n=m
  10   char_int(i:i)=char(id+48)

       return
       end
