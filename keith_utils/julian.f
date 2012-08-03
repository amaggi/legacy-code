        Subroutine julian(year,month,day,julday)
c
c	Subroutine to convert year, month, and day into julian day, taking
c	account of leap year.
c
        integer year, day

        leap=1
        itest=mod(year,4)
        if (itest .ne. 0) leap=0

        goto(5,10,15,20,25,30,35,40,45,50,55,60) month
c
5       julday=day
        goto 100
10      julday=31+day
        goto 100
15      julday=59+day+leap
        goto 100
20      julday=90+day+leap
        goto 100
25      julday=120+day+leap
        goto 100
30      julday=151+day+leap
        goto 100
35      julday=181+day+leap
        goto 100
40      julday=212+day+leap
        goto 100
45      julday=243+day+leap
        goto 100
50      julday=273+day+leap
        goto 100
55      julday=304+day+leap
        goto 100
60      julday=334+day+leap
100     return
        end
