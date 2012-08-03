c----------------------------------------------------------------------
c Dumb program to reformat vsfin and anvsfin files from 20 columns to
c 10 columns.
c----------------------------------------------------------------------

      program reformat
      real*8 c(360,720), a1(360,720), a2(360,720)
      real tt
      character*80 vs1, vs2, anvs1, anvs2
      integer ili, icol, idep

      write(*,*) 'This program reformats vsfin and anvsfin files from '
      write(*,*) '20 columns to 10 columns.'
      write(*,*) 'Name of input vsfin file :'
      read (*,'(a)') vs1
      write(*,*) 'Name of input anvsfin file :'
      read (*,'(a)') anvs1
      write(*,*) 'Name of output vsfin file :'
      read (*,'(a)') vs2
      write(*,*) 'Name of output anvsfin file :'
      read (*,'(a)') anvs2
      write(*,*)'Number of lines in input files (latitude points): '
      read (*,*) ili
      write(*,*)'Number of columns in input files (longitude points): '
      read (*,*) icol
      write(*,*)'Number of depths in input files: '
      read (*,*) idep

c----------------------------------------------------------------------
c     start process of reading old files and writing newfiles
c----------------------------------------------------------------------
      open(unit=11,file=vs1)
      open(unit=12,file=anvs1)
      open(unit=13,file=vs2)
      open(unit=14,file=anvs2)

      do 10, id=1,idep

        read(11,*) tt
        read(12,*) tt
        write(13,*) tt
        write(14,*) tt

        do 30 i=1,ili
30        read(11,'(20(f10.6,1x))')(c(i,j),j=1,icol)
        do 50 i=1,ili
50        read(12,'(20(f10.6,1x))')(a1(i,j),j=1,icol)
        do 60 i=1,ili
60        read(12,'(20(f10.6,1x))')(a2(i,j),j=1,icol)


        do 300 i=1,ili
300       write(13,'(10(f10.6,1x))')(c(i,j),j=1,icol)
        do 500 i=1,ili
500       write(14,'(10(f10.6,1x))')(a1(i,j),j=1,icol)
        do 600 i=1,ili
600       write(14,'(10(f10.6,1x))')(a2(i,j),j=1,icol)

10    continue


      close(11)
      close(12)
      close(13)
      close(14)

c----------------------------------------------------------------------
c     end of process
c----------------------------------------------------------------------


      end
