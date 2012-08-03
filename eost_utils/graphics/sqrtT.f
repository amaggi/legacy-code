       real bstart, bend, bstep
       real factor, offset
       real depth, age
       character filename*80

       write(*,*) 'Start, end and step values of T : '
       read(*,*) bstart, bend, bstep
 
       write(*,*) 'Proportionality factor and vertical offset : '
       read(*,*) factor, offset

       write(*,*) 'Output filename'
       read(*,*) filename

       nsteps=(bend-bstart)/bstep + 1
       
       open(12, file=filename)
       do 10, i=1,nsteps
         age=bstart + (i-1)*bstep
         depth=offset+factor*sqrt(age)
         write(12,*) age, depth
10     continue
       close(12)

       end
