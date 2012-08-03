c     bins and averages velocities according to ages
c     both velocities and ages are masked with -999=NAN
c     velocity and age grids are IDENTICAL

      parameter (MAXNPTS=360*720, MAXBINS=1000)
      real lat(MAXNPTS), lon(MAXNPTS), age(MAXNPTS), vs(MAXNPTS)
      real binsize, binstep, binsum(MAXBINS)
      real agemin, agemax
      integer npts, nbins, nan, bincount(MAXBINS)

      character*80 agename, vsname, outname

      iage=10
      ivs=11
      iout=12
      nan=-999

c     get relevant informtion from user
      write(*,*) 'Number of points in .xyz files: '
      read (*,*) npts
      if (npts.gt.MAXNPTS) goto 1001

      write(*,*) 'Name of age .xyz file :'
      read (*,*) agename
      write(*,*) 'Name of masked vs .xyz file :'
      read (*,*) vsname

      write(*,*) 'Size of age bin (Ma) : '
      read (*,*) binsize
      write(*,*) 'Bin step (Ma, less than bin size) : '
      read (*,*) binstep

      write(*,*) 'Name of output file :'
      read (*,*) outname

c     read the two data files into memory
      open(iage,file=agename)
      open(ivs,file=vsname)

      agemin=1000
      agemax=0
      do 10 i=1,npts  
        read(iage, *) lat(i), lon(i), age(i)
        if (age(i).gt.nan.and.age(i).lt.agemin) agemin = age(i)
        if (age(i).gt.agemax) agemax = age(i)
        read(ivs, *) bidlat, bidlon, vs(i)
10    continue

      close(iage)
      close(ivs)

c     calculate number of bins
      nbins=int((agemax-agemin-binsize)/binstep)
      if (nbins.gt.MAXBINS) goto 1002

c     write(*,*) '# Min age : ',agemin
c     write(*,*) '# Max age : ',agemax
c     write(*,*) '# Num bins: ',nbins
c     write(*,*) '# Binsize : ',binsize
c     write(*,*) '# Binstep : ',binstep

c     initialise bins
      do 15, ibin=1,nbins
        binsum(ibin)=0
        bincount(ibin)=0
15    continue

c     write(*,*) 'Starting loop'

c     start loop over points
      do 20, i=1,npts
        if(vs(i).gt.nan.and.age(i).gt.nan) then
c         start loop over bins
          do 21, ibin=1,nbins
            bmin=agemin+(ibin-1)*binstep
            bmax=agemin+binsize+(ibin-1)*binstep
            if(age(i).ge.bmin.and.age(i).lt.bmax) then
              binsum(ibin) = binsum(ibin) + vs(i)
              bincount(ibin) = bincount(ibin) + 1
            endif
21        continue
        endif
20    continue

c     write(*,*) 'Finished loop'

c     output
      open(iout,file=outname)
      write(iout,*) '# Min age : ',agemin
      write(iout,*) '# Max age : ',agemax
      write(iout,*) '# Num bins: ',nbins
      write(iout,*) '# Binsize : ',binsize
      write(iout,*) '# Binstep : ',binstep
      write(iout,*) '# Bincenter, binmin, binmax, average, nvalues'
      write(iout,*) '---------------------------------------------'
      do 30, ibin=1,nbins
        bmin=agemin+(ibin-1)*binstep
        bmax=agemin+binsize+(ibin-1)*binstep
        write(iout,*) (bmax+bmin)/2, bmin, bmax, 
     *                binsum(ibin)/bincount(ibin), bincount(ibin)
30    continue
      close(iout)


      return

c     ERROR MESSAGES
1001  write(*,*) 'npts ',npts, ' > MAXNPTS ', MAXNPTS
      stop 'Increase MAXNPTS in bin_ages.f'
1002  write(*,*) 'nbins',nbins,' > MAXBINS ', MAXBINS
      stop 'Increase MAXBINS in bin_ages.f'
      end
