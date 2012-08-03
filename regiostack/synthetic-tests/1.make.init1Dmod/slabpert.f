      subroutine slab_read(filename,ndepths,depths,npts,lats,lons)

      include 'slabs.h'

      integer ndepths,npts(NSLABDEPTHS)
      real lats(NSLABDEPTHS,NSLABPTS),lons(NSLABDEPTHS,NSLABPTS)
      real depths(NSLABDEPTHS)
      character*256 filename

      character*40 dummy
      
      iunit=17

      open(unit=iunit,file=filename,status='old')
      read(iunit,'(a)') dummy
      ndepths=1
      while (ndepths.lt.NSLABDEPTHS) do

c       read number of lat/lon pairs at this depth
        read(iunit,*,err=666) npts(ndepths)

c       sanity check on size of contour
        if (npts(ndepths).gt.NSLABPTS) then
          write(*,*) 'Too many lat/lon points per depth in ',filename
          stop
        endif

c       read lat/lon/depth
        do i=1,npts(ndepths)
          read(iunit,*) lons(ndepths,i),lats(ndepths,i),depths(ndepths)
        enddo
        
        ndepths=ndepths+1
      enddo

c     control reaches here if we have hit NSLABDEPTHS contours or
c     the file has ended
666   ndepths=ndepths-1
      write(*,*) 'Read ',ndepths,' contours from file ',filename

      close(iunit)
      
      end subroutine
