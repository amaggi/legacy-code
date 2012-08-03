      ! Test how much memory can be allocated
      program alloctest2
      implicit none
      integer :: i
      real(8), allocatable :: a(:,:)
      do i = 10000, huge(i), 1000
      print *, 'Trying to allocate: ', i**2*8._8 / 1024._8/1024_8,' MB.'
      allocate (a(i,i))
      deallocate (a)
      end do
      end program alloctest2
