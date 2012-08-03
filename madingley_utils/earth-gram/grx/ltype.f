c
c*******************************************************************************
c
c    Subroutine ltype
c
c    Author - Danny Harvey
c
c*******************************************************************************
c
      subroutine ltype(n)
c
c    routine ltype sets the line type for graphing.
c
c     n  =  0 - solid
c        =  1 - dotted
c        =  2 - dot-dashed
c        =  3 - short-dashed
c        =  4 - long-dashed
c
      common  /ocflag/  imflag,iocflg,iltp
C
      iltp = n
      if (n .lt. 0 .or. n .gt. 10)  iltp = 0
      return
      end
