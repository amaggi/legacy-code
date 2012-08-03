c
c*******************************************************************************
c
c    Subroutine sclsiz
c
c*******************************************************************************
c
      subroutine sclsiz (scale)
c
      real*4             scale
c
c    sclsiz will scale character sizes either up, if scale is greater than 1.0,
c    or down if scale is less than 1.0. scale = 1.0 corresponds to normal size.
c
      common /npchr2/ sc
c
      common /partxt/ sinang,cosang,tansln,hite,rat
c
      hite = hite/sc
      sc = scale
      hite = hite*sc
c
      return
      end
