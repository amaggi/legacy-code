c
c*******************************************************************************
c
c    Subroutine chrsiz
c
c    Author - Danny Harvey
c
c*******************************************************************************
c
      subroutine chrsiz(height, ratio, slant)
c
c    Subroutine CHRSIZ will set the character size and slant of text
c    labelling.
c
c    inputs  - height = character height always in inches regardless of
c                       scaling.
c              ratio  = ratio of character cell width to height relative
c                       to "normal" (i.e. ratio=1. gives a "normal"
c                       looking character)
c              slant  = slant angle in degrees of character.  A positive
c                       slant will cause the top of the character to be
c                       displaced forward relative to the bottom of the
c                       character.  Slant=0. gives no slant.
c
      common /spc/ xll,yll,xur,yur,ASPECT,xc,xcl,fl
c
      common /partxt/ sinang,cosang,tansln,hite,rat
c
      common /npchr2/ scale
c
      data  rad     /  0.0174532925  /
      data  scale   /  1.0 /
c
      a = slant*rad
      tansln = tan(a)
      hite = height*xc*scale
      rat = ratio
c
      return
c
c    end of Subroutine CHRSIZ
c
      end
      subroutine gethit (hit)
      common /partxt/ sinang,cosang,tansln,hite,rat
      hit = hite
      return
      end
