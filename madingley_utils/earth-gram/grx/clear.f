c
c*******************************************************************************
c
c    Subroutine clear
c
c    Author - Danny Harvey
c
c*******************************************************************************
c
      subroutine clear
c
c     routine clear causes the screen to clear.
c
c
      common /spc/ xll,yll,xur,yur,ASPECT,xc,xcl,fl
c
      common /drwclr/ idraw
c
      data  idraw  / 0 /
c
c     if (idraw .eq. 1) then
	call hdclr
	ix = xll*fl + 0.5
	iy = yll*fl + 0.5
	iw = (xur-xll)*fl + 1.5
	ih = (yur-yll)*fl + 1.5
	call hdclrl(ix, iy, iw, ih)
	idraw = 0
c     end if
c
      return
      end
