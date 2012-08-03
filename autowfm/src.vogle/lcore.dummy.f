*************************************************************
*
*   Subroutine bidons pour ne pas dessiner
*
*
       subroutine SetLinewidth(w)
       return
       end
       subroutine SetLinestyle(i)
       return
       end
       subroutine DelAllRetainSegs()
       return
       end
       subroutine DelRetainSegment()
       return
       end
       subroutine DelRetSegment()
       return
       end
       subroutine CloseRetainseg()
       return
       end
       subroutine NewFrame
       return
       end
       subroutine CreateRetainSeg(i)
       return
       end
       subroutine CreateRetSeg(i)
       return
       end
       subroutine MoveAbs2(x,y)
       return
       end
       subroutine LineAbs2(x,y)
       return
       end
       subroutine MoveRel2(x,y)
       return
       end
       subroutine LineRel2(x,y)
       return
       end
       subroutine PolyLineAbs2(x,y,n)
       return
       end
       subroutine Text
       return
       end
       subroutine SetFont
       return
       end
       subroutine SetCharPrecision
       return
       end
       subroutine SetCharPrecis
       return
       end
       subroutine SetCharSize(sx,sy)
	common/jjl/xmin,xmax,ymin,ymax,sizex,sizey
	sizex=sx
	sizey=sy
       return
       end

       subroutine inqwindow(xmi,xma,ymi,yma)
	common/jjl/xmin,xmax,ymin,ymax,sizex,sizey
	xmi=xmin
	xma=xmax
	ymi=ymin
	yma=ymax
       return
       end
       subroutine deselectvwsurf
       return
       end
       subroutine selectvwsurf
       return
       end
       subroutine bw2dd
       return
       end
       subroutine initializevwsurf
       return
       end
       subroutine setechosurface
       return
       end
       subroutine initializedevice
       return
       end
       subroutine pixwindd
       return
       end
       subroutine setecho
       return
       end
       subroutine setviewport2
       return
       end
       subroutine terminatecore
       return
       end
       subroutine setwindow(xmi,xma,ymi,yma)
	common/jjl/xmin,xmax,ymin,ymax,sizex,sizey
	xmin=xmi
	xmax=xma
	ymin=ymi
	ymax=yma
       return
       end
       subroutine closetempseg
       return
       end
       subroutine createtempseg
       return
       end
       subroutine initializecore
       return
       end
       subroutine inqtextextent2(txt,x,y)
	common/jjl/xmin,xmax,ymin,ymax,sizex,sizey
	character *(*) txt
	x=lnblnk(txt)*sizex
	y=sizey
       return
       end
       subroutine exit
       return
       end
       subroutine inline_st_int
       return
       end
       subroutine inline_st_float
       return
       end
       subroutine inline_ld_int
       return
       end
       subroutine inline_ld_float
       return
       end
       subroutine initdes
       return
       end
       subroutine findes
       return
       end
