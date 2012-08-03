LIBDIR=/software/apps/seismo/lib

LIB= libgrx.a
DATA = fonts.bin

default:${LIB}

LIBSOURCES = axis.f  box.f  cfont.f  chrdir.f  chrsiz.f  circle.f  clear.f  clrrgn.f  cpyrgn.f  cursor.f  drawv.f  drawx.f  finitt.f  fplot.f  getchr.f  getsiz.f  imgutl.f  initt.f  line.f  laxis.f  llaxis.f  ltype.f  movev.f  movex.f  nalpha.f  nclose.f  ncontour.f  ncrimg.f  nfill.f  ngrid.f  nhist.f  ninimg.f  nlstyl.f  nopen.f  npause.f  nplot.f  npnimg.f  ntype.f  plot1.f  sclsiz.f  setbac.f  setbg.f  selpen.f  setfg.f  setfor.f  setdim.f  setscl.f  stack.f  tdraw.f  tek.f  text.f  tics.f  vt102.f  xmap.f   contrd_.c  hdback_.c  hdbatoff_.c  hdbaton_.c  hdbg_.c  hdcircle_.c  hdclr_.c  hdclrg_.c  hdcprg_.c  hdcrimg_.c  hddash_.c  hddraw_.c  hdend_.c  hdfill_.c  hdfore_.c  hdfplt_.c  hdgptr_.c  hdinit_.c  hdkild_.c  hdmove_.c  hdmptr_.c  hdopen_.c  hdptimg_.c  hdsetclp_.c  hdstrk_.c  hdtext_.c  hdthik_.c  my_system.c  nxplotdlib.c  nxplotdlsv.c  parse_on_1char.c

INCLUDE = nxplotdlib_types.h

MAN3       = axis.3f  box.3f  cfont.3f  chrdir.3f  chrsiz.3f  circle.3f  clear.3f  clrrgn.3f  cpyrgn.3f  cursor.3f  finitt.3f  fplot.3f  getchr.3f  getsiz.3f  initt.3f  laxis.3f  line.3f  llaxis.3f  ltype.3f  ncontour.3f  niceplot.3f  nplot.3f  ntype.3f  plot1.3f  sclsiz.3f  setbac.3f  setbg.3f  setdim.3f  setfg.3f  setfor.3f  setscl.3f  text.3f  tics.3f  xmap.3f  ymap.3f 
EXCLUDE=test 
#include	$(JSPCMAKE) 
CFLAGS =	-O -I/usr/openwin/share/include
 FFLAGS =	-O
CC=cc
DIRS=

OBJS=	${LIBSOURCES:%.f=${LIB}(%.o)} ${LIBSOURCES:%.c=${LIB}(%.o)}
$(LIB) : ${OBJS}

clean:
	rm ${LIB}


install: ${LIB}
	install -c -m 0644 ${LIB} ${LIBDIR}
	ranlib ${LIBDIR}/${LIB}

#install_solaris: ${LIB}
#	install -c ${LIBDIRS} -m 0644 ${LIB} 
#	ranlib ${LIBDIRS}/${LIB}

