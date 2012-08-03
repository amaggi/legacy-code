FC=g77
#FC=gfortran
FFLAGS = -Wall

#LIBDIR=/software/apps/seismo/lib
LIBDIR=${CODE}/lib
INCDIR=${CODE}/include
BINDIR=${CODE}/bin
LIB=libeost_utils.a
LIBSRC =

PROGS=diagnrj

LIBOBJ = ${LIBSRC:%.f=$(LIB)(%.o)}

INSTALL = install
RANLIB = ranlib

default: $(LIB) $(PROGS)

SUFFIXES = .f .o
.SUFFIXES: $(SUFFIXES)

.f.o:
	$(FC) ${CFLAGS} -c -o $@ $< 

${LIB}:${LIBOBJ} 

diagnrj:diagnrj.f
	$(FC) -o $@ $<  -L${LIBDIR} -L/usr/X11R6/lib -L/sw/lib -lsac -lX11 -lps

install:
	${INSTALL} -c -m 0644 $(LIB) $(LIBDIR)
	${INSTALL} -c -m 0644 $(INC) $(INCDIR)
	${RANLIB} ${LIBDIR}/$(LIB)
	for i in ${PROGS}; do ${INSTALL} -c -m 0755 $$i $(BINDIR); done

clean:
	rm $(LIB)
