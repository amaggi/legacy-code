#FC=g77
FC=gfortran
FFLAGS = -Wall -O5

LIBDIR=${CODE}/lib
INCDIR=${CODE}/include
BINDIR=${CODE}/bin
LIB=
LIBSRC =

PROGS=minos

LIBOBJ = ${LIBSRC:%.f=$(LIB)(%.o)}

INSTALL = install
RANLIB = ranlib

default: $(LIB) $(PROGS)

SUFFIXES = .f .o
.SUFFIXES: $(SUFFIXES)

.f.o:
	$(FC) ${CFLAGS} -c -o $@ $< 

${LIB}:${LIBOBJ} 

minos:minos.f
	$(FC) -o $@ $< -L${LIBDIR} -lm

install:
	for i in ${PROGS}; do ${INSTALL} -c -m 0755 $$i $(BINDIR); done
#	${INSTALL} -c -m 0644 $(LIB) $(LIBDIR)
#	${INSTALL} -c -m 0644 $(INC) $(INCDIR)
#	${RANLIB} ${LIBDIR}/$(LIB)

clean:
	rm $(LIB)
