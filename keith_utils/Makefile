FC=gfortran
FFLAGS = -m32

#LIBDIR=/software/apps/seismo/lib
CODE=${HOME}/code
LIBDIR=${CODE}/lib
INCDIR=${CODE}/include
BINDIR=${CODE}/bin
LIB=libutil.a
INC=kutils.h
LIBSRC =asktxt.f caldat.f calendar.f distaz.f julday.f julian.f openf.f remove.f yesno.f

PROGS=distaz jday

LIBOBJ = ${LIBSRC:%.f=$(LIB)(%.o)}

INSTALL = install
RANLIB = ranlib

default: $(LIB) $(PROGS)

SUFFIXES = .f .o
.SUFFIXES: $(SUFFIXES)

.f.o:
	$(FC) ${FFLAGS} -c -o $@ $<

${LIB}:${LIBOBJ} 

distaz:call_distaz.f
	$(FC) ${FFLAGS} -o $@ $< ${LIB}

jday:call_jday.f
	$(FC) ${FFLAGS} -o $@ $< ${LIB}

install:
	${INSTALL} -c -m 0644 $(LIB) $(LIBDIR)
	${INSTALL} -c -m 0644 $(INC) $(INCDIR)
	${RANLIB} ${LIBDIR}/$(LIB)
	for i in ${PROGS}; do ${INSTALL} -c -m 0755 $$i $(BINDIR); done

clean:
	rm $(LIB) ${PROGS}
