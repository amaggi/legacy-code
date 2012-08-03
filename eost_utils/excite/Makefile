FC=gfortran
FFLAGS =  -g -O5

LIBDIR=${CODE}/lib
INCDIR=${CODE}/include
BINDIR=${CODE}/bin
LIB=libexcite.a
LIBSRC = nrtype.f90 nrutil.f90 module_read_mode2.f90 module_spline.f90 legendre7.f90 distaz3.f90 modes_R2.f90 modes_S2.f90 modes_T2.f90 read_model2.f90


PROGS=excite_disp

LIBOBJ = ${LIBSRC:%.f90=$(LIB)(%.o)}

INSTALL = install
RANLIB = ranlib

default: $(LIB) $(PROGS)

SUFFIXES = .f90 .o
.SUFFIXES: $(SUFFIXES)

.f90.o:
	$(FC) ${FFLAGS} -c -o $@ $<

${LIB}:${LIBOBJ} 

excite_disp:excit_ter2_disp.f90
	$(FC) -o $@ $< ${LIB}

install:
	${INSTALL} -c -m 0644 $(LIB) $(LIBDIR)
	${RANLIB} ${LIBDIR}/$(LIB)
	for i in ${PROGS}; do ${INSTALL} -c -m 0755 $$i $(BINDIR); done

clean:
	rm $(LIB) *.mod
