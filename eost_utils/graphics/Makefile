OPT=-g -Wall
FC=g77
CFLAGS = $(OPT)

LIBDIR=${CODE}/lib
BINDIR=${CODE}/bin
INCDIR=${CODE}/include
LIBSRC=io.f math.f 
LIB=libgraphics.a

LIBPATH=-L${LIBDIR}

LIBOBJ = ${LIBSRC:%.f=$(LIB)(%.o)}
INSTALL = install
RANLIB = ranlib

PROGS = ray_density pourcentvs bin_ages sqrtT sta_ev rays_in_zone request_in_zone

default: $(LIB) ${PROGS}

SUFFIXES = .f .o
.SUFFIXES: $(SUFFIXES)

.f.o:
	$(FC) ${CFLAGS} -c -o $@ $< 

${LIB}:${LIBOBJ} 

ray_density:ray_density.f ${LIB}
	$(FC) ${FFLAGS} -o $@ $< ${LIB} $(LIBPATH) -lm

rays_in_zone:rays_in_zone.f ${LIB}
	$(FC) ${FFLAGS} -o $@ $< ${LIB} $(LIBPATH) -lm

request_in_zone:request_in_zone.f ${LIB}
	$(FC) ${FFLAGS} -o $@ $< ${LIB} $(LIBPATH) -lm

sta_ev:sta_ev.f ${LIB}
	$(FC) ${FFLAGS} -o $@ $< ${LIB} $(LIBPATH) -lm

pourcentvs:pourcentvs.f ${LIB}
	$(FC) ${FFLAGS} -o $@ $< ${LIB} $(LIBPATH) -lm


install:
	for i in ${PROGS}; do ${INSTALL} -m 0755 -s $$i $(BINDIR); done

clean: 
	rm $(LIB) ${PROGS} 

