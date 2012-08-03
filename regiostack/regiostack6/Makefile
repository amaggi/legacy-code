F90 = gfortran
PROGS = xregiostack6 xregiostack6dt
BIN = ${HOME}/code/bin

all:${PROGS}

xregiostack6: regiostack6_ompV3d_alessia.f
	${F90} -O2 regiostack6_ompV3d_alessia.f -o xregiostack6

xregiostack6dt: regiostack6_ompV3d_alessia_delaytime.f
	${F90} -O2 regiostack6_ompV3d_alessia_delaytime.f -o xregiostack6dt

install: ${PROGS}
	for i in ${PROGS}; do install -m 0755 -s $$i $(BIN) ; done

clean: 
	rm ${PROGS}
# for compilation under OMP on IDRIS machine zahir:
#F77 = xlf_r -q64
#all:regiostack6_ompV3d
#DEST = /homegpfs/rech/ing/ring884/bin/
#
#regiostack6_ompV3d : regiostack6_ompV3d_alessia.o 
#	$(F77) -qsmp=omp -O3  -qnosave  -o $(DEST)regiostack6_ompV3d regiostack6_ompV3d_alessia.o
#	rm regiostack6_ompV3d_alessia.o
