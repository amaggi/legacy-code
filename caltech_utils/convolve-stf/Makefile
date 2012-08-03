#FC=g77
FC=gfortran
CC=gcc
FFLAGS = -Wall
CFLAGS = -Wall

CODE=${HOME}/code
LIBDIR=${CODE}/lib
INCDIR=${CODE}/include
BINDIR=${CODE}/bin

PROGS=convolve_stf


INSTALL = install

default: $(PROGS)

SUFFIXES = .f .f90 .c .o
.SUFFIXES: $(SUFFIXES)

.f.o:
	$(FC) ${FFLAGS} -c -o $@ $< 

.f90.o:
	$(FC) ${FFLAGS} -c -o $@ $< 

.c.o:
	$(CC) ${CFLAGS} -c -o $@ $< 

convolve_stf:convolve_source_timefunction.f90
	$(FC) -o $@ $< -lm

install:
	for i in ${PROGS}; do ${INSTALL} -c -m 0755 $$i $(BINDIR); done

clean:
	rm $(PROGS) 
