OPT= -g -O
CC = gcc
CFLAGS = ${OPT}
FC = gfortran
FFLAGS= ${OPT} 

BINDIR= ${CODE}/bin

INSTALL=install

PROGS= displ dispr  

default: ${PROGS}

SUFFIXES= .c .f .o

.SUFFIXES= ${SUFFIXES}

.c.o:
	${CC} ${CFLAGS} -c -o $@ $<
.f.o:
	${FC} ${FFLAGS} -c -o $@ $<

displ:displ.o 
	${FC} ${FFLAGS} -o $@ $@.o 

dispr:dispr.o 
	${FC} ${FFLAGS} -o $@ $@.o 

plot_derivatives:plot_derivatives.o
	${CC} -o $@ $@.o 

clean:
	rm -f *.o  ${PROGS}

install:
	for i in ${PROGS}; do ${INSTALL} -c -m 0755 -s $$i $(BINDIR); done

