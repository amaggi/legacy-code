#
# Makefile
#
# $Id: Makefile,v 1.6 2006/02/24 09:07:57 alessia Exp $
#
#
OPT=-g -Wall
CC=gcc
FC=gfortran
CFLAGS = $(OPT)
FFLAGS = $(OPT)

CODE=${HOME}/code

LIBDIR=${CODE}/lib
BINDIR=${CODE}/bin
INCDIR=${CODE}/include
LIBSRC=hilbert.c envelope.c cluster.f read_txt.f maxima.c interpolate.f sac_stats.c
INC=amutils.h
LIB=libamutils.a
PROGS = run_sac_stats

LIBOBJ = ${LIBSRC:%.c=$(LIB)(%.o)} ${LIBSRC:%.f=$(LIB)(%.o)}
INSTALL = install
RANLIB = ranlib

default: $(LIB)  $(PROGS)

SUFFIXES = .c .f .o
.SUFFIXES: $(SUFFIXES)

.c.o:
	$(CC) ${CFLAGS} -c -o $@ $< 
.f.o:
	$(FC) ${FFLAGS} -c -o $@ $< 

${LIB}:${LIBOBJ} ${INC}

run_sac_stats: run_sac_stats.c ${LIB}
	$(CC) $< -o $@ -L${LIBDIR} -I${INCDIR} ${LIB} -lsac -lm -lnumrec -lg2c -lgfortran

install:
	${INSTALL} -c -m 0644 $(LIB) $(LIBDIR)
	${INSTALL} -c -m 0644 $(INC) $(INCDIR)
	${RANLIB} ${LIBDIR}/$(LIB)
	for i in ${PROGS} ; do ${INSTALL} -c -m 0755 -s $$i ${BINDIR}; done

clean: 
	rm $(LIB) $(PROGS) 
