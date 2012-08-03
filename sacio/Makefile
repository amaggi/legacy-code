#
# Makefile
#
# $Id: Makefile,v 1.2 2005/09/07 22:30:33 alessia Exp $
#
#
OPT=-g -Wall
CC=gcc
CFLAGS = $(OPT) -m32

CODE=${HOME}/code

LIBDIR=${CODE}/lib
BINDIR=${CODE}/bin
INCDIR=${CODE}/include
LIBSRC=sacio.c
LIB=libsac.a

LIBOBJ = ${LIBSRC:%.c=$(LIB)(%.o)}
INSTALL = install
RANLIB = ranlib

default: $(LIB) 

SUFFIXES = .c .o 
.SUFFIXES: $(SUFFIXES)

.c.o:
	$(CC) ${CFLAGS} -c -o $@ $< 

${LIB}:${LIBOBJ} sacio.h

install:
	${INSTALL} -c -m 0644 $(LIB) $(LIBDIR)
	${INSTALL} -c -m 0644 sacio.h $(INCDIR)
	${RANLIB} ${LIBDIR}/$(LIB)

clean: 
	rm $(LIB) 
