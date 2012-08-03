#
# Makefile
#
# $Id: Makefile,v 1.5 2002/07/10 16:24:16 maggi Exp $
#
#
OPT=-g -Wall
CC=gcc
CFLAGS = $(OPT)

LIBDIR=${HOME}/code/lib
INCDIR=${HOME}/code/include
LIBSRC= geom.c globals.c io.c mem.c poly.c qhull.c set.c
LIB=libqhull.a

LIBOBJ = ${LIBSRC:%.c=$(LIB)(%.o)}
INSTALL = install
RANLIB = ranlib

default: $(LIB) 

SUFFIXES = .c .o
.SUFFIXES: $(SUFFIXES)

.c.o:
	$(CC) ${CFLAGS} -c -o $@ $< 

${LIB}:${LIBOBJ} 

install:
	${INSTALL} -c -m 0644 $(LIB) $(LIBDIR)
	${RANLIB} ${LIBDIR}/$(LIB)

clean: 
	rm $(LIB) *.o
