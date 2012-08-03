OPT= -g -O
CC = cc
CFLAGS = ${OPT}
FC =gfortran
FFLAGS = ${OPT} 
# SACLIBDIR would usually be inheritied from the environment...

BINDIR= ${CODE}/bin
LIBDIR= ${CODE}/lib


LIBSRC = earthsubs.f peripheral.f irisresp.f slresp.f char_int.f respc.c evresp.c mkhetero.f mkhomog.f plt_rtns.f fftl.f pitsa.f in_sac.f out_sac.f datetime.f 

LIB=libearthgram.a
LIBS=../grx/libgrx.a /usr/X11R6/liblibX11.a
#SOL_LIBS=${SOL_LIBDIR}/libgrx.a /usr/X/lib/libX11.a

LIBOBJ = ${LIBSRC:%.f=${LIB}(%.o)} ${LIBSRC:%.c=${LIB}(%.o)}

INSTALL=install

PROGS= earth gram

# any sane make should have these already but hey this is an example

default: ${PROGS}

SUFFIXES = .c .f .o

.SUFFIXES: ${SUFFIXES}

.c.o:
	$(CC) ${CFLAGS} -c -o $@ $<
.f.o:
	$(FC) ${FFLAGS} -c -o $@ $<

${LIB}:${LIBOBJ} 

${LIBC}:${LIBOBJC} 

earth:earth.o ${LIB} 
	${FC} ${FFLAGS} -o $@ $@.o ${LIB} 

gram:gram.o ${LIB}
	${FC} ${FFLAGS} -o $@ $@.o ${LIB} ${LIBS} ${SACLIB}

#gram_solaris:gram.o ${LIB}
#	${FC} ${FFLAGS} -o gram gram.o ${LIB} ${SOL_LIBS} ${SOL_SACLIB}

clean:
	rm -f *.o *.spg *.smb spag.log ${LIB} ${PROGS}

install: 
	for i in ${PROGS}; do ${INSTALL} -c -m 0755 -s $$i $(BINDIR); done

install_solaris: 
	for i in ${PROGS}; do ${INSTALL} -f $(SOL_BINDIR) -m 0755 -s $$i ; done
