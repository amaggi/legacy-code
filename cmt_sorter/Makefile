OPT= -g -Wall
CC = gcc
FC = gfortran
CFLAGS = ${OPT}
FFLAGS = ${OPT}

CODE=${HOME}/code

BINDIR= ${CODE}/bin
INSTALL=install
LIBDIR=${CODE}/lib 
LIBS=-L${LIBDIR} -lutil -lm 
#LIBS=-L${LIBDIR} -L${INTEL_LIBDIR} -lutil -lm 

PROGS= cmt_sorter
default:${PROGS} 

LIBSRC=parser.c delta_angle.c fixed_len_atonumber.c print_help.c file_io.c criteria.c  julian.c

LIB=libcmt_sorter.a
LIBOBJ = ${LIBSRC:%.f=${LIB}(%.o)} ${LIBSRC:%.c=${LIB}(%.o)}

${LIB}:${LIBOBJ} header.h


cmt_sorter: ${LIB} cmt_sorter.o  
	$(CC) ${FFLAGS}  -o $@ cmt_sorter.o ${LIB} ${LIBS} 

clean:
	rm -f  ${LIB} ${PROGS} *.o

install: 
	for i in ${PROGS}; do ${INSTALL} -c -m 0755 -s $$i $(BINDIR); done

