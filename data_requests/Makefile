OPT= 
CC = gcc
FC = g77
CFLAGS = ${OPT}
FFLAGS = ${OPT}


BINDIR= ${HOME}/code/bin
INSTALL=install

PROGS= requestdata geoscope_request
default:${PROGS} 

LIBSRC= requestdata.c read_station_file.c read_station_time_file.c \
read_event_file.c distaz.f general_request_routines.c mt5_request.c \
station_request.c range_request.c read_STATIONS_file.c \
read_CMTSOLUTION_file.c

LIB=libauto.a
LIBOBJ = ${LIBSRC:%.f=${LIB}(%.o)} ${LIBSRC:%.c=${LIB}(%.o)}

${LIB}:${LIBOBJ} 


requestdata: requestdata.o ${LIB}
	$(FC) ${FFLAGS}  -o $@ requestdata.o ${LIB} ${LIBS}

geoscope_request: geoscope_request.o ${LIB}
	$(FC) ${FFLAGS}  -o $@ geoscope_request.o ${LIB} ${LIBS}

clean:
	rm -f  ${LIB} ${PROGS} *.o

install: 
	for i in ${PROGS}; do ${INSTALL} -c -m 0755 -s $$i $(BINDIR); done

