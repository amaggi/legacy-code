OPT=-g -m32
FC=gfortran
FFLAGS=$(OPT)

CODE=${HOME}/code

BINDIR=${CODE}/bin
SHAREDIR=${CODE}/share

LIBDIR=${CODE}/lib

LIBSRC=libtau.f libsun.f emiasp91.f
LIB=libtau.a
LIBOBJ = ${LIBSRC:%.f=${LIB}(%.o)}

INSTALL = install

PROGS= remodl setbrn ttimes findit pick_index

default: ${PROGS} iaspei

SUFFIXES= .f .o
.SUFFIXES: ${SUFFIXES}

.f.o:
	${FC} ${FFLAGS} -c -o $@ $<

${LIB}:${LIBOBJ}

remodl:remodl.o limits.inc ${LIB}
	${FC} ${FFLAGS} -o $@ $@.o ${LIB} 

setbrn:setbrn.o ${LIB}
	${FC} ${FFLAGS} -o $@ $@.o  ${LIB}

ttimes:ttimes.o ttlim.inc ${LIB}
	${FC} ${FFLAGS} -o $@ $@.o ${LIB} 

findit:findit.o 
	${FC} ${FFLAGS} -o $@ $@.o 

pick_index:pick_index.o 
	${FC} ${FFLAGS} -o $@ $@.o -L${LIBDIR} -lsac 

iaspei: remodl setbrn
	./remodl
	./setbrn

clean:
	rm -f *.o ${LIB} ${PROGS} iasp91.hed iasp91.tbl

install:
	for i in ${PROGS}; do ${INSTALL} -m 0755 -s $$i $(BINDIR); done
	for i in SEED_prep.sh data_prep.sh iasp_prep.sh run_rdseed ; do ${INSTALL} -m 0755 $$i $(BINDIR); done
	${INSTALL} -c -m 0644 iasp91.hed iasp91.tbl resp_wwlpbn $(SHAREDIR)
	${INSTALL} -c -m 0644 $(LIB) $(LIBDIR)
	ranlib  $(LIBDIR)/${LIB}

	
