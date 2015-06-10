OPT=-g -Wall
FC=gfortran
FFLAGS = $(OPT)

LIBDIR=${HOME}/lib
BINDIR=${HOME}/bin
INCDIR=${HOME}/include
LIBSRC=synt_io.f synt_math.f choices.f Tcammarano.f lithcooling.f
LIBSRCV3=synt_ioV3.f synt_math.f choices.f Tcammarano.f lithcooling.f choices_dt.f kaminski_ribe.f
LIB=libsyntmodel.a
LIBV3=libsyntmodelV3.a

LIBPATH=-L${LIBDIR}

LIBOBJ = ${LIBSRC:%.f=$(LIB)(%.o)}
LIBOBJV3 = ${LIBSRCV3:%.f=$(LIBV3)(%.o)}
INSTALL = install
RANLIB = ranlib

PROGS = IntSyntV3 IntSyntDT intomocov
#PROGS = IntSYNTPREM-choice IntSYNTPREM-choice-clustered IntSYNTSMAC \
#IntSYNTSMAC-clustered IntSYNTSMAC-clusteredV3 IntSYNTPREM-choice-clusteredV3


default: $(LIB) $(LIBV3) ${PROGS}

SUFFIXES = .f .o
.SUFFIXES: $(SUFFIXES)

.f.o:
	$(FC) ${FFLAGS} -c -o $@ $< 

${LIB}:${LIBOBJ} 
${LIBV3}:${LIBOBJV3} 

IntSYNTPREM-choice:IntSYNTPREM-choice.f ${LIB}
	$(FC) ${FFLAGS} -o $@ $< ${LIB} $(LIBPATH) -lutil -lm

IntSYNTPREM-choice-clustered:IntSYNTPREM-choice-clustered.f ${LIB}
	$(FC) ${FFLAGS} -o $@ $< ${LIB} $(LIBPATH) -lutil -lm

IntSYNTPREM-choice-clusteredV3:IntSYNTPREM-choice-clusteredV3.f ${LIBV3}
	$(FC) ${FFLAGS} -o $@ $< ${LIBV3} $(LIBPATH) -lutil -lm

IntSyntV3:IntSyntV3.f ${LIBV3}
	$(FC) ${FFLAGS} -o $@ $< ${LIBV3} $(LIBPATH) -lutil -lnumrec -lm

IntSyntDT:IntSyntDT.f ${LIBV3}
	$(FC) ${FFLAGS} -o $@ $< ${LIBV3} $(LIBPATH) -lutil -lm

intomocov:intomo_coverage.f ${LIBV3}
	$(FC) ${FFLAGS} -o $@ $< ${LIBV3} $(LIBPATH) -lutil -lm


IntSYNTSMAC:IntSYNTSMAC.f ${LIB}
	$(FC) ${FFLAGS} -o $@ $< ${LIB} $(LIBPATH) -lutil -lm

IntSYNTSMAC-clustered:IntSYNTSMAC-clustered.f ${LIB}
	$(FC) ${FFLAGS} -o $@ $< ${LIB} $(LIBPATH) -lutil -lm

IntSYNTSMAC-clusteredV3: IntSYNTSMAC-clusteredV3.f $(LIBV3)
	$(FC) ${FFLAGS} -o $@ $< ${LIBV3} $(LIBPATH) -lutil -lm


install:
	for i in ${PROGS}; do ${INSTALL} -m 0755 -s $$i $(BINDIR); done

clean: 
	rm $(LIB) ${PROGS}  ${LIBV3}

cleanrun:
	rm *.xyz des.S* intomoPREM gc.cross.xy
