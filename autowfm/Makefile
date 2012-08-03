FFLAGS = -g -Wall
CFLAGS = -g -Wall -O
FC = gfortran
CC = gcc
SAC_LIB = ${HOME}/code/lib/libsac.a
BINDIR = ${HOME}/bin
LIBS = ${SAC_LIB} -lm

PROG=chheadz1.sac read_resp sortdatavgz2.sac preselec.tz3 autosac2wfm \
     dat2nomcoo 3SMACGCavCrust.spremQcst maklst.mod saito2002 \
     initsaito2 autowfm.lovray convdp.lrlog10 foray95 folov94 \
     choose_resp estampsynt gktv.morelobes dgk12.morelobes estlin2.lrlog10 \
     tv.cross12.lrlog10.phase autowfm.lovray.shortperiod \
     autowfm.lovray.only-50-70 radpattern des2intomo regiostack \
     pourcentvs syntsmac forcbrut GROUPintomodesVs errstatintomodesVs \
     deletenoderand

LIBSRC = daymo.f lpyr.f  lecdgk12.phase.lrlog10.f tvinvma22b.f \
         sordgk3.lrlog10.f ma22b.f stack.c
LIB=libautowfm.a

LIBOBJ = ${LIBSRC:%.f=${LIB}(%.o)} ${LIBSRC:%.c=${LIB}(%.o)}

default: ${LIB} ${PROG}

${LIB}:${LIBOBJ}

.f.o:
	$(FC) ${FFLAGS} -c -o $@ $< 

.c.o:
	$(CC) ${CFLAGS} -c -o $@ $< 

chheadz1.sac:chheadz1.sac.f 
	${FC} ${FFLAGS} -o $@ $< ${SAC_LIB}

read_resp:read_resp.f 
	${FC} ${FFLAGS} -o $@ $< 

sortdatavgz2.sac:sortdatavgz2.sac.f 
	${FC} ${FFLAGS} -o $@ $< ${SAC_LIB}

preselec.tz3:preselec.tz3.c 
	${CC} ${CFLAGS} -o $@ $< ${LIBS}

autosac2wfm:autosac2wfm.f 
	${FC} ${FFLAGS} -o $@ $< ${SAC_LIB} ${LIB}

dat2nomcoo:dat2nomcoo.f 
	${FC} ${FFLAGS} -o $@ $< 

3SMACGCavCrust.spremQcst:3SMACGCavCrust.spremQcst.c 
	${CC} ${CFLAGS} -o $@ $< -lm

maklst.mod:maklst.mod.c 
	${CC} ${CFLAGS} -o $@ $< -lm

saito2002:saito2002.f 
	${FC} ${FFLAGS} -o $@ $< 

initsaito2:initsaito2.f 
	${FC} ${FFLAGS} -o $@ $< 

autowfm.lovray:autowfm.lovray.sylvana.linux.c
	${CC} ${CFLAGS} -o $@ $<  -lm

autowfm.lovray.shortperiod:autowfm.lovray.sylvana.linux.shortperiod.c
	${CC} ${CFLAGS} -o $@ $<  -lm

autowfm.lovray.only-50-70:autowfm.lovray.sylvana.linux.only-50-70.c
	${CC} ${CFLAGS} -o $@ $<  -lm

convdp.lrlog10:convdp.lrlog10.f 
	${FC} ${FFLAGS} -o $@ $< 

foray95:foray95.f 
	${FC} ${FFLAGS} -o $@ $< 

folov94:folov94.f 
	${FC} ${FFLAGS} -o $@ $< 

choose_resp:choose_resp.f 
	${FC} ${FFLAGS} -o $@ $< 

estampsynt:estampsynt.f 
	${FC} ${FFLAGS} -o $@ $< 

gktv.morelobes:gktv.morelobes.f 
	${FC} ${FFLAGS} -o $@ $< 

dgk12.morelobes:dgk12.morelobes.f 
	${FC} ${FFLAGS} -o $@ $< 

estlin2.lrlog10:estlin2.lrlog10.f 
	${FC} ${FFLAGS} -o $@ $< 

tv.cross12.lrlog10.phase:tv.cross.lr.f ${LIB}
	${FC} ${FFLAGS} -o $@ $< ${LIB}

radpattern:rayrad94e.f
	${FC} ${FFLAGS} -o $@ $< 

des2intomo:des2intomoGLOBV3coord.f
	${FC} ${FFLAGS} -o $@ $< 

regiostack:regiostack6.f ${LIB}
	${FC} ${FFLAGS} -o $@ $<  ${LIB}

pourcentvs:pourcenT.VS.aus.f 
	${FC} ${FFLAGS} -o $@ $<  

syntsmac:IntSYNTSMAC.f 
	${FC} ${FFLAGS} -o $@ $<  

forcbrut:forcbrut2V3.c  nnsphereric.o
	${CC} ${CFLAGS} -o $@ $<  nnsphereric.o -Lqhull_src -lqhull -lm -lg2c 

deltenoderand:deltenoderand.c  
	${CC} ${CFLAGS} -o $@ $<  -lm 


GROUPintomodesVs:GROUPintomodesVSV3.f
	${FC} ${FFLAGS} -o $@ $<  

errstatintomodesVs:errstatintomodesVsV3.f
	${FC} ${FFLAGS} -o $@ $<  


clean:
	rm -f ${PROG} *.o

install:
	for i in ${PROG} ; do install  -m 0755 -s $$i $(BINDIR) ; done
