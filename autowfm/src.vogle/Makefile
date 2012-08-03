# MAKEFILE POUR COMPILATION SOUS LINUX
VOGLE   =      -L. -lcore2vogle libvogle.a -L/usr/X11R6/lib -lX11 -lefence

OPTIONS		= -Wall
DEST	= ${HOME}/bin
BNDCHK      =     -C
LIBSTD      =     /usr/local/sac/lib/sac.a  -lm
F77     = g77
CORE = $(VOGLE)

all: synt5.lvogl cross12.lausol dessin.lvogl

%: %.f
	$(F77) $(OPTIONS) $*.f $(LIBSTD) -o $(DEST)/$@

%: %.c
	cc $(OPTIONS) $*.c $(LIBSTD) -o $(DEST)/$@

libcore2vogle.a: suncore2vogle.f
	$(F77) -c suncore2vogle.f
	ar cr libcore2vogle.a suncore2vogle.o
	ranlib libcore2vogle.a

clean:
	make clean1
	make clean2
	make cleanps

clean1:
	yes|rm error dgk*.it0 *.it[1-9]* *.it1[0-9]*
clean2:
	yes|rm att.init.ray disp.init.ray dpa.init.ray tv.init \
	gk.init.ray.[1-9] sour.init.ray.[1-9] tf.init.ray.[1-9] \
	pick.ray.[1-9]
	yes|rm att.init.lov disp.init.lov dpa.init.lov tv.init \
	gk.init.lov.[1-9] sour.init.lov.[1-9] tf.init.lov.[1-9] \
	pick.lov.[1-9]

cleanps:
	yes|rm *.postscript

cleanlr:
	make cleanlr1
	make cleanlr2
	make cleanps

cleanlr1:
	yes|rm error? *.it* des.so.? des.so.??

cleanlr2:
	yes|rm *.init*

cleandgk:
	yes|rm dgk.lov* dgk.ray*

clean.o: 
	\rm *.o

etudedephasage: etudedephasage.f
	$(F77) $(OPTIONS) etudedephasage.f -o $(DEST)/etudedephasage

synt5.lvogl: synt5.lr.f  libcore2vogle.a
	$(F77) synt5.lr.f \
	-o $(DEST)/synt5.lvogl $(VOGLE)

cross13.lr: cross13.lr.o des.cross.o \
	suncore2vogle.o lcore.dummy.o
	$(F77) $(OPTIONS)  cross13.lr.o des.cross.o \
	suncore2vogle.o coredum.o -o $(DEST)/cross13.lr $(CORE)
	$(F77) $(OPTIONS)  cross13.lr.o initdes.ps.o \
	suncore2vogle.o lcore.dummy.o coredum.o \
	-o $(DEST)/cross13.lr.nodess

cross13.lr.o:cross13.lr.F convsuncore.h
	$(F77) $(OPTIONS)  -c cross13.lr.F -o cross13.lr.o

cross12.lvogl: cross12.lr.o lcore.dummy.o libcore2vogle.a
	$(F77) $(OPTIONS)  cross12.lr.o -o $(DEST)/cross12.lvogl $(VOGLE)
	$(F77) $(OPTIONS)  cross12.lr.o lcore.dummy.o \
	-o $(DEST)/cross12.lvogl.nodess

cross12.lr.o:cross12.lr.f 
	$(F77) $(OPTIONS)  -c cross12.lr.f -o cross12.lr.o

cross12.ausol: cross12.ausol.o lcore.dummy.o libcore2vogle.a 
	$(F77) $(OPTIONS)  cross12.ausol.o  -o $(DEST)/cross12.ausol $(VOGLE)
	$(F77) $(OPTIONS)  cross12.ausol.o lcore.dummy.o \
	-o $(DEST)/cross12.ausol.nodess

cross12.lausol: cross12.lausol.o lcore.dummy.o libcore2vogle.a 
	$(F77) $(OPTIONS)  cross12.lausol.o  -o $(DEST)/cross12.lausol $(VOGLE)
	$(F77) $(OPTIONS)  cross12.lausol.o lcore.dummy.o \
	-o $(DEST)/cross12.lausol.nodess

cross12.lausol.o:cross12.lausol.f 
	$(F77) $(OPTIONS)  -c cross12.lausol.f -o cross12.lausol.o

cross12.ausol2: cross12.ausol2.o lcore.dummy.o libcore2vogle.a 
	$(F77) $(OPTIONS)  cross12.ausol2.o  -o $(DEST)/cross12.ausol2 $(VOGLE)
	$(F77) $(OPTIONS)  cross12.ausol2.o lcore.dummy.o \
	-o $(DEST)/cross12.ausol2.nodess

cross12.ausol2.o:cross12.ausol2.f 
	$(F77) $(OPTIONS)  -c cross12.ausol2.f -o cross12.ausol2.o

lcore.dummy.o: lcore.dummy.f
	$(F77) $(OPTIONS)  -c lcore.dummy.f -o lcore.dummy.o

cross12.ausol.o:cross12.ausol.f 
	$(F77) $(OPTIONS)  -c cross12.ausol.f -o cross12.ausol.o

cross12.ausol2: cross12.ausol2.o lcore.dummy.o libcore2vogle.a 
	$(F77) $(OPTIONS)  cross12.ausol2.o  -o $(DEST)/cross12.ausol2 $(VOGLE)
	$(F77) $(OPTIONS)  cross12.ausol2.o lcore.dummy.o \
	-o $(DEST)/cross12.ausol2.nodess

cross12.ausol2.o:cross12.ausol2.f 
	$(F77) $(OPTIONS)  -c cross12.ausol2.f -o cross12.ausol2.o

lcore.dummy.o: lcore.dummy.f
	$(F77) $(OPTIONS)  -c lcore.dummy.f -o lcore.dummy.o

desresol.lvogl: desresol.lvogl.f libcore2vogle.a
	$(F77) $(OPTIONS)  desresol.lvogl.f -o $(DEST)/desresol.lvogl $(VOGLE)

dessin.lvogl: dessin.f libcore2vogle.a
	$(F77) $(OPTIONS)  dessin.f -o $(DEST)/dessin.lvogl $(VOGLE)

dessbar.lvogl: dessbar.f libcore2vogle.a
	$(F77) $(OPTIONS)  dessbar.f -o $(DEST)/dessbar $(VOGLE)

gradue.o: gradue.f libcore2vogle.a
	$(F77) $(OPTIONS) -c gradue.f -o gradue.o

desdepl.lvogl: desdepl.lvogl.o gradue.o
	$(F77) $(OPTIONS)  desdepl.lvogl.o \
	gradue.o -o $(DEST)/desdepl.lvogl $(VOGLE)

desdepl.lvogl.o: desdepl.lvogl.f libcore2vogle.a
	$(F77) $(OPTIONS)  -c desdepl.lvogl.f -o desdepl.lvogl.o
