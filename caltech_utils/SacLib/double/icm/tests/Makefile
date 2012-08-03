CC = /usr/bin/gcc
F90 = f90
OPTF = -O2 
OPTC = 
# options to link with c functions
PORT =# -YEXT_NAMES=LCS  -YEXT_SFX=_ # for Absoft f90 compiler
#PORT = -f -s -B108 -N90  # for Absoft f77 compiler

testsaclib : testsaclib.o ../../libDSacLib.a
	$(F90) $(OPTF) $(PORT) -o testsaclib  testsaclib.o \
        -L/opt/seismo-util/lib/  -lDSacio  -lm ../../libDSacLib.a

testsaclib.o : testsaclib.f90
	$(F90) $(OPTF) $(PORT) -c testsaclib.f90
clean: 
	\rm -f *.o testsaclib core *~ 

