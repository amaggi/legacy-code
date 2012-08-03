CC = /usr/bin/gcc
F90 = f90
OPTF = -O2 
OPTC = 

testsaclib : testsaclib.o ../../libSacLib.a
	$(F90) $(OPTF) $(PORT) -o testsaclib  testsaclib.o \
        -L/opt/seismo-util/lib/  -lSacLib -lRWFiles -lSacio -lm

testsaclib.o : testsaclib.f90
	$(F90) $(OPTF) $(PORT) -c testsaclib.f90
clean: 
	\rm -f *.o testsaclib core *~ 

