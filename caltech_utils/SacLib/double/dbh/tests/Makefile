CC = /usr/bin/gcc
F90 = f90
OPTF = -O2 
OPTC = 

testsaclib : testsaclib.o 
	$(F90) $(OPTF) $(PORT) -g -o testsaclib  testsaclib.o \
        -L/opt/seismo-util/lib/ -lDRWFiles -lDSacLib -lDSacio -lm 

testsaclib.o : testsaclib.f90
	$(F90) $(OPTF) $(PORT) -g -c testsaclib.f90
clean: 
	\rm -f *.o testsaclib core *~ *.sgf

