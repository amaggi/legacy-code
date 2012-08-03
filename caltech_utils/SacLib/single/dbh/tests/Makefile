CC = /usr/bin/gcc
F90 = f90
OPTF = -O2 
OPTC = 
LIB =  -L/opt/seismo-util/lib/  -lSacLib -lRWFiles -lSacio  -lm
# options to link with c functions

testsaclib : testsaclib.o 
	$(F90) $(OPTF) -o testsaclib  testsaclib.o $(LIB)

testsaclib.o : testsaclib.f90
	$(F90) $(OPTF) $(PORT) -c testsaclib.f90
clean: 
	\rm -f *.o testsaclib core *~ *.sgf temp.ps sac.tmp sac.out 

