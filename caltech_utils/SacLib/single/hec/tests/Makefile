CC = /usr/bin/gcc
CFLAG =  
F90 = f90
OPTF = -O2

LIB = -L/opt/seismo-util/lib/ -lRWFiles -lSacLib -lSacio -lm 
all : testcorrf90 testcorrc
testcorrf90 : testcorrf90.f90 
	$(F90)   $(OPTF) -o testcorrf90 testcorrf90.f90 \
	         $(LIB) -lRWFiles -lSacio -lm $(LIB)

testcorrc: testcorrc.c  
	$(CC) -o testcorrc testcorrc.c $(LIB) -I/opt/seismo-util/lib/includes
                

clean : 
	\rm -f $(OBJS) *~  *.o testcorrc  testcorrf90 *.ascii *.sgf

